clear %%THIS USES CG
%% Constants CHECKED
Re=10000000; %20 and 400 
H=1; 
U=1; %inputs
eps=10e-5; %error for CG solver
dt_guess=0.01;
tf=100;

Ta=1500;
Tb=1000;
A=1000;

T_init=300;

D=10^-4;

L = 15*H;
W = 5*H;

dx=H/20; dy=H/20; lambda=1/dx^2; 
nx=L/dx; 
ny=W/dy; %number of cells in x,y

Niter=ceil(200/dt_guess); %number of iterations

CFL=0.8;
%% Initialize Arrays

p=zeros(nx,ny);

u=zeros(nx+1,ny); %u array
u(1,:)=U;
ustar=u;
ug=zeros(nx+2,ny+2); %ghost cell array %adding 1 col to right, 1 above and below
ug(1,2:end-1)=U; %set left boundary to U=1 
ug(1:end-1,2:end-1) = u;%update ug with inner u
ug(end,:) = ug(end-2,:); %to impose neumann on RHS

v=zeros(nx,ny+1); %v array
vstar=v;
vg=zeros(nx+2,ny+1); %ghost cell array %adding 1 col to left and right, nothing above/below
vg(2:end-1,:)=v;
vg(end,:)=vg(end-1,:); %update vg with inner v

H_u_new=zeros(nx,ny);
H_v_new=zeros(nx,ny-1);
RHS1u=zeros(nx+1,ny);
RHS1v=zeros(nx,ny+1);

% delete these once solver works
y=linspace(-2.5*H,2.5*H,ny);
phi_inlet=exp(-(y/(0.5*H)).^2);
H_phi_old=zeros(nx,ny);
phi=zeros(nx,ny);
phig=zeros(nx+3,ny+4); %augmented 2 on top and bottom, 2 on right, 1 on left
phi(1,:)=phi_inlet;
phig(2,3:end-2)=phi(1,:);
phig(1,3:end-2)=phi(1,:);
phi_x=zeros(nx+1,ny);
phi_x(1,:)=phi_inlet;
phi_y=zeros(nx,ny+1); %not caring about vertical phi at top and bottom since that gets zerod out


%T and xphi should have same ghost configs, since dirichlet at left,
%neumann elsewhere
H_T_old=zeros(nx,ny);
T = ones(nx,ny);
T = T_init*T;
Tg=ones(nx+3,ny+4); %augmented 2 on top and bottom, 2 on right, 1 on left
Tg = T_init*Tg;

T_x=zeros(nx+1,ny);
T_x(1,:)=T_init;
T_y=zeros(nx,ny+1); %not caring about vertical phi at top and bottom since that gets zerod out

H_xphi_old=zeros(nx,ny);
xphi = ones(nx,ny);
xphig=ones(nx+3,ny+4); %augmented 2 on top and bottom, 2 on right, 1 on left
xphi_x=zeros(nx+1,ny);
xphi_x(1,:)=1;
xphi_y=zeros(nx,ny+1); %not caring about vertical phi at top and bottom since that gets zerod out


%% Square geometry CHECKED

L_s=H;
W_s=H;
nx_s=L_s/dx;
ny_s=W_s/dy;

l_edge_s=5*H/dx;
r_edge_s=l_edge_s+nx_s+1;
b_edge_s=2*H/dy;
t_edge_s=b_edge_s+ny_s+1;

l_edge_f=r_edge_s;
r_edge_f=l_edge_f+nx_s-1;
b_edge_f=b_edge_s+1;
t_edge_f=t_edge_s-1;

%impose velocity = 0 inside and on boundary of square for no slip and 0
%pressure inside
u(l_edge_s+1:r_edge_s,b_edge_s+1:t_edge_s-1)=0;
v(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s)=0;
p(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;

%impose ignition flame condition T = 2500 in a H x H region behind square
T(l_edge_f:r_edge_f,b_edge_f:t_edge_f)=2500;

%impose T at square = 0
T(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0; %unnecessary?
Tg(2:end-2,3:end-2)=T;

%impose xphi at square = 0
xphi(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0; %unnecessary?
xphig(2:end-2,3:end-2)=xphi;

%laplace multiplier to make sure pressure in the square evaluates to 0
pmultiplier=ones(nx,ny);
pmultiplier(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;


%%  Laplacian coefficient array CHECKED
lapcoeff=4*ones(nx,ny);
lapcoeff(:,1)=3; lapcoeff(:,end)=3; lapcoeff(1,:)=3; lapcoeff(end,:)=3;
lapcoeff(1,1)=2; lapcoeff(1,end)=2; lapcoeff(end,1)=2; lapcoeff(end,end)=2;

lapcoeff(l_edge_s,b_edge_s+1:t_edge_s-1)=3; %left side
lapcoeff(r_edge_s,b_edge_s+1:t_edge_s-1)=3; %right side
lapcoeff(l_edge_s+1:r_edge_s-1,t_edge_s)=3; %top side
lapcoeff(l_edge_s+1:r_edge_s-1,b_edge_s)=3; %bottom side


%% Iterator CHECKED
%%% REMEMBER TO PLOT FLIP(U') or rot90(U). Default u is i-1/2,j. Default v
%%% is i, j-1/2
pg_GS=zeros(nx+2,ny+2); %adding 1 column on left/right, 1 row on top/bottom
lap_GS=zeros(nx,ny);
c1=clock;
t=zeros(Niter*5,1);
C_l=zeros(Niter*5,1);
C_d=zeros(Niter*5,1);
t(1)=0; %change to 0

%%
for iter=1:2*Niter
    %tic
    
    dt_constrained = CFL*dx/(max(max(abs(u))) + max(max(abs(v))));
    %dt=min(dt_constrained,dt_guess);
    dt=dt_constrained;
    t(iter+1) = t(iter) + dt;
    t(iter+1)
%%%%%%%%%%%%%%%%%%%%%%%%%% SQUARE NEUMANN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %bottom neumann
    T(l_edge_s+1:r_edge_s-1,b_edge_s+1)=T(l_edge_s+1:r_edge_s-1,b_edge_s);
    T(l_edge_s+1:r_edge_s-1,b_edge_s+2)=T(l_edge_s+1:r_edge_s-1,b_edge_s);

    %top neumann
    T(l_edge_s+1:r_edge_s-1,t_edge_s-1)=T(l_edge_s+1:r_edge_s-1,t_edge_s);
    T(l_edge_s+1:r_edge_s-1,t_edge_s-2)=T(l_edge_s+1:r_edge_s-1,t_edge_s);
    
    %left neumann
    T(l_edge_s+1,b_edge_s+1:t_edge_s-1)=T(l_edge_s,b_edge_s+1:t_edge_s-1); 
    T(l_edge_s+2,b_edge_s+1:t_edge_s-1)=T(l_edge_s,b_edge_s+1:t_edge_s-1); 

    %right neumann
    T(r_edge_s-1,b_edge_s+1:t_edge_s-1)=T(r_edge_s,b_edge_s+1:t_edge_s-1);
    T(r_edge_s-2,b_edge_s+1:t_edge_s-1)=T(r_edge_s,b_edge_s+1:t_edge_s-1);
%%%%%%%%%%%%%%%%%%%%%%%%%% SQUARE NEUMANN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%% wxphi and wT CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wdot=phi.*((T/Tb).^3).*exp(-Ta./T);
    wxphi=-wdot;
    wT=A*wdot;
%%%%%%%%%%%%%%%% wxphi and wT CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%%%%%%%%%%%%%%%% PHI CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %implement neumann boundary
    phig(end,:)=phig(end-2,:); %neumann on right
    phig(end-1,:)=phig(end-2,:); %neumann on right
    
    phig(:,end)=phig(:,end-2); %neumann on top
    phig(:,end-1)=phig(:,end-2); %neumann on top
    
    phig(:,2)=phig(:,3); %neumann on bottom
    phig(:,1)=phig(:,3); %neumann on bottom
    
    %first work out phi bar for x at each cell
    for i = 2:nx+1 %note this is indexing x 1 in from left since we know it at i = 1 = inlet
        for j = 1:ny
            if u(i,j)>=0
                phi_x(i,j)=(-1/6)*phig(i-1,j+2) + (5/6)*phig(i,j+2) + (2/6)*phig(i+1,j+2);
            else
                phi_x(i,j)=(2/6)*phig(i,j+2) + (5/6)*phig(i+1,j+2) + (-1/6)*phig(i+2,j+2);
            end
        end
    end
    
    %then work out phi bar for y at each cell
    for i = 1:nx %note this is indexing x 1 in from left
        for j = 1:ny+1
            if v(i,j)>=0
                phi_y(i,j)=(-1/6)*phig(i+1,j) + (5/6)*phig(i+1,j+1) + (2/6)*phig(i+1,j+2);
            else
                phi_y(i,j)=(2/6)*phig(i+1,j+1) + (5/6)*phig(i+1,j+2) + (-1/6)*phig(i+1,j+3);
            end
        end
    end
    %then put together with u to get convective term
    
    %set phi to 0 in square
    
    fu=(u.*phi_x);
    fv=(v.*phi_y);
    H_phi_new= (fu(2:end,:)-fu(1:end-1,:))/dx + (fv(:,2:end)-fv(:,1:end-1))/dy;
    
    CONV = 1.5*H_phi_new -0.5*H_phi_old;
    
    %then calc diffusion terms
    d2phidx2=(phig(1:end-3,3:end-2) - 2*phig(2:end-2,3:end-2) + phig(3:end-1,3:end-2))/dx^2;
    d2phidy2=(phig(2:end-2,2:end-3) - 2*phig(2:end-2,3:end-2) + phig(2:end-2,4:end-1))/dy^2;
    
    
    %then do ADI to calc phi n+1
    
    phiRHS=-dt*CONV + dt*D*(d2phidx2 + d2phidy2); %is this negative?
    phiRHS(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;
    
    %dphi**
    Cphix = D*dt/(2*dx^2);
    main_temp = (1+2*Cphix)*ones(nx-1,1);
    lower_temp = -Cphix*ones(nx-1,1);
    %lower_temp(end)=(-2*Cphix); %this is for neumann condition
    main_temp(end)=(1+Cphix); %this is for neumann condition
    upper_temp = -Cphix*ones(nx-1,1);
    dphi_ss=zeros(nx,ny);
    for j = 1:ny
        dphi_ss(2:end,j) = ADRtridiag(main_temp,lower_temp,upper_temp,phiRHS(2:end,j));
    end
    
    %dphi THIS IS CORRECT SINCE dphi zero above and below
    Cphiy = D*dt/(2*dy^2);
    main_temp = (1+2*Cphiy)*ones(ny,1);
    lower_temp = -Cphiy*ones(ny,1);
    upper_temp = -Cphiy*ones(ny,1);
    dphi=zeros(nx,ny);
    for i = 2:nx
        dphi(i,:) = ADRtridiag(main_temp,lower_temp,upper_temp,dphi_ss(i,:));
    end
    
    phi = phi + dphi;
    phi(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0; %unnecessary?
    phig(2:end-2,3:end-2)=phi;
    
    %impose phi at square = 0
    
    %WHERE WOULD I IMPLEMENT BC?

    H_phi_old=H_phi_new;
    
%%%%%%%%%%%%%%%% PHI CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% T CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %implement neumann boundary
    Tg(end,:)=Tg(end-2,:); %neumann on right
    Tg(end-1,:)=Tg(end-2,:); %neumann on right
    
    Tg(:,end)=Tg(:,end-2); %neumann on top
    Tg(:,end-1)=Tg(:,end-2); %neumann on top
    
    Tg(:,2)=Tg(:,3); %neumann on bottom
    Tg(:,1)=Tg(:,3); %neumann on bottom
    
    %first work out T bar for x at each cell
    for i = 2:nx+1 %note this is indexing x 1 in from left since we know it at i = 1 = inlet
        for j = 1:ny
            if u(i,j)>=0
                T_x(i,j)=(-1/6)*Tg(i-1,j+2) + (5/6)*Tg(i,j+2) + (2/6)*Tg(i+1,j+2);
            else
                T_x(i,j)=(2/6)*Tg(i,j+2) + (5/6)*Tg(i+1,j+2) + (-1/6)*Tg(i+2,j+2);
            end
        end
    end
    
    %then work out T bar for y at each cell
    for i = 1:nx %note this is indexing x 1 in from left
        for j = 1:ny+1
            if v(i,j)>=0
                T_y(i,j)=(-1/6)*Tg(i+1,j) + (5/6)*Tg(i+1,j+1) + (2/6)*Tg(i+1,j+2);
            else
                T_y(i,j)=(2/6)*Tg(i+1,j+1) + (5/6)*Tg(i+1,j+2) + (-1/6)*Tg(i+1,j+3);
            end
        end
    end
    %then put together with u to get convective term
    
    %set phi to 0 in square
    
    T_fu=(u.*T_x);
    T_fv=(v.*T_y);
    H_T_new= (T_fu(2:end,:)-T_fu(1:end-1,:))/dx + (T_fv(:,2:end)-T_fv(:,1:end-1))/dy;
    
    T_CONV = 1.5*H_T_new -0.5*H_T_old;
    
    %then calc diffusion terms
    d2Tdx2=(Tg(1:end-3,3:end-2) - 2*Tg(2:end-2,3:end-2) + Tg(3:end-1,3:end-2))/dx^2;
    d2Tdy2=(Tg(2:end-2,2:end-3) - 2*Tg(2:end-2,3:end-2) + Tg(2:end-2,4:end-1))/dy^2;
    
    
    %then do ADI to calc T n+1
    
    TRHS=-dt*T_CONV + dt*(D)*(d2Tdx2 + d2Tdy2) + dt*wT;
    TRHS(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;
    
    %dT**
    CTx = (D)*dt/(2*dx^2);
    main_temp = (1+2*CTx)*ones(nx-1,1);
    lower_temp = -CTx*ones(nx-1,1);
    %lower_temp(end)=(-2*Cphix); %this is for neumann condition
    main_temp(end)=(1+CTx); %this is for neumann condition
    upper_temp = -CTx*ones(nx-1,1);
    dT_ss=zeros(nx,ny);
    for j = 1:ny
        dT_ss(2:end,j) = ADRtridiag(main_temp,lower_temp,upper_temp,TRHS(2:end,j));
    end
    
    %dT THIS IS CORRECT SINCE dphi zero above and below
    CTy = (D)*dt/(2*dy^2);
    main_temp = (1+2*CTy)*ones(ny,1);
    lower_temp = -CTy*ones(ny,1);
    upper_temp = -CTy*ones(ny,1);
    dT=zeros(nx,ny);
    for i = 2:nx
        dT(i,:) = ADRtridiag(main_temp,lower_temp,upper_temp,dT_ss(i,:));
    end
    
    T = T + dT; 
    T=max(T,0); %zero out negative temperatures for stability
    T(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0; %unnecessary?
    Tg(2:end-2,3:end-2)=T;
    
    %impose phi at square = 0
    
    %WHERE WOULD I IMPLEMENT BC?

    H_T_old=H_T_new;
    
%%%%%%%%%%%%%%%% T CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% xphi CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %implement neumann boundary
    xphig(end,:)=xphig(end-2,:); %neumann on right
    xphig(end-1,:)=xphig(end-2,:); %neumann on right
    
    xphig(:,end)=xphig(:,end-2); %neumann on top
    xphig(:,end-1)=xphig(:,end-2); %neumann on top
    
    xphig(:,2)=xphig(:,3); %neumann on bottom
    xphig(:,1)=xphig(:,3); %neumann on bottom
    
    %first work out xphi bar for x at each cell
    for i = 2:nx+1 %note this is indexing x 1 in from left since we know it at i = 1 = inlet
        for j = 1:ny
            if u(i,j)>=0
                xphi_x(i,j)=(-1/6)*xphig(i-1,j+2) + (5/6)*xphig(i,j+2) + (2/6)*xphig(i+1,j+2);
            else
                xphi_x(i,j)=(2/6)*xphig(i,j+2) + (5/6)*xphig(i+1,j+2) + (-1/6)*xphig(i+2,j+2);
            end
        end
    end
    
    %then work out T bar for y at each cell
    for i = 1:nx %note this is indexing x 1 in from left
        for j = 1:ny+1
            if v(i,j)>=0
                xphi_y(i,j)=(-1/6)*xphig(i+1,j) + (5/6)*xphig(i+1,j+1) + (2/6)*xphig(i+1,j+2);
            else
                xphi_y(i,j)=(2/6)*xphig(i+1,j+1) + (5/6)*xphig(i+1,j+2) + (-1/6)*xphig(i+1,j+3);
            end
        end
    end
    %then put together with u to get convective term
    
    %set phi to 0 in square
    
    xphi_fu=(u.*xphi_x);
    xphi_fv=(v.*xphi_y);
    H_xphi_new= (xphi_fu(2:end,:)-xphi_fu(1:end-1,:))/dx + (xphi_fv(:,2:end)-xphi_fv(:,1:end-1))/dy;
    
    xphi_CONV = 1.5*H_xphi_new -0.5*H_xphi_old;
    
    %then calc diffusion terms
    d2xphidx2=(xphig(1:end-3,3:end-2) - 2*xphig(2:end-2,3:end-2) + xphig(3:end-1,3:end-2))/dx^2;
    d2xphidy2=(xphig(2:end-2,2:end-3) - 2*xphig(2:end-2,3:end-2) + xphig(2:end-2,4:end-1))/dy^2;
    
    
    %then do ADI to calc xphi n+1
    
    xphiRHS=-dt*xphi_CONV + dt*D*(d2xphidx2 + d2xphidy2) + dt*wxphi;
    xphiRHS(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;
    
    %dxphi**
    Cxphix = D*dt/(2*dx^2);
    main_temp = (1+2*Cxphix)*ones(nx-1,1);
    lower_temp = -Cxphix*ones(nx-1,1);
    %lower_temp(end)=(-2*Cphix); %this is for neumann condition
    main_temp(end)=(1+Cxphix); %this is for neumann condition
    upper_temp = -Cxphix*ones(nx-1,1);
    dxphi_ss=zeros(nx,ny);
    for j = 1:ny
        dxphi_ss(2:end,j) = ADRtridiag(main_temp,lower_temp,upper_temp,xphiRHS(2:end,j));
    end
    
    %dxphi THIS IS CORRECT SINCE dphi zero above and below
    Cxphiy = D*dt/(2*dy^2);
    main_temp = (1+2*Cxphiy)*ones(ny,1);
    lower_temp = -Cxphiy*ones(ny,1);
    upper_temp = -Cxphiy*ones(ny,1);
    dxphi=zeros(nx,ny);
    for i = 2:nx
        dxphi(i,:) = ADRtridiag(main_temp,lower_temp,upper_temp,dxphi_ss(i,:));
    end
    
    xphi = xphi + dxphi;
    xphi(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0; %unnecessary?
    xphig(2:end-2,3:end-2)=xphi;
    
    %impose phi at square = 0
    
    %WHERE WOULD I IMPLEMENT BC?

    H_xphi_old=H_xphi_new;
    
%%%%%%%%%%%%%%%% xphi CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% %Step 1: compute viscous, convective terms in u and v, then compute ustar    
    %assumes u is nx+1 by ny, v is nx by ny+1, no ghost cells anywhere, 
    %works out u momentum terms for nx, ny to add to r right inner region of u domain
    %works out v momentum terms for nx, ny-1 to add to vertically inner region of v
    %domain
    
    %1DONE here need a weird condition, or not with ghost cell
    duudx=((ug(2:end-1,2:end-1) + ug(3:end,2:end-1)).^2 - (ug(1:end-2,2:end-1) + ug(2:end-1,2:end-1)).^2) /(4*dx);
    
    %1DONE here need ug and vg
    duvdy=((ug(2:end-1,3:end) + ug(2:end-1,2:end-1)).*(vg(3:end,2:end) + vg(2:end-1,2:end)) - (ug(2:end-1,2:end-1) + ug(2:end-1,1:end-2)).*(vg(3:end,1:end-1) + vg(2:end-1,1:end-1)))/(4*dy);
    
    %1DONE here need a weird condition, or not with ghost cell
    d2udx2=(ug(3:end,2:end-1) - 2*ug(2:end-1,2:end-1) + ug(1:end-2,2:end-1))/dx^2;
    
    %1DONE need ug
    d2udy2=(ug(2:end-1,3:end) - 2*ug(2:end-1,2:end-1) + ug(2:end-1,1:end-2))/dy^2;
    
    
    %1DONE only need plain v, no changes
    dvvdy=((v(:,2:end-1) + v(:,3:end)).^2 - (v(:,1:end-2) + v(:,2:end-1)).^2)/(4*dy);
    
    %1DONE need vg, u
    duvdx=((u(2:end,2:end) + u(2:end,1:end-1)).*(vg(3:end,2:end-1) + vg(2:end-1,2:end-1)) - (u(1:end-1,2:end) + u(1:end-1,1:end-1)).*(vg(2:end-1,2:end-1) + vg(1:end-2,2:end-1)))/(4*dx);
    
    %1DONE need vg, no changes
    d2vdx2=(vg(3:end,2:end-1) - 2*vg(2:end-1,2:end-1) + vg(1:end-2,2:end-1))/dx^2;
    
    %1DONE need v, no changes
    d2vdy2=(v(:,3:end) - 2*v(:,2:end-1) + v(:,1:end-2))/dy^2;
    
    H_u_old=H_u_new;
    H_u_new=duudx + duvdy;
    
    H_v_old=H_v_new;
    H_v_new=duvdx + dvvdy;
    
    RHS1u(2:end,:) = 0.5*dt*(-3*H_u_new + H_u_old) + (dt/Re)*(d2udx2 + d2udy2);
    RHS1v(:,2:end-1) = 0.5*dt*(-3*H_v_new + H_v_old) + (dt/Re)*(d2vdx2 + d2vdy2);

  
%%%IMPOSE 0 SQUARE VELOCITY    
    RHS1u(l_edge_s+1:r_edge_s,b_edge_s+1:t_edge_s-1)=0;
    RHS1v(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s)=0;
    
    %du**
    Cux = dt/(2*Re*dx^2);
    main_temp = (1+2*Cux)*ones(nx,1);
    %main_temp(end)=(1+Cux); %this is for neumann condition
    lower_temp = -Cux*ones(nx,1);
    lower_temp(end)=(-2*Cux); %this is for neumann condition
    upper_temp = -Cux*ones(nx,1);
    du_ss=zeros(size(RHS1u(2:end,:)));
    for j = 1:ny
        du_ss(:,j) = ADRtridiag(main_temp,lower_temp,upper_temp,RHS1u(2:end,j));
    end
    
    %du* %imposing 0 for neumann conditions since velocity doesn't change
    %outside
    Cuy = dt/(2*Re*dy^2);
    main_temp = (1+2*Cuy)*ones(ny,1);
    lower_temp = -Cuy*ones(ny,1);
    upper_temp = -Cuy*ones(ny,1);
    du_s=zeros(size(RHS1u(2:end,:)));
    for i = 1:nx
        du_s(i,:) = ADRtridiag(main_temp,lower_temp,upper_temp,du_ss(i,:));
    end
    
    %dv**
    Cvx = dt/(2*Re*dx^2);
    main_temp = (1+2*Cvx)*ones(nx,1);
    main_temp(end) = (1+Cvx);
    lower_temp = -Cvx*ones(nx,1);
    upper_temp = -Cvx*ones(nx,1);
    dv_ss=zeros(size(RHS1v(:,2:end-1)));
    for j = 1:ny-1
        dv_ss(:,j) = ADRtridiag(main_temp,lower_temp,upper_temp,RHS1v(:,j+1)); %make nx, ny-1
    end
    
    %dv*
    Cvy = dt/(2*Re*dy^2);
    main_temp = (1+2*Cvy)*ones(ny-1,1);
    lower_temp = -Cvy*ones(ny-1,1);
    upper_temp = -Cvy*ones(ny-1,1);
    dv_s=zeros(size(RHS1v(:,2:end-1)));
    for i = 1:nx
        dv_s(i,:) = ADRtridiag(main_temp,lower_temp,upper_temp,dv_ss(i,:));
    end
    
    %CHANGE THIS FOR STEP 1
    %ustar(2:end,:) = RHS1u(2:end,:);
    ustar(2:end,:) = u(2:end,:) + du_s;
    %vstar(:,2:end-1) = RHS1v(:,2:end-1);
    vstar(:,2:end-1) = v(:,2:end-1) + dv_s;
    
%     ustar(end,:) = max(0, ustar(end,:));

    %Step 3.5: do mass conservation on right
    curr=sum(ustar(end,:)); %current velocity out
    diff=sum(ustar(1,:))-curr; %difference between velocity in and out
    add=diff/ny; %how much to add per cell
    ustar(end,:)=ustar(end,:)+add; %adding missing velocity
    
    %impose nothing in/around square
    ustar(l_edge_s+1:r_edge_s,b_edge_s+1:t_edge_s-1)=0;
    vstar(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s)=0;
    p(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;
    pg_GS(2:end-1,2:end-1)=p;
    
%Step 2: use ustar to compute p_n+1 field
    tic
    f = ((ustar(2:end,:) - ustar(1:end-1,:))/dx + (vstar(:,2:end) - vstar(:,1:end-1))/dy)/dt;
    %~~~~DO I UPDATE F HERE TO IMPOSE CONDITION AGAIN? I DONT THINK SO, F
    %SHOULD BE ZERO FROM EQ 4 IF IMPLEMENTED CORRECTLY
    %sum(sum(f)) %check this is zero each iteration

    fnorm=norm(f,inf);
    
    dg_k=zeros(nx+2,ny+2);
    e_k=zeros(nx,ny);
    p_GS=zeros(nx,ny);

    k=0;
    r_0=f;
    rho_k=norm(r_0,'fro')^2;
    criterion=norm(r_0,'fro')*eps;
    while sqrt(rho_k) > criterion
        k=k+1;
        if k == 1
            d_k=r_0;
            p_km1=p_GS;
            r_km1=r_0;
            rho_km1=rho_k;
        else
            beta_km1=rho_km1/rho_km2;
            d_k=r_km1 + beta_km1*d_km1;
        end
        dg_k(2:end-1,2:end-1)=d_k;
        e_k=zeros(nx,ny);
        for i = 1:nx
            for j = 1:ny
                e_k(i,j)=pmultiplier(i,j)*(lambda*(dg_k(i+2,j+1)+dg_k(i,j+1)+dg_k(i+1,j+2)+dg_k(i+1,j)-lapcoeff(i,j)*d_k(i,j)));
            end
        end
        dtemp=reshape(d_k,1,[]);
        etemp=reshape(e_k,1,[]);
        a_k=rho_km1/(dtemp*etemp'); %step length in new direction
        p_k=p_km1 + a_k*d_k; %update solution
        r_k=r_km1 - a_k*e_k; %update residual
        rho_k=norm(r_k,'fro')^2; %2 norm of residual


        %update stuff for next loop
        r_km1=r_k;
        rho_km2=rho_km1;
        rho_km1=rho_k;
        p_km1=p_k;
        d_km1=d_k;

        %err_CG(k)=sqrt(rho_k); for debugging
    end
    %err_CG=err_CG/norm(r_0,2); for debugging
    p=p_k;
    pg_GS(2:end-1,2:end-1)=p;
    toc
    
%Step 3: use p_n+1 field and ustar to compute u_n+1
    pg_GS(end,:)=pg_GS(end-1,:);
    u(2:end,:) = ustar(2:end,:) - dt*(pg_GS(3:end,2:end-1) - pg_GS(2:end-1,2:end-1))/dx;
    v(:,2:end-1) = vstar(:,2:end-1) - dt*(p(:,2:end) - p(:,1:end-1))/dy;
    pg_GS(end,:)=0;
    
    %ZERO OUT NEGATIVE VELOCITIES
    u(end,:)=max(u(end,:),0);
    
    %impose square BC
    u(l_edge_s+1:r_edge_s,b_edge_s+1:t_edge_s-1)=0;
    v(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s)=0;
    p(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;
    
%Step 4: update augmented domain 
    ug(1:end-1,2:end-1)=u;
    vg(2:end-1,:)=v;
    vg(end,:)=vg(end-1,:); %update vg with inner v to impose Neumann on RHS
    ug(end,:)=ug(end-1,:); %to impose neumann on RHS
    
    C_l(iter)=sum(p(l_edge_s+1:r_edge_s-1,t_edge_s)-p(l_edge_s+1:r_edge_s-1,b_edge_s))/nx_s;
    C_d(iter)=sum(p(l_edge_s,b_edge_s+1:t_edge_s-1)-p(r_edge_s,b_edge_s+1:t_edge_s-1))/ny_s;
    
    %THEN UPDATE VELOCITY ON RHS TO ENSURE MASS CONSERVATION ?
    disp(iter)
    %pause
    %toc
    if mod(iter,100)==0
        filename = ['T_' 'Re_' num2str(Re) 'Iter_' num2str(iter,'%05d') '.mat'];
        save(filename , 'T' , 't')
    end	

    if mod(iter,500)==0
        filename = ['backupoutputs_' 'Re_' num2str(Re) 'Iter_' num2str(iter,'%05d') '.mat'];
        save(filename)
    end	
    if t(iter+1) > tf
        break
    end
end
c2=clock;

%% plot stuff
%close all
Vx=linspace(0,L,nx);
Vy=linspace(-0.5*W,0.5*W,ny);
[X,Y] = meshgrid(Vx,Vy);
%% Contour Dye Plot
figure()
contourf(X,Y,phi',40,'LineStyle','none')
axis equal
xlabel('$x/H$','Interpreter','Latex','FontSize',18)
ylabel('$y/H$','Interpreter','Latex','FontSize',18)
saveas(gcf,['phi_' 'Re_' num2str(Re) 'grid_' num2str(nx) 'Iter_' num2str(iter) 'H_' num2str(H) 'U_' num2str(U) 'dt_' num2str(dt) 'eps_' num2str(eps) '.png'])



%% Lift and Drag Coefficients

figure()
plot(t(1:iter+1),C_l(1:iter+1));
xlabel('$tH/U$','Interpreter','Latex','FontSize',18)
ylabel('$C_L$','Interpreter','Latex','FontSize',18)
xlim([0 200])
saveas(gcf,['CL_' 'Re_' num2str(Re) 'grid_' num2str(nx) 'Iter_' num2str(iter) 'H_' num2str(H) 'U_' num2str(U) 'dt_' num2str(dt) 'eps_' num2str(eps) '.png'])


figure()
plot(t(1:iter-2),C_d(2:iter-1));
xlabel('$tH/U$','Interpreter','Latex','FontSize',18)
ylabel('$C_D$','Interpreter','Latex','FontSize',18)
xlim([0 200])
saveas(gcf,['CD_' 'Re_' num2str(Re) 'grid_' num2str(nx) 'Iter_' num2str(iter) 'H_' num2str(H) 'U_' num2str(U) 'dt_' num2str(dt) 'eps_' num2str(eps) '.png'])


%% Velocity Field (not saved)
figure()
fig=quiver(u(1:end-1,:)',v(:,1:end-1)');
xlim([1 nx])
ylim([1 ny])
set(fig,'AutoScale','on','AutoScaleFactor',10);

%% Streamlines (not saved)
figure()

[StartX,StartY] = meshgrid(linspace(0,L,15),linspace(0,W,15));
streamline(X,Y,u(1:end-1,:)',v(:,1:end-1)',StartX,StartY)
xlabel('$x/H$ (Normalized x-position)','Interpreter','Latex','FontSize',18)
ylabel('$y/H$ (Normalized y-position)','Interpreter','Latex','FontSize',18)
%saveas(gcf,['Streamlines_' 'Re_' num2str(Re) 'grid_' num2str(nx) 'Iter_' num2str(iter) 'H_' num2str(H) 'U_' num2str(U) 'dt_' num2str(dt) 'eps_' num2str(eps) '.png'])

%% save outputs
filename=['Outputs_' 'Re_' num2str(Re) 'grid_' num2str(nx) 'Iter_' num2str(iter) 'H_' num2str(H) 'U_' num2str(U) 'dt_' num2str(dt) 'eps_' num2str(eps) '.mat'];
save(filename)