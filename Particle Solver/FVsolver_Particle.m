clear %%THIS USES CG
%% Constants CHECKED
Re=1000; %20 and 400 
H=1; 
U=1; %inputs
eps=10e-5; %error for CG solver
dt_guess=0.01;
tf=200;

R=0.05;
St=100;

n_max=100;
max_inject=25;
n_curr=0;
part_data=zeros(n_max,4); %x pos, y pos, u, v
zeroarray=1:n_max;

D=10^-4;
L = 2*H;
W = H;

dx=H/64; dy=H/64; lambda=1/dx^2; 
nx=L/dx; 
ny=W/dy; %number of cells in x,y

Niter=ceil(200/dt_guess); %number of iterations

CFL=0.8;
%% Initialize Arrays
y=linspace(-0.5*W,0.5*W,ny);

p=zeros(nx,ny);

u=zeros(nx+1,ny); %u array

u_inlet=exp(-((y)/R).^2);

u(1,:)=u_inlet;

ustar=u;
ug=zeros(nx+2,ny+2); %ghost cell array %adding 1 col to right, 1 above and below
ug(1,2:end-1)=U; %set left boundary to U=1 
ug(1:end-1,2:end-1) = u;%update ug with inner u
ug(end,:) = ug(end-2,:); %to impose neumann on RHS

% %Step 3.5: do mass conservation on right to ensure f in first timestep is 0
% curr=sum(u(end,:)); %current velocity out
% diff=sum(u(1,:))-curr; %difference between velocity in and out
% add=diff/ny; %how much to add per cell
% u(end,:)=u(end,:)+add; %adding missing velocity


v=zeros(nx,ny+1); %v array
vstar=v;
vg=zeros(nx+2,ny+1); %ghost cell array %adding 1 col to left and right, nothing above/below
vg(2:end-1,:)=v;
vg(end,:)=vg(end-1,:); %update vg with inner v

H_u_new=zeros(nx,ny);
H_v_new=zeros(nx,ny-1);
RHS1u=zeros(nx+1,ny);
RHS1v=zeros(nx,ny+1);

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

%% Square geometry CHECKED
% 
% L_s=H;
% W_s=H;
% nx_s=L_s/dx;
% ny_s=W_s/dy;
% 
% l_edge_s=5*H/dx;
% r_edge_s=l_edge_s+nx_s+1;
% b_edge_s=2*H/dy;
% t_edge_s=b_edge_s+ny_s+1;
% 
% %impose velocity = 0 inside and on boundary of square for no slip and 0
% %pressure inside
% u(l_edge_s+1:r_edge_s,b_edge_s+1:t_edge_s-1)=0;
% v(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s)=0;
% p(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;

%laplace multiplier to make sure pressure in the square evaluates to 0
pmultiplier=ones(nx,ny);
% pmultiplier(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;


%%  Laplacian coefficient array CHECKED
lapcoeff=4*ones(nx,ny);
lapcoeff(:,1)=3; lapcoeff(:,end)=3; lapcoeff(1,:)=3; lapcoeff(end,:)=3;
lapcoeff(1,1)=2; lapcoeff(1,end)=2; lapcoeff(end,1)=2; lapcoeff(end,end)=2;

% lapcoeff(l_edge_s,b_edge_s+1:t_edge_s-1)=3; %left side
% lapcoeff(r_edge_s,b_edge_s+1:t_edge_s-1)=3; %right side
% lapcoeff(l_edge_s+1:r_edge_s-1,t_edge_s)=3; %top side
% lapcoeff(l_edge_s+1:r_edge_s-1,b_edge_s)=3; %bottom side

%% Grid definition for scalars
Vx=linspace(0,2*H,nx);
Vy=linspace(-0.5*H,0.5*H,ny);
[X,Y] = meshgrid(Vx',Vy');


%% Grid definition for u velocity
Vx_u=linspace(0,2*H,nx+1);
Vy_u=linspace(-0.5*H + dy/2,0.5*H - dy/2,ny);
[X_u,Y_u] = meshgrid(Vx_u',Vy_u');

%% Grid definition for v velocity
Vx_v=linspace(0+dx/2,2*H-dx/2,nx);
Vy_v=linspace(-0.5*H,0.5*H,ny+1);
[X_v,Y_v] = meshgrid(Vx_v,Vy_v);

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
    
% %%%%%%%%%%%%%%%% PHI CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%     phiRHS(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;
    
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
%     phi(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0; %unnecessary?
    phig(2:end-2,3:end-2)=phi;
    
    %impose phi at square = 0
    
    %WHERE WOULD I IMPLEMENT BC?

    H_phi_old=H_phi_new;
    
%%%%%%%%%%%%%%%% PHI CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    RHS1v(1,:)=0; %IMPOSING dirichlet V
  
%%%IMPOSE 0 SQUARE VELOCITY    
%     RHS1u(l_edge_s+1:r_edge_s,b_edge_s+1:t_edge_s-1)=0;
%     RHS1v(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s)=0;
    
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
    
    dv_s(1,:)=0; %IMPOSING DIRICHLET V
    
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
    
%     %impose nothing in/around square
%     ustar(l_edge_s+1:r_edge_s,b_edge_s+1:t_edge_s-1)=0;
%     vstar(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s)=0;
%     p(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;
    pg_GS(2:end-1,2:end-1)=p;
    
%Step 2: use ustar to compute p_n+1 field
    tic
    f = ((ustar(2:end,:) - ustar(1:end-1,:))/dx + (vstar(:,2:end) - vstar(:,1:end-1))/dy)/dt;
    %~~~~DO I UPDATE F HERE TO IMPOSE CONDITION AGAIN? I DONT THINK SO, F
    %SHOULD BE ZERO FROM EQ 4 IF IMPLEMENTED CORRECTLY
    %sum(sum(f)) %check this is zero each iteration
    
    
    
%     for i = 1:nx
%         for j = 1:ny
%             pg_GS(i+1,j+1)=pmultiplier(i,j)*((pg_GS(i+2,j+1)+pg_GS(i,j+1)+pg_GS(i+1,j+2)+pg_GS(i+1,j))-f(i,j)/lambda)/lapcoeff(i,j);
%         end
%     end
%     p_GS=pg_GS(2:end-1,2:end-1);
%     p_GS(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;
%     pg_GS(2:end-1,2:end-1)=p_GS;
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
%     while check ~=1 %&& k < 20000
%         k=k+1;
%         for i = 1:nx
%             for j = 1:ny
%                 pg_GS(i+1,j+1)=pmultiplier(i,j)*((1-w)*pg_GS(i+1,j+1) + w*((pg_GS(i+2,j+1)+pg_GS(i,j+1)+pg_GS(i+1,j+2)+pg_GS(i+1,j))-f(i,j)/lambda)/lapcoeff(i,j));
%             end
%         end
% 
%         for i = 1:nx
%             for j = 1:ny
%                 lap_GS(i,j)=pmultiplier(i,j)*(lambda*(pg_GS(i+2,j+1)+pg_GS(i,j+1)+pg_GS(i+1,j+2)+pg_GS(i+1,j)-lapcoeff(i,j)*pg_GS(i+1,j+1)));
%             end
%         end
%         
%         %impose no pressure in square
%         lap_GS(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;
% 
%         err_GS(k)=norm(lap_GS-f,inf)/fnorm;
% 
%         if err_GS(k) < eps
%             check = 1;
%         end 
%         
%     end
%    p=pg_GS(2:end-1,2:end-1);
    
%Step 3: use p_n+1 field and ustar to compute u_n+1
    pg_GS(end,:)=pg_GS(end-1,:);
    u(2:end,:) = ustar(2:end,:) - dt*(pg_GS(3:end,2:end-1) - pg_GS(2:end-1,2:end-1))/dx;
    v(:,2:end-1) = vstar(:,2:end-1) - dt*(p(:,2:end) - p(:,1:end-1))/dy;
    v(1,:)=0; %IMPOSING DIRICHLET V
    pg_GS(end,:)=0;
    
    %impose square BC
%     u(l_edge_s+1:r_edge_s,b_edge_s+1:t_edge_s-1)=0;
%     v(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s)=0;
%     p(l_edge_s+1:r_edge_s-1,b_edge_s+1:t_edge_s-1)=0;
    
%Step 4: update augmented domain 
    ug(1:end-1,2:end-1)=u;
    vg(2:end-1,:)=v;
    vg(end,:)=vg(end-1,:); %update vg with inner v to impose Neumann on RHS
    ug(end,:)=ug(end-1,:); %to impose neumann on RHS
    
%     C_l(iter)=sum(p(l_edge_s+1:r_edge_s-1,t_edge_s)-p(l_edge_s+1:r_edge_s-1,b_edge_s))/nx_s;
%     C_d(iter)=sum(p(l_edge_s,b_edge_s+1:t_edge_s-1)-p(r_edge_s,b_edge_s+1:t_edge_s-1))/ny_s;
%     
    %THEN UPDATE VELOCITY ON RHS TO ENSURE MASS CONSERVATION ?
    
    
    
    %%%%%%%% PARTICLE CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %determine number of particles to be injected into domain
    n_inject=max(n_max-n_curr,0);
    
    %finds index with first zero value
%     b=part_data(:,3)==0; %looking at u velocity because this shouldnt ever be 0?
%     [whatever, idx]=max(b'./zeroarray);
    
    %create new particle array based on current particles (note we're just
    %adding to zero entries of current array
    y_part=0.1*(rand(n_inject,1)-0.5);
    u_part=exp(-(y_part/R).^2);
    part_data(n_curr+1:end,2)=y_part;
    part_data(n_curr+1:end,3)=u_part;
%     part_data(idx:end,2)=y_part;
%     part_data(idx:end,3)=u_part;
%     for i = 1:n_inject
%         y_in=y_part(i);
%         v_in=0;
%         u_in=u_part(i);
%         x_in=0;
%         
%         part(i,1)=x_in;
%         part(i,2)=y_in;
%         part(i,3)=u_in;
%         part(i,4)=v_in;
%         
%         
%         
%         
%         
%     end
    
    for i = 1:n_max
    %interpolate fluid velocity to location of each particle
    
%         Xp=1+floor(part_data(i,1)/dx); %maybe double check this
%         Yp=ny/2 + ceil(part_data(i,2)/dy);
        U_f = interp2(X_u,Y_u,u',part_data(i,1),part_data(i,2),'linear'); %have i got this right?
        V_f = interp2(X_v,Y_v,v',part_data(i,1),part_data(i,2),'linear',0); %have i got this right?
        
        U_p=part_data(i,3);
        V_p=part_data(i,4);
    %advance each particle in time using RK O(3)
        k1 = (1/St)*(U_p-U_f);
        k2 = (1/St)*((U_p + 0.5*dt*k1)-U_f);
        k3 = (1/St)*((U_p + 0.5*dt*k2)-U_f);
        k4 = (1/St)*((U_p + dt*k3)-U_f);
        U_p = U_p + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        
        k1 = (1/St)*(V_p-V_f);
        k2 = (1/St)*((V_p + 0.5*dt*k1)-V_f);
        k3 = (1/St)*((V_p + 0.5*dt*k2)-V_f);
        k4 = (1/St)*((V_p + dt*k3)-V_f);
        V_p = V_p + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        
        part_data(i,1) = part_data(i,1) + dt*U_p;
        part_data(i,2) = part_data(i,2) + dt*V_p;
        part_data(i,3) = U_p;
        part_data(i,4) = V_p;
        
        
        
    end
    
    
    %determine if any particles have left the domain
    
    remove = (part_data(:,1) > L) + (part_data(:,1) < 0) + (part_data(:,2) > 0.5*W) + (part_data(:,2) < -0.5*W); %remove will give 1s for particles outside of domain
    
    %remove their data and collapse to leave zeros on end, basically sort
    %the array
    for i = n_max:-1:1
        if remove(i) > 0
            part_data(i,:)=[];
        end
    end
    
    %update current number of particles and n_curr
    n_curr=sum(~remove);
    part_data(end+1:n_max,:)=zeros(n_max-n_curr,4);
    
    
    if mod(iter,5)==0
        
        
        contourf(X,Y,phi',40,'LineStyle','none');
        hold on
        xlim([0 L])
        ylim([-0.5*W 0.5*W])
        plot(part_data(:,1),part_data(:,2),'ro');
        drawnow
        hold off
    end	 
    
    %%%%%%%% PARTICLE CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(iter)
    %pause
    %toc
    if mod(iter,5)==0
        filename = ['phiParticles_' 'Re_' num2str(Re) '_St_' num2str(St) 'Iter_' num2str(iter,'%05d') '.mat'];
        save(filename , 'phi' , 't','part_data')
    end	

    if mod(iter,2500)==0
        filename = ['backupoutputs_' 'Re_' num2str(Re) 'Iter_' num2str(iter,'%05d') '.mat'];
        save(filename)
    end	
%     if fuck == 1
%         break
%     end
    if t(iter+1) > 200
        break
    end
end
c2=clock;

%% plot stuff
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

[StartX,StartY] = meshgrid(linspace(0,25,15),linspace(0,5,15));
streamline(X,Y,u(1:end-1,:)',v(:,1:end-1)',StartX,StartY)
xlabel('$x/H$ (Normalized x-position)','Interpreter','Latex','FontSize',18)
ylabel('$y/H$ (Normalized y-position)','Interpreter','Latex','FontSize',18)
%saveas(gcf,['Streamlines_' 'Re_' num2str(Re) 'grid_' num2str(nx) 'Iter_' num2str(iter) 'H_' num2str(H) 'U_' num2str(U) 'dt_' num2str(dt) 'eps_' num2str(eps) '.png'])

%% save outputs
filename=['Outputs_' 'Re_' num2str(Re) 'grid_' num2str(nx) 'Iter_' num2str(iter) 'H_' num2str(H) 'U_' num2str(U) 'dt_' num2str(dt) 'eps_' num2str(eps) '.mat'];
save(filename)
