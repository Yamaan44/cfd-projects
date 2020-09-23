%% only need to run this once, because i fucked up with naming files. use this snippet to remove things on the end of file names
% files = dir;
% directoryNames = {files.name};
% directoryNamesStr = directoryNames(3:end)';
% %directoryNamesArray = zeros(size(directoryNamesStr));
% for i = 1:length(directoryNamesStr)
%     name=directoryNamesStr{i,1}
%     s=[name]
%     ix=strfind(s,'_')  % get the underscore locations
%     t=s(1:ix(4)-3)
%     newname=[t '.mat']
%     movefile(s,newname)
% end

%% only need to run this once, use this snippet to make all doubles 5 long
% files = dir;
% directoryNames = {files.name};
% directoryNamesStr = directoryNames(3:end)';
% %directoryNamesArray = zeros(size(directoryNamesStr));
% for i = 1:length(directoryNamesStr)
%     name=directoryNamesStr{i,1}
%     s=[name(1:end-4)]
%     ix=strfind(s,'_')  % get the underscore locations
%     t=s(ix(3)+1:end)
%     t=str2num(t)
%     t=num2str(t,'%05d')
%     newname=[s(1:ix(3)) t '.mat']
%     if movefile(name,newname)
%         1
%     end
% end

%% Colour definitions
cobalt=[19/255,56/255,189/255]; %2500
arctic=[130/255,237/255,253/255]; %2200
butter=[254/255,226/255,39/255]; %2000
butterscotch=[252/255,188/255,2/255]; %1800
tiger=[252/255,107/255,2/255]; %1500
carrot=[237/255,113/255,23/255]; %1200
rosered=[226/255,37/255,43/255]; %1000
candy=[210/255,21/255,2/255]; %800
blood=[113/255,12/255,4/255]; %500
grey=[108/255,98/255,109/255]; %300

lo=[0; 0; 0]';
hi=[255; 255; 255]';

map = [cobalt; arctic; butter; butterscotch; tiger; carrot; grey ;rosered; candy; blood];
map = flipud(map);

%clev=[150; 400; 500; 800; 1000; 1200; 1500; 1800; 2000; 2200];
%clev=[150; 200; 250; 300; 350; 400; 500; 800; 1200; 2200];
clev=[150; 200; 250; 297; 303; 350; 400; 500; 800; 1200; 2200];


cmap=map;
%% Visualizer

sSuffix='00.mat';
iFileStart=1;
iFileEnd=250;
sPath='/Users/YA1/Documents/MATLAB/MAE 6230/LES Turbulence Solver/Outputs_Re_10000/';
sPrefix='phi_Re_10000Iter_';

h=figure()
Vx=linspace(0,25*H,nx);
Vy=linspace(-2.5*H,2.5*H,ny);
[X,Y] = meshgrid(Vx,Vy);

i=iFileStart;
sName=num2str(i,'%03d');
fname=[sPath sPrefix sName sSuffix]
load(fname);
time=t(find(t,1,'last'));


contourf(X,Y,phi',40,'LineStyle','none')
%colormap(map)
%caxis([0 2500])
%h=contourfcmap(X,Y,T',clev,cmap,lo,hi);

set(gca,'fontsize',18)
axis image
axis equal
axis xy
xlabel('$x/H$','Interpreter','Latex','FontSize',18)
ylabel('$y/H$','Interpreter','Latex','FontSize',18)
title(['tU/H = ' num2str(time,'%2d')],'Interpreter','Latex','FontSize',18)
drawnow
gif('LESSimulationRe10000.gif','DelayTime',0.1,'LoopCount',5,'frame',gcf)

for i=iFileStart+1:iFileEnd
   sName=num2str(i,'%03d');
   fname=[sPath sPrefix sName sSuffix]
   load(fname);
   time=t(find(t,1,'last'));
   
   
   contourf(X,Y,phi',40,'LineStyle','none')
   %colormap(map)
   %caxis([0 2500])
   %h=contourfcmap(X,Y,T',clev,cmap,lo,hi); %use this for better color but
   %with lines

   set(gca,'fontsize',18)
   %colorbar
   axis image
   axis equal
   axis xy
   xlabel('$x/H$','Interpreter','Latex','FontSize',18)
   ylabel('$y/H$','Interpreter','Latex','FontSize',18)
   title(['$tU/H = $' num2str(time,'%3d')],'Interpreter','Latex','FontSize',18)
   drawnow
   gif
   %pause(0.002)
   
end