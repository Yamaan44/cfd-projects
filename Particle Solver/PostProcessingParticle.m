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

%% Visualizer

sSuffix='.mat';
iFileStart=5;
iFileEnd=1140;
sPath='/Users/YA1/Documents/MATLAB/MAE 6230/Particle-Laden Jet Solver/Re_1000_St_100/';
sPrefix='phiParticles_Re_1000_St_100Iter_';

h=figure()
Vx=linspace(0,2*H,nx);
Vy=linspace(-0.5*H,0.5*H,ny);
[X,Y] = meshgrid(Vx,Vy);

i=iFileStart;
sName=num2str(i,'%05d');
fname=[sPath sPrefix sName sSuffix]
load(fname);
time=t(find(t,1,'last'));


contourf(X,Y,phi',40,'LineStyle','none')
set(gca,'fontsize',18)
axis image
axis equal
axis xy
hold on 
plot(part_data(:,1),part_data(:,2),'ro')
hold off

xlabel('$x/H$','Interpreter','Latex','FontSize',18)
ylabel('$y/H$','Interpreter','Latex','FontSize',18)
title(['tU/H = ' num2str(time,'%2d')],'Interpreter','Latex','FontSize',18)
drawnow
gif('ParticleSimulationRe1000St100.gif','DelayTime',0.1,'LoopCount',5,'frame',gcf)

for i=iFileStart:5:iFileEnd
   sName=num2str(i,'%05d');
   fname=[sPath sPrefix sName sSuffix]
   load(fname);
   time=t(find(t,1,'last'));
   
   
   contourf(X,Y,phi',40,'LineStyle','none')
   set(gca,'fontsize',18)
   %colorbar
   axis image
   axis equal
   axis xy
   hold on
   plot(part_data(:,1),part_data(:,2),'ro')
   hold off
   xlabel('$x/H$','Interpreter','Latex','FontSize',18)
   ylabel('$y/H$','Interpreter','Latex','FontSize',18)
   title(['$tU/H = $' num2str(time,'%3d')],'Interpreter','Latex','FontSize',18)
   drawnow
   gif
   %pause(0.002)
   
end
