%%% Melt rate figures for paper

Nt = 60;
t1min = 60;
t1hour = 60*t1min;
t1day = 24*t1hour;
m1km = 1000;
rhoShelfIce = 920; 

scrsz = get(0,'ScreenSize');
fontsize = 26;
framepos = [0 scrsz(4)/2 scrsz(3) scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(8);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');
hold on;





%%% for drag 
cd /Users/Kayliec23/Desktop/MISSI/experiments/drag

load('Cd2.5e-4/fwFlux_Cd2.5e-4.mat'); load('Cd2.5e-4/Diag_Cd2.5e-4.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0 1 0.25])

load('Cd7.5e-4/fwFlux_Cd7.5e-4.mat'); load('Cd7.5e-4/Diag_Cd7.5e-4.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0 1 0.5])


load('Cd7.5e-3/fwFlux_Cd7.5e-3.mat'); load('Cd7.5e-3/Diag_Cd7.5e-3.mat'); 
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0 1 0.75])


load('Cd2.5e-2/fwFlux_Cd2.5e-2.mat'); load('Cd2.5e-2/Diag_Cd2.5e-2.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0 1 1]);





%%% for dS
cd /Users/Kayliec23/Desktop/MISSI/experiments/dS

load('dS0.4/fwFlux_dS0.4.mat'); load('dS0.4/Diag_dS0.4.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0.5 0 0.25]);


load('dS0.6/fwFlux_dS0.6.mat'); load('dS0.6/Diag_dS0.6.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0.5 0 0.5]);


load('dS0.8/fwFlux_dS0.8.mat'); load('dS0.8/Diag_dS0.8.mat'); 
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0.5 0 0.75]);


load('dS1/fwFlux_dS1.mat'); load('dS1/Diag_dS1.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0.5 0 1]);





%%% for Tforcing
cd /Users/Kayliec23/Desktop/MISSI/experiments/Tforcing

%load('/Users/Kayliec23/Desktop/MISSI/experiments/Tforcing/Tneg1/fwFlux_Tneg1.mat')
%m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
%plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0.5 0 0.5 0.25]); 

load('Tneg0.5/fwFlux_Tneg0.5.mat'); load('Tneg0.5/Diag_Tneg05.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0.5 0 0.5 0.25]); 

load('T0/fwFlux_T0.mat'); load('T0/Diag_T0.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0.5 0 0.5 0.5]); 

load('T1/fwFlux_T1.mat'); load('T1/Diag_T1.mat'); 
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0.5 0 0.5 0.75]); 

load('T1.5/fwFlux_T1.5.mat'); load('T1.5/Diag_T15.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0.5 0 0.5 1]);




%%% for slope size
cd /Users/Kayliec23/Desktop/MISSI/experiments/slope

load('Ly10km/fwFlux_Ly10km.mat'); load('Ly10km/Diag_Ly10km.mat'); 
m = squeeze(sum(fwFlux(1,:,1:30),3)*t1day/rhoShelfIce);
plot(yy(1:100)/m1km,-smooth(m(1,1:100)/max(abs(m)),7),'linewidth',2,'color',[1 0.6445 0 0.25]);  

load('Ly126km_dy100/fwFlux_Ly126km_dy100.mat'); load('Ly126km_dy100/Diag_Ly126km_dy100.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[1 0.6445 0 0.5]);

load('Ly200km_dy100/fwFlux_Ly200km_dy100.mat'); load('Ly200km_dy100/Diag_Ly200km_dy100.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[1 0.6445 0 0.75]);

load('Ly400km_highres/fwFlux_Ly400km_highres.mat'); load('Ly400km_highres/Diag_Ly400km_highres.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[1 0.6445 0 1]);




%%% and shape
% cd /Users/Kayliec23/Desktop/MISSI/experiments/slope
% load('linear_ref/fwFlux_linear_ref.mat'); load('linear_ref/Diag_linear_ref.mat');
% m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
% plot(yy/m1km,-smooth(m/max(abs(m))),7),'linewidth',2,'color',[0.8203    0.7031    0.5469]);  

cd /Users/Kayliec23/Desktop/MISSI/experiments/resolution
load('dz5/fwFlux_dz5.mat'); load('dz5/Diag_dz5.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'--k','linewidth',2); 


cd /Users/Kayliec23/Desktop/MISSI/experiments/
load('ref2D/fwFlux_ref2D.mat'); load('ref2D/Diag_ref2D.mat');
m = squeeze(sum(fwFlux(1,:,1:Nt),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'k','linewidth',3.5); 

xlim([0 40])
set(gca,'fontsize',42);
xlabel('y (km)');
ylabel('normalized total melt');

legend('C_d = 2.5e-4','C_d = 7.5e-4','C_d = 7.5e-3','C_d = 2.5e-2',...
'\DeltaS = 0.4','\DeltaS = 0.6','\DeltaS = 0.8','\DeltaS = 1',...
'T_a = -1.5','T_a = 0','T_a = 1','T_a = 1.5',...
'L_y = 10km','L_y = 120km','L_y = 200km','L_y = 400km',...
'\Deltaz = 5m','ref',...
'Orientation','horizontal','NumColumns',9,'fontsize',26)

set(gca,'Box','on','linewidth',2)












