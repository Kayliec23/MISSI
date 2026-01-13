%%% Melt rate figures for paper

nsimtime = 25;
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




subplot(5,1,1); hold on;
%%% for drag 
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/Drag

load('cd2_5e-4/fwFlux_cd2_5e-4.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,1); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0 1 0.25])
subplot(6,1,6); hold on;
scatter(1,mmax,'filled','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.25)
mtot = -sum(m(:));
scatter(1,mtot,'filled','diamond','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.25)

load('cd1_4e-3/fwFlux_cd1_4e-3.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,1); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0 1 0.5])
subplot(6,1,6); hold on;
scatter(2,mmax,'filled','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.5)
mtot = -sum(m(:));
scatter(2,mtot,'filled','diamond','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.5)


load('cd1_4e-2/fwFlux_cd1_4e-2.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,1); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0 1 0.75])
subplot(6,1,6); hold on;
scatter(3,mmax,'filled','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.75)
mtot = -sum(m(:));
scatter(3,mtot,'filled','diamond','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.75)


load('cd2_5e-2/fwFlux_cd2_5e-2.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,1); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0 1 1]);
subplot(6,1,6); hold on;
scatter(4,mmax,'filled','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',1)
mtot = -sum(m(:));
scatter(4,mtot,'filled','diamond','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',1)

% Ref
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/
load('Ref2D/fwFlux_ref2D.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
subplot(5,1,1); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'--k','linewidth',3.5); 


legend('C_d = 2.5e-4','C_d = 1.4e-3','C_d = 1.4e-2','C_d = 2.5e-2','numColumns',4,'fontsize',18)
xlim([0 40])
ylim([0 1.1])
set(gca,'fontsize',24);
xticklabels([]);
set(gca,'Box','on','linewidth',2)
set(gca,'TickLabelInterpreter','latex')


subplot(5,1,2); hold on;
%%% for dS
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/dS

load('dS04/fwFlux_dS04.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,2); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0.5 0 0.25]);
subplot(6,1,6); hold on;
scatter(5,mmax,'filled','MarkerFaceColor',[0 0.5 0],'MarkerFaceAlpha',0.25)


load('dS06/fwFlux_dS06.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,2); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0.5 0 0.5]);
subplot(6,1,6); hold on;
scatter(6,mmax,'filled','MarkerFaceColor',[0 0.5 0],'MarkerFaceAlpha',0.5)


load('dS08/fwFlux_dS08.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,2); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0.5 0 0.75]);
subplot(6,1,6); hold on;
scatter(7,mmax,'filled','MarkerFaceColor',[0 0.5 0],'MarkerFaceAlpha',0.75)


load('dS1/fwFlux_dS1.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,2); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0 0.5 0 1]);
subplot(6,1,6); hold on;
scatter(8,mmax,'filled','MarkerFaceColor',[0 0.5 0],'MarkerFaceAlpha',1)


% Ref
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/
load('Ref2D/fwFlux_ref2D.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
subplot(6,1,2); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'--k','linewidth',3.5); 


legend('\DeltaS = 0.4','\DeltaS = 0.6','\DeltaS = 0.8','\DeltaS = 1','numColumns',4,'fontsize',18)
xlim([0 40])
ylim([0 1.1])
set(gca,'fontsize',24);
xticklabels([]);
set(gca,'Box','on','linewidth',2)
set(gca,'TickLabelInterpreter','latex')



subplot(5,1,3); hold on;
%%% for Tforcing
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/Tf

load('Tfneg1/fwFlux_Tfneg1.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,3); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0.5 0 0.5 0.2]); 
subplot(6,1,6); hold on;
scatter(9,mmax,'filled','MarkerFaceColor',[0.5 0 0.5],'MarkerFaceAlpha',0.2)

load('Tfneg05/fwFlux_Tfneg05.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,3); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0.5 0 0.5 0.4]); 
subplot(6,1,6); hold on;
scatter(10,mmax,'filled','MarkerFaceColor',[0.5 0 0.5],'MarkerFaceAlpha',0.4)

load('Tf0/fwFlux_Tf0.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,3); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0.5 0 0.5 0.6]); 
subplot(6,1,6); hold on;
scatter(11,mmax,'filled','MarkerFaceColor',[0.5 0 0.5],'MarkerFaceAlpha',0.6)

load('Tf1/fwFlux_Tf1.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,3); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0.5 0 0.5 0.8]);
subplot(6,1,6); hold on;
scatter(12,mmax,'filled','MarkerFaceColor',[0.5 0 0.5],'MarkerFaceAlpha',0.8)

load('Tf1.5/fwFlux_Tf1.5.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,3); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[0.5 0 0.5 1]);
subplot(6,1,6); hold on;
scatter(13,mmax,'filled','MarkerFaceColor',[0.5 0 0.5],'MarkerFaceAlpha',1)

% Ref
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/
load('Ref2D/fwFlux_ref2D.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
subplot(6,1,3); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'--k','linewidth',3.5); 


legend('T_a = -1','T_a = -0.5','T_a = 0','T_a = 1','T_a = 1.5','numColumns',5,'fontsize',18)
xlim([0 40])
ylim([0 1.1])
set(gca,'fontsize',24);
xticklabels([]);
ylabel('normalized total melt','interpreter','latex');
set(gca,'Box','on','linewidth',2)
set(gca,'TickLabelInterpreter','latex')


subplot(5,1,4); hold on;
%%% for slope size
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/Slope

load('Ly20/fwFlux_Ly20.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,4); hold on;
plot(yy(1:100)/m1km,-smooth(m(1,1:100)/max(abs(m)),7),'linewidth',2,'color',[1 0.6445 0 0.25]);  
subplot(6,1,6); hold on;
scatter(14,mmax,'filled','MarkerFaceColor',[1 0.6445 0],'MarkerFaceAlpha',0.25)

load('Ly80/fwFlux_Ly80.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,4); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[1 0.6445 0 0.5]);
subplot(6,1,6); hold on;
scatter(15,mmax,'filled','MarkerFaceColor',[1 0.6445 0],'MarkerFaceAlpha',0.5)

load('Ly200/fwFlux_Ly200.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,4); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[1 0.6445 0 0.75]);
subplot(6,1,6); hold on;
scatter(16,mmax,'filled','MarkerFaceColor',[1 0.6445 0],'MarkerFaceAlpha',0.75)

load('Ly400/fwFlux_Ly400.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
mmax = max(abs(m));
subplot(6,1,4); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[1 0.6445 0 1]);
subplot(6,1,6); hold on;
scatter(17,mmax,'filled','MarkerFaceColor',[1 0.6445 0],'MarkerFaceAlpha',1)

% Ref
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/
load('Ref2D/fwFlux_ref2D.mat');
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
subplot(6,1,4); hold on;
plot(yy/m1km,-smooth(m/max(abs(m)),7),'--k','linewidth',3.5); 


legend('L_y = 20km','L_y = 80km','L_y = 200km','L_y = 400km','numColumns',4,'fontsize',18)
xlim([0 40])
ylim([0 1.1])
set(gca,'fontsize',24);
xticklabels([]);
set(gca,'Box','on','linewidth',2)
set(gca,'TickLabelInterpreter','latex')


subplot(5,1,5); hold on;
% dz =5
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/altCases/
load('dy25_dz5/fwFlux_dy25_dz5.mat'); 
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[1 0 0 0.3]); 

% dy = 12.5
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/altCases/
load('Ref2D_dy12_5/fwFlux_Ref2D_dy12_5.mat'); 
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[1 0 0 0.65]); 

% dy = 100
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/altCases/
load('Ref2D_dy100/fwFlux_Ref2D_dy100.mat'); 
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'linewidth',2,'color',[1 0 0 1]); 


% Ref
cd /Users/kalyliec23/Desktop/Research/MISSI/experiments/
load('Ref2D/fwFlux_ref2D.mat');
nt = 30;
m = squeeze(sum(fwFlux(1,:,1:nsimtime),3)*t1day/rhoShelfIce);
plot(yy/m1km,-smooth(m/max(abs(m)),7),'--k','linewidth',3.5); 

xlim([0 40])
ylim([0 1.1])
set(gca,'fontsize',24);
xlabel('y (km)');
legend('\Deltaz = 5m','\Deltay=12.5','\Deltay=100m','numColumns',3,'fontsize',18)
set(gca,'TickLabelInterpreter','latex')
set(gca,'Box','on','linewidth',2)













