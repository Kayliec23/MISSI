%%% Calculating PV tendency using MITgcm model diagnostics
 
%%% need u,b and udot, bdot (dot = tendency).
%%% Um_Diss - U momentum tendency from Dissipation (Explicit part)
%%% Um_Impl - U momentum tendency from Dissipation (Implicit part)
%%% Um_Advec - U momentum tendency from Advection terms
%%% Um_Cori - U momentum tendency from Coriolis term
%%% Um_dPhiX - U momentum tendency from Pressure/Potential grad
%%% USidDrag - U momentum tendency from Side Drag
%%% ADVx_TH - Zonal Advective Flux of Pot.Temperature
%%% DFxE_TH - Zonal Diffusive Flux of Pot.Temperature
%%% ADVx_SLT - Zonal Advective Flux of Salinity
%%% DFxE_SLT - Zonal Diffusive Flux of Salinity
%%% TOTUTEND - Tendency of Zonal Component of Velocity
%%% TOTTTEND - Tendency of Potential Temperature
%%% TOTSTEND - Tendency of Salinity
 
addpath /Users/Kayliec23/Desktop/analysis
addpath /Users/Kayliec23/Desktop/MISSI/
 
% Constants
Nx = length(delX);
Ny = length(delY);
Nr = length(delR);
Nt = size(T,4);
dy = delY(1);
dx = delX(1);
dz = delR(1);
Omega = 2*pi*366/365/86400;
lat0 = -80; %%% Latitude at southern boundary
f0 = 2*Omega*sind(lat0); %%% Coriolis parameter
m1km = 1000;
g = 9.81;
rho0 = 1027;
t1min = 60;
t1hour = 60*t1min;
t1day = 24*t1hour;
dt = t1day;
[YY,ZZ] = meshgrid(yy,zz);
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% define b and its components
% find alpha and beta
SS = iceNaN(S,Nx,Ny,Nr,icetopo,zz); SS = squeeze(mean(SS(1,:,:,:),4,'omitnan'));
TT = iceNaN(T,Nx,Ny,Nr,icetopo,zz); TT = squeeze(mean(TT(1,:,:,:),4,'omitnan'));
dPT = 1e-3;
dSP = 1e-3;
alpha = - (densjmd95(SS,TT+0.5*dPT,0) - densjmd95(SS,TT-0.5*dPT,0))./densjmd95(SS,TT,0)/dPT;
beta = (densjmd95(SS+0.5*dSP,TT,0) - densjmd95(SS-0.5*dSP,TT,0))./densjmd95(SS,TT,0)/dSP;
 
%%% calculate b
for n = 1:Nt
    b_(1,:,:,n) = g.*(alpha.*squeeze(T(1,:,:,n)) - beta.*squeeze(S(1,:,:,n))) - g.*(alpha.*squeeze(T(1,:,:,1)) - beta.*squeeze(S(1,:,:,1)));
end
%b = (b(:,:,:,1:Nt-1)+b(:,:,:,2:Nt))./2;
b = iceNaN(b_,Nx,Ny,Nr,icetopo,zz);
b = (b_(:,:,:,1:Nt-1)+b_(:,:,:,2:Nt))/2;
bt = diff(b,1,4)./dt;
 
% calculate bdot_tot
[bdot bdotdot] = buoyancy(TOTSTEND./t1day,TOTTTEND./t1day,alpha,beta,icetopo,zz,dt);
bdot = iceNaN(bdot,Nx,Ny,Nr,icetopo,zz);
 
vol = dx.*dy.*dz; %Lx.*Ly.*flip(-icetopo);
% calculate bAdv
advY_S(:,1:Ny-1,:,:) = diff(ADVy_S,1,2); advY_S(:,Ny,:,:) = 2*advY_S(:,Ny-1,:,:) - advY_S(:,Ny-2,:,:); advY_S = advY_S./vol;
advY_T(:,1:Ny-1,:,:) = diff(ADVy_T,1,2); advY_T(:,Ny,:,:) = 2*advY_T(:,Ny-1,:,:) - advY_T(:,Ny-2,:,:); advY_T = advY_T./vol;
advR_S(:,:,1:Nr-1,:) = (ADVr_S(:,:,1:Nr-1,:) - ADVr_S(:,:,2:Nr,:))./vol; advR_S(:,:,Nr,:) = 2*advR_S(:,:,Nr-1,:) - advR_S(:,:,Nr-2,:);
%advR_S(:,:,2:Nr,:) = (advR_S(:,:,2:Nr,:) + advR_S(:,:,1:Nr-1,:))./2; advR_S(:,:,1,:) = 2*advR_S(:,:,2,:) - advR_S(:,:,3,:);
advR_T(:,:,1:Nr-1,:) = (ADVr_T(:,:,1:Nr-1,:) - ADVr_T(:,:,2:Nr,:))./vol; advR_T(:,:,Nr,:) = 2*advR_T(:,:,Nr-1,:) - advR_T(:,:,Nr-2,:);
%advR_T(:,:,2:Nr,:) = (advR_T(:,:,2:Nr,:) + advR_T(:,:,1:Nr-1,:))./2; advR_T(:,:,1,:) = 2*advR_T(:,:,2,:) - advR_T(:,:,3,:);
 
[bAdvy bAdvydot] = buoyancy(advY_S,advY_T,alpha,beta,icetopo,zz,dt);
[bAdvz bAdvzdot] = buoyancy(advR_S,advR_T,alpha,beta,icetopo,zz,dt);
bAdv = bAdvy + bAdvz;
bAdv = iceNaN(bAdv,Nx,Ny,Nr,icetopo,zz);
 
% calculate bDiff (horizontal component equals zero in set-up)
diffY_S(:,1:Ny-1,:,:) = diff(DFyE_S,1,2); diffY_S(:,Ny,:,:) = 2*diffY_S(:,Ny-1,:,:) - diffY_S(:,Ny-2,:,:); diffY_S = diffY_S./vol;
diffY_T(:,1:Ny-1,:,:) = diff(DFyE_T,1,2); diffY_T(:,Ny,:,:) = 2*diffY_T(:,Ny-1,:,:) - diffY_T(:,Ny-2,:,:); diffY_T = diffY_T./vol;
DFr_T = DFrI_T + DFrE_T;
DFr_S = DFrI_S + DFrE_S;
diffR_S(:,:,1:Nr-1,:) = (DFr_S(:,:,1:Nr-1,:) - DFr_S(:,:,2:Nr,:))./vol; diffR_S(:,:,Nr,:) = 2*diffR_S(:,:,Nr-1,:) - diffR_S(:,:,Nr-2,:);
%diffR_S(:,:,2:Nr,:) = (diffR_S(:,:,2:Nr,:) + diffR_S(:,:,1:Nr-1,:))./2; diffR_S(:,:,1,:) = 2*diffR_S(:,:,2,:) - diffR_S(:,:,3,:);
diffR_T(:,:,1:Nr-1,:) = (DFr_T(:,:,1:Nr-1,:) - DFr_T(:,:,2:Nr,:))./vol; diffR_T(:,:,Nr,:) = 2*diffR_T(:,:,Nr-1,:) - diffR_T(:,:,Nr-2,:);
%diffR_T(:,:,2:Nr,:) = (diffR_T(:,:,2:Nr,:) + diffR_T(:,:,1:Nr-1,:))./2; diffR_T(:,:,1,:) = 2*diffR_T(:,:,2,:) - diffR_T(:,:,3,:);
 
[bdiffy bdiffydot] = buoyancy(diffY_S,diffY_T,alpha,beta,icetopo,zz,dt);
[bdiffz bdiffzdot] = buoyancy(diffR_S,diffR_T,alpha,beta,icetopo,zz,dt);
bDiff = bdiffy + bdiffz;
bDiff = iceNaN(bDiff,Nx,Ny,Nr,icetopo,zz);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% calculate PV flux-form terms %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% define u and its components
u = (U(:,:,:,1:Nt-1)+U(:,:,:,2:Nt))/2;
 
%%% define udot
udot = TOTUTEND./t1day;
udot = iceNaN(udot,Nx,Ny,Nr,icetopo,zz);
 
% calculate uvisc
% viscY(:,1:Ny-1,:,:) = diff(VISCy_Um,1,2); viscY(:,Ny,:,:) = 2*viscY(:,Ny-1,:,:) - viscY(:,Ny-2,:,:); viscY = viscY./vol;
% VISr = VISrI_Um + VISrE_Um;
% viscR(:,:,1:Nr-1,:) = (VISr(:,:,1:Nr-1,:) - VISr(:,:,2:Nr,:))./vol; viscR(:,:,Nr,:) = 2*viscR(:,:,Nr-1,:) - viscR(:,:,Nr-2,:);
% uvisc = viscY + viscR;
 
% calculate udiss
udiss = Um_Diss + Um_ImplD;
 
qdot = PVtend(udot,bdot,u,b,f0,dz,dy);
qAdv = PVtend(-Um_Advec,bAdv,u,b,f0,dz,dy);
qF = PVtend(-udiss,zeros(size(udiss)),u,b,f0,dz,dy);
qB = PVtend(zeros(size(bDiff)),bDiff,u,b,f0,dz,dy);
qtot = PVtend(udot - (Um_Advec + udiss),bdot + bAdv + bDiff,u,b,f0,dz,dy);
 
%%% Calculate Q and take time derivative
q = calculateQ(S,T,U,V,W,Nx,Ny,Nr,Nt,dy,dz,g,rho0,f0);
q = iceNaN(q,Nx,Ny,Nr,icetopo,zz);
qt = diff(q,1,4)./dt;
 
%%% tendency residual
qt_residual = qt - qdot;
 
% calculate cumulative sum anomaly
qdot_sum = cumsum(qdot,4).*dt; qdot_sum_anom = qdot_sum - qdot(:,:,:,1);
qt_sum = cumsum(qt,4).*dt; qt_sum_anom = qt_sum - qt(:,:,:,1);
q_sum = cumsum(q,4).*dt; q_sum_anom = q_sum - q(:,:,:,1);
qB_sum = cumsum(qB,4).*dt; qB_sum_anom = qB_sum - qB(:,:,:,1);
qF_sum = cumsum(qF,4).*dt; qF_sum_anom = qF_sum - qF(:,:,:,1);
qAdv_sum = cumsum(qAdv,4).*dt; qAdv_sum_anom = qAdv_sum - qAdv(:,:,:,1);
qt_residual_sum = cumsum(qt-qdot,4).*dt;  qt_residual_sum_anom = qt_residual_sum - qt_residual(:,:,:,1);
% calculate cumulative sum anomaly - positive regions only
qB_sum_anom_pos = posVal(qB_sum_anom);
qF_sum_anom_pos = posVal(qF_sum_anom);
qt_sum_anom_pos = posVal(qt_sum_anom);
qAdv_sum_anom_pos = posVal(qAdv_sum_anom);
% calculate cumulative sum anomaly - regions of positive qt_sum_anom only
qB_sum_anom_posPV = posPVval(qB_sum_anom,qt_sum_anom);
qF_sum_anom_posPV = posPVval(qF_sum_anom,qt_sum_anom);
qt_sum_anom_posPV = posPVval(qt_sum_anom,qt_sum_anom);
qAdv_sum_anom_posPV = posPVval(qAdv_sum_anom,qt_sum_anom);
 
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures fq components %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /Users/Kayliec23/Desktop/plottingKit/
 
fontsize = 18;
scrsz = get(0,'ScreenSize');
fontsize = 28;
framepos = [0 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(30);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');
t1 = 2; t2 = 8; t3 = 60;
gap = [0.03 0.04];
marg_h = [0.065 0.075];
marg_w = [0.075 0.1];
hplot = 0.27;
wplot = 0.14;
clim = [-1 1]*1e-6;
slope_exp = 'linear';
switch slope_exp
    case 'linear'
        x_lim = [15 25];
        y_lim = [-330 -200];
    case 'tanh'
        x_lim =  [10 20]; %[0 7];
        y_lim = [-250 -90]; %[-500 -300];
end
 
curve1 = icetopo;
curve2 = zeros(1,length(curve1));
x1 = yy/m1km; %linspace(xlim_min, xlim_max,length(curve1));
x2 = [x1, fliplr(x1)];
inBetween = [curve1, fliplr(curve2)];
c_ice = [0.9 0.9 0.9];
 
% qt
h1 = subplot('position',[0.08   0.6874    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(qt_sum_anom(:,:,:,t1))'); shading flat;
mymap = colormap('redblue');
colormap(h1,mymap)
caxis([-5 5]*1e-10)
hold on;
ylabel('z (m)')
xticklabels([]);
yticks([-310 -280 -250 -220])
title('q_t')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(a)']};
h_annot=annotation('textbox',...
    [0.082  0.84 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','tex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
txt={['t=',num2str(t1),'day']};
h_annot=annotation('textbox',...
    [0.14  0.63 1 0.11],...
    'String',txt,...
    'FontSize',26,...
    'FontName','latex',...
    'EdgeColor','none',...
    'LineWidth',1.5,...
    'BackgroundColor','none',...
    'FitBoxToText','on');
set(gca,'fontsize',fontsize)
 
h2 = subplot('position',[0.08    0.4034    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(qt_sum_anom(:,:,:,t2))'); shading flat;
mymap = colormap('redblue');
colormap(h2,mymap);
cb = colorbar('fontsize',20,'box','on','linewidth',2,'location','eastoutside','position',[0.23    0.415    0.01   hplot*2-0.015]);
caxis([-5 5]*1e-10)
hold on;
set(gca,'fontsize',fontsize)
xticklabels([]);
yticks([-310 -280 -250 -220])
ylabel('z (m)')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(e)']};
h_annot=annotation('textbox',...
    [0.082  0.56 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
txt={['t=',num2str(t2),'day']};
h_annot=annotation('textbox',...
    [0.14  0.34 1 0.11],...
    'String',txt,...
    'FontSize',26,...
    'FontName','latex',...
    'EdgeColor','none',...
    'LineWidth',1.5,...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
 
h3 = subplot('position',[0.08    0.1195    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(qt_sum_anom(:,:,:,t3))'); shading flat;
mymap = colormap('redblue');
colormap(h3,mymap);
cb = colorbar('fontsize',20,'box','on','linewidth',2,'location','eastoutside','position',[0.23    0.13    0.01   hplot-0.025]);
caxis([-5 5]*1e-9)
hold on;
set(gca,'fontsize',fontsize)
ylabel('z (m)');
xlabel('y (km)');
xlim(x_lim);
ylim(y_lim);
xticks([17 22])
yticks([-310 -280 -250 -220])
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(i)']};
h_annot=annotation('textbox',...
    [0.082  0.27 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
txt={['t=',num2str(t3),'day']};
h_annot=annotation('textbox',...
    [0.13  0.05 1 0.11],...
    'String',txt,...
    'FontSize',26,...
    'FontName','tex',...
    'EdgeColor','none',...
    'LineWidth',1.5,...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
h4 = subplot('position',[0.29    0.6874    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(-qAdv_sum_anom(:,:,:,t1))'); shading flat;
%colorbar('FontSize',20,'FontName','latex','linewidth',1.5,'position',[0.93    0.7350    0.0116    0.19000]);
mymap = colormap('redblue');
colormap(h4,mymap)
caxis([-5 5]*1e-10)
hold on;
set(gca,'fontsize',fontsize)
title('q_{adv}')
xticklabels([]);
yticklabels([]);
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(b)']};
h_annot=annotation('textbox',...
    [0.295  0.84 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
h5 = subplot('position',[0.29    0.4034    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(-qAdv_sum_anom(:,:,:,t2))'); shading flat;
mymap = colormap('redblue');
colormap(h5,mymap)
caxis([-5 5]*1e-7)
hold on;
set(gca,'fontsize',fontsize)
xticklabels([]);
yticklabels([])
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(f)']};
h_annot=annotation('textbox',...
    [0.295  0.56 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
h6 = subplot('position',[0.29    0.1195    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(-qAdv_sum_anom(:,:,:,t3))'); shading flat;
mymap = colormap('redblue');
colormap(h6,mymap)
caxis([-5 5]*1e-6)
hold on;
set(gca,'fontsize',fontsize)
yticklabels([])
xlabel('y (km)')
xticks([17 22])
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(j)']};
h_annot=annotation('textbox',...
    [0.295  0.27 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
h7 = subplot('position',[0.445    0.6874    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(-qF_sum_anom(:,:,:,t1))'); shading flat;
mymap = colormap('redblue');
colormap(h7,mymap)
caxis([-5 5]*1e-10)
hold on;
set(gca,'fontsize',fontsize)
xticklabels([]);
yticklabels([]);
title('qF')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(c)']};
h_annot=annotation('textbox',...
    [0.45  0.84 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
 
h8 = subplot('position',[0.445    0.4034    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(-qF_sum_anom(:,:,:,t2))'); shading flat;
%colorbar('FontSize',20,'FontName','latex','linewidth',1.5,'position',[0.93    0.074    0.0116    0.6204]);
mymap = colormap('redblue');
colormap(h8,mymap)
caxis([-5 5]*1e-7)
hold on;
set(gca,'fontsize',fontsize)
xticklabels([]);
yticklabels([]);
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(g)']};
h_annot=annotation('textbox',...
    [0.45  0.56 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
h9 = subplot('position',[0.445    0.1195    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(-qF_sum_anom(:,:,:,t3))'); shading flat;
mymap = colormap('redblue');
colormap(h9,mymap)
caxis([-5 5]*1e-6)
hold on;
set(gca,'fontsize',fontsize)
xlabel('y (km)')
xlim(x_lim);
ylim(y_lim);
xticks([17 22])
yticklabels([]);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(k)']};
h_annot=annotation('textbox',...
    [0.45  0.27 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
 
 
h10 = subplot('position',[0.601    0.6874    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(-qB_sum_anom(:,:,:,t1))'); shading flat;
mymap = colormap('redblue');
colormap(h10,mymap)
caxis([-5 5]*1e-10)
cb = colorbar('fontsize',20,'box','on','linewidth',2,'location','eastoutside','position',[0.75    0.7    0.01   hplot-0.025]);
hold on;
set(gca,'fontsize',fontsize);
title('qB')
yticklabels([]);
xticklabels([])
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(d)']};
h_annot=annotation('textbox',...
    [0.606  0.84 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
 
h11 = subplot('position',[0.601    0.4034    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(-qB_sum_anom(:,:,:,t2))'); shading flat;
mymap = colormap('redblue');
colormap(h11,mymap)
caxis([-5 5]*1e-7);
cb = colorbar('fontsize',20,'box','on','linewidth',2,'location','eastoutside','position',[0.75    0.415    0.01   hplot-0.025]);
hold on;
set(gca,'fontsize',fontsize)
yticklabels([]);
xticklabels([]);
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(h)']};
h_annot=annotation('textbox',...
    [0.606  0.56 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
 
h12 = subplot('position',[0.601    0.1195   wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(-qB_sum_anom(:,:,:,t3))'); shading flat;
mymap = colormap('redblue');
colormap(h12,mymap)
caxis([-5 5]*1e-6);
cb = colorbar('fontsize',20,'box','on','linewidth',2,'location','eastoutside','position',[0.75    0.13    0.01   hplot-0.025]);
hold on;
set(gca,'fontsize',fontsize)
yticklabels([]);
xticks([17 22])
xlabel('y (km)')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(l)']};
h_annot=annotation('textbox',...
    [0.606  0.27 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures PV closure %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
addpath /Users/Kayliec23/Desktop/plottingKit/
 
fontsize = 18;
scrsz = get(0,'ScreenSize');
fontsize = 28;
framepos = [0 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(30);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');
 
curve1 = icetopo;
curve2 = zeros(1,length(curve1));
x1 = yy/m1km; %linspace(xlim_min, xlim_max,length(curve1));
x2 = [x1, fliplr(x1)];
inBetween = [curve1, fliplr(curve2)];
c_ice = [0.9 0.9 0.9];
 
h1 = subplot('position',[0.1300    0.55    0.24    0.4]);
pcolor(YY/m1km,ZZ,squeeze(qAdv_sum_anom(:,:,:,end))'); shading flat;
mymap = colormap('redblue');
colormap(h1,mymap)
caxis([-1 1]*1e-5);
hold on;
ylabel('z (m)')
xticklabels([]);
title('$\mathrm{q_{adv}}$','interpreter','latex')
fill(x2, inBetween, c_ice,'EdgeColor','none');
xlim([0 40])
set(gca,'fontsize',24)
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(a)']};
h_annot=annotation('textbox',...
    [0.14  0.81 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','tex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
h2 = subplot('position',[0.4    0.55    0.24    0.4]);
pcolor(YY/m1km,ZZ,squeeze(qF_sum_anom(:,:,:,end))'); shading flat;
mymap = colormap('redblue');
colormap(h2,mymap)
caxis([-1 1]*1e-5)
hold on;
xticklabels([]);
yticklabels([])
title('$\mathrm{qF}$','interpreter','latex')
fill(x2, inBetween, c_ice,'EdgeColor','none');
xlim([0 40])
set(gca,'fontsize',24)
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(b)']};
h_annot=annotation('textbox',...
    [0.41  0.81 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','tex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
h3 = subplot('position',[0.67    0.55    0.24    0.4]);
pcolor(YY/m1km,ZZ,squeeze(qB_sum_anom(:,:,:,end))'); shading flat;
mymap = colormap('redblue');
colormap(h3,mymap)
caxis([-1 1]*1e-5);
colorbar('position',[0.93    0.55    0.015    0.4],'linewidth',1.5);
hold on;
yticklabels([]);
xticklabels([]);
title('$\mathrm{qB}$','interpreter','latex')
xlim([0 40])
fill(x2, inBetween, c_ice,'EdgeColor','none');
set(gca,'fontsize',24)
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(c)']};
h_annot=annotation('textbox',...
    [0.68  0.81 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','tex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
h4 = subplot('position',[0.13    0.08    0.24    0.4]);
pcolor(YY/m1km,ZZ,squeeze(qdot_sum_anom(:,:,:,end))'); shading flat;
mymap = colormap('redblue');
colormap(h4,mymap)
caxis([-5 5]*1e-7)
hold on;
ylabel('z (m)')
xlabel('y (km)')
set(gca,'fontsize',24)
title('$\mathrm{\dot{q}}$','interpreter','latex')
xlim([0 40])
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(d)']};
h_annot=annotation('textbox',...
    [0.14  0.36 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','tex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
h5 = subplot('position',[0.4    0.08    0.24    0.4]);
pcolor(YY/m1km,ZZ,squeeze(qt_residual_sum(:,:,:,end))'); shading flat;
mymap = colormap('redblue');
colormap(h5,mymap)
caxis([-5 5]*1e-7)
hold on;
xlabel('y (km)')
yticklabels([]);
set(gca,'fontsize',24)
title('$\mathrm{q_t^{residual}}$','interpreter','latex')
xlim([0 40])
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(e)']};
h_annot=annotation('textbox',...
    [0.41  0.36 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','tex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
 
h6 = subplot('position',[0.67    0.08    0.24    0.4]);
pcolor(YY/m1km,ZZ,squeeze(qtot_sum_anom(:,:,:,end))'); shading flat;
mymap = colormap('redblue');
colormap(h6,mymap);
colorbar('position',[0.93 0.08 0.015 0.4],'linewidth',2);
caxis([-5 5]*1e-7);
hold on;
yticklabels([]);
xlabel('y (km)')
xlim([0 40])
title('$\mathrm{q_{tot}}$','interpreter','latex')
fill(x2, inBetween, c_ice,'EdgeColor','none');
set(gca,'fontsize',24)
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(f)']};
h_annot=annotation('textbox',...
    [0.68  0.36 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','tex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [b bdot] = buoyancy(S,T,alpha,beta,icetopo,zz,dt)
rho0 = 1027;
g = 9.81;
Nx = size(S,1);
Ny = size(S,2);
Nr = size(S,3);
Nt = size(S,4);
for n = 1:Nt
    b(1,:,:,n) = g.*(alpha.*squeeze(T(1,:,:,n))-beta.*squeeze(S(1,:,:,n)))- g.*(alpha.*squeeze(T(1,:,:,1)) - beta.*squeeze(S(1,:,:,1)));
end
b = iceNaN(b,Nx,Ny,Nr,icetopo,zz);
bdot = (b(:,:,:,2:Nt)-b(:,:,:,1:Nt-1)) ./ dt;
end
 
function BB = interp_yz(AA)
Nx = size(AA,1);
Ny = size(AA,2);
Nr = size(AA,3);
Nt = size(AA,4);
BB = NaN(Nx,Ny+1,Nr+1,Nt);
BB(:,1:Ny,2:Nr+1,:) = AA;
BB(:,:,1,:) = 2*BB(:,:,2,:) - BB(:,:,3,:);   % interpolate in z for boundary point
BB(:,Ny+1,:,:) = 2*BB(:,Ny,:,:) - BB(:,Ny-1,:,:);    % interpolate in y for boundary point
end
 
 
function Q = calculateQ(S,T,U,V,W,Nx,Ny,Nr,Nt,dy,dz,g,rho0,f0)
% Calculate potential density (at grid center)
sigmatheta = zeros(Nx,Ny,Nr,Nt);
for n = 1:Nt
    SS = squeeze(S(1,:,:,n));
    TT = squeeze(T(1,:,:,n));
    theta(:,:) = densjmd95(SS,TT,0);
    sigmatheta(1,:,:,n) = theta - 1000;
end
 
 
%     % calculate ZETAx
%     dWdy = NaN(Nx,Ny,Nr,Nt); dWdy(:,1:Ny-1,:,:) = diff(W(:,:,:,1:Nt),1,2)./dy; % derivative wrt y, moves to bottom right grid
%     dWdy(:,Ny,:,:) = 2*dWdy(:,Ny-1,:,:) - dWdy(:,Ny-2,:,:); % add boundary at the end
%     dVdz = NaN(Nx,Ny,Nr,Nt); dVdz(:,:,2:Nr,:) = diff(V(:,:,:,1:Nt),1,3)./dz; % derivative wrt z, moves to bottom left grid, shifted to bottom right
%     dVdz(:,:,1,:) = 2*dVdz(:,:,2,:) - dVdz(:,:,3,:); % add boundary at the beginning
%     ZETAx = dWdy - dVdz;
 
% Calculate ZETAy
dUdz = NaN(Nx,Ny,Nr,Nt); dUdz(:,:,2:Nr,:) = (U(:,:,[1:Nr-1],1:Nt) - U(:,:,[2:Nr],1:Nt))./dz; % derivative wrt z, shifts to top face center
dUdz(:,1:Ny-1,:,:) = (dUdz(:,[2:Ny],:,:)+dUdz(:,[1:Ny-1],:,:))/2; % average to top right corner
dUdz(:,:,1,:) = 2*dUdz(:,:,2,:)-dUdz(:,:,3,:); dUdz(:,Ny,:,:) = 2*dUdz(:,Ny-1,:,:)-dUdz(:,Ny-2,:,:); % add boundaries
ZETAy = dUdz;
 
% Calcualte ZETAz
dUdy = NaN(Nx,Ny,Nr,Nt); dUdy(:,1:Ny-1,:,:) = diff(U(:,:,:,1:Nt),1,2)./dy; % derivative wrt y, shifts to right face center
dUdy(:,:,2:Nr,:) = (dUdy(:,:,[2:Nr],:)+dUdy(:,:,[1:Nr-1],:))/2; % average to top right corner
dUdy(:,:,1,:) = 2*dUdy(:,:,2,:)-dUdy(:,:,3,:); dUdy(:,Ny,:,:) = 2*dUdy(:,Ny-1,:,:)-dUdy(:,Ny-2,:,:); % add boundaries
ZETAz = -dUdy;
 
% Density derivatives
dRHOdy = NaN(Nx,Ny,Nr,Nt); dRHOdy(:,1:Ny-1,:,:) = diff(sigmatheta,1,2)./dy; % derivative wrt y, shifts to right face center
dRHOdy(:,:,2:Nr,:) = (dRHOdy(:,:,[2:Nr],:)+dRHOdy(:,:,[1:Nr-1],:))/2; % move to top right corner
dRHOdy(:,Ny,:,:) = 2*dRHOdy(:,Ny-1,:,:)-dRHOdy(:,Ny-2,:,:); dRHOdy(:,:,1,:) = 2*dRHOdy(:,:,2,:)-dRHOdy(:,:,3,:); % add boundaries
dRHOdz = NaN(Nx,Ny,Nr,Nt); dRHOdz(:,:,2:Nr,:) = (sigmatheta(:,:,[1:Nr-1],1:Nt) - sigmatheta(:,:,[2:Nr],1:Nt))./dz; % derivative wrt z, shifts to top face center
dRHOdz(:,1:Ny-1,:,:) = (dRHOdz(:,[2:Ny],:,:)+dRHOdz(:,[1:Ny-1],:,:))/2; % move to top right corner
dRHOdz(:,:,1,:) = 2*dRHOdz(:,:,2,:)-dRHOdz(:,:,3,:); dRHOdz(:,Ny,:,:) = 2*dRHOdz(:,Ny-1,:,:)-dRHOdz(:,Ny-2,:,:); % add boundaries
 
dbdy = -g/rho0.*dRHOdy;
dbdz = -g/rho0.*dRHOdz;
 
 
% Q and its symmetric, convecitve and centrifugal components
%Qsym = dbdy.*ZETAy; % symmetric
%Qconv = dbdz.*f0; % convective
%Qcent = dbdz.*ZETAz; % centrifugal
Q = dbdy.*ZETAy + (f0 + ZETAz).*dbdz; % total
 
end
 
 
function Apos = posVal(A)
Ny = size(A,2);
Nr = size(A,3);
Nt = size(A,4);
for j = 1:Ny
    for k = 1:Nr
        for n = 1:Nt
            if A(1,j,k,n)>0
                Apos(1,j,k,n) = A(1,j,k,n);
            else
                Apos(1,j,k,n) = NaN;
            end
        end
    end
end
end
 
function Apos = posPVval(A,qt)
Ny = size(A,2);
Nr = size(A,3);
Nt = size(A,4);
for j = 1:Ny
    for k = 1:Nr
        for n = 1:Nt
            if qt(1,j,k,n)>0
                Apos(1,j,k,n) = A(1,j,k,n);
            else
                Apos(1,j,k,n) = NaN;
            end
        end
    end
end
end
 
function AA = PVtend(udot_comp,bdot_comp,u,b,f0,dz,dy)
Ny = size(u,2);
Nr = size(u,3);
Nt = size(u,4);
 
%%% term1
term1 = f0*(bdot_comp(:,:,[1:Nr-1],:) - bdot_comp(:,:,[2:Nr],:))./dz;    % diff in z
term1 = (term1(:,1:Ny-1,:,:) + term1(:,2:Ny,:,:)) / 2;    % average locally in y
term1 = interp_yz(term1);    % interpolate in z and y
 
%%% term2
A = diff(udot_comp,1,2); % diff in y
A = (A(:,:,1:Nr-1,:) + A(:,:,2:Nr,:)) / 2; % averaging locally in z
A = interp_yz(A);    % interpolate in z and y
B = (b(:,:,[1:Nr-1],:) - b(:,:,[2:Nr],:))./dz;   % diff in z
B = (B(:,1:Ny-1,:,:) + B(:,2:Ny,:,:)) / 2;    % average locally in y
B = interp_yz(B);    % interpolate in z and y
term2 = -(A.*B)./dy./dz;
 
%%% term3
A = diff(u,1,2); % diff in y
A = (A(:,:,1:Nr-1,:) + A(:,:,2:Nr,:)) / 2; % averaging locally in z
A = interp_yz(A);    % interpolate in z and y
B = (bdot_comp(:,:,[1:Nr-1],:) - bdot_comp(:,:,[2:Nr],:))./dz;    % diff in z
B = (B(:,1:Ny-1,:,:) + B(:,2:Ny,:,:)) / 2;    % average locally in y
B = interp_yz(B);    % interpolate in z and y
term3 = -(A.*B)./dy./dz;
 
%%% term4
A = (udot_comp(:,:,[1:Nr-1],:) - udot_comp(:,:,[2:Nr],:))./dz; % diff in z
A = (A(:,1:Ny-1,:,:) + A(:,2:Ny,:,:)) / 2;    % average locally in y
A = interp_yz(A);    % interpolate in z and y
B = diff(b,1,2); % diff in y
B = (B(:,:,1:Nr-1,:) + B(:,:,2:Nr,:)) / 2;    % average locally in z
B = interp_yz(B);    % interpolate in z and y
term4 = (A.*B)./dy./dz;
 
%%% term5
A = (u(:,:,[1:Nr-1],:) - u(:,:,[2:Nr],:))./dz;    % diff in z
A = (A(:,1:Ny-1,:,:) + A(:,2:Ny,:,:)) / 2;    % average locally in y
A = interp_yz(A);    % interpolate in z and y
B = diff(bdot_comp,1,2);    % diff in y
B = (B(:,:,1:Nr-1,:) + B(:,:,2:Nr,:)) / 2;    % average locally in z
B = interp_yz(B);    % interpolate in z and y
term5 = (A.*B)./dy./dz;
 
AA = term1 + term2 + term3 + term4 + term5;
end

