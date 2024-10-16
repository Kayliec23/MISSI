 

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
viscY(:,1:Ny-1,:,:) = diff(VISCy_Um,1,2); viscY(:,Ny,:,:) = 2*viscY(:,Ny-1,:,:) - viscY(:,Ny-2,:,:); viscY = viscY./vol;
VISr = VISrI_Um + VISrE_Um;
viscR(:,:,1:Nr-1,:) = (VISr(:,:,1:Nr-1,:) - VISr(:,:,2:Nr,:))./vol; viscR(:,:,Nr,:) = 2*viscR(:,:,Nr-1,:) - viscR(:,:,Nr-2,:);
uvisc = viscY + viscR;

udiss = Um_Diss + Um_ImplD;


qdot = PVtend(udot,bdot,u,b,f0,dz,dy);
qAdv = PVtend(-Um_Advec,bAdv,u,b,f0,dz,dy);
qB = PVtend(-udiss,zeros(size(udiss)),u,b,f0,dz,dy);
qF = PVtend(zeros(size(bDiff)),bDiff,u,b,f0,dz,dy);
qtot = PVtend(udot - (Um_Advec + udiss),bdot + bAdv + bDiff,u,b,f0,dz,dy);

%%% Calculate Q and take time derivative
q = calculateQ(S,T,U,V,W,Nx,Ny,Nr,Nt,dy,dz,g,rho0,f0);
q = iceNaN(q,Nx,Ny,Nr,icetopo,zz);
qt = diff(q,1,4)./dt;

% calculate cumulative sum anomaly
qdot_sum = cumsum(qdot,4).*dt; qdot_sum_anom = qdot_sum - qdot(:,:,:,1);
qt_sum = cumsum(qt,4).*dt; qt_sum_anom = qt_sum - qt(:,:,:,1);
q_sum = cumsum(q,4).*dt; q_sum_anom = q_sum - q(:,:,:,1);
qB_sum = cumsum(qB,4).*dt; qB_sum_anom = qB_sum - qB(:,:,:,1);
qF_sum = cumsum(qF,4).*dt; qF_sum_anom = qF_sum - qF(:,:,:,1);
qAdv_sum = cumsum(qAdv,4).*dt; qAdv_sum_anom = qAdv_sum - qAdv(:,:,:,1);
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

%%% Calculate Q and take time derivative
q = calculateQ(S,T,U,V,W,Nx,Ny,Nr,Nt,dy,dz,g,rho0,f0);
q = iceNaN(q,Nx,Ny,Nr,icetopo,zz);
qt = diff(q,1,4)./dt;

% calculate cumulative sum anomaly
qdot_sum = cumsum(qdot,4).*dt; qdot_sum_anom = qdot_sum - qdot(:,:,:,1);
qt_sum = cumsum(qt,4).*dt; qt_sum_anom = qt_sum - qt(:,:,:,1);
q_sum = cumsum(q,4).*dt; q_sum_anom = q_sum - q(:,:,:,1);
qB_sum = cumsum(qB,4).*dt; qB_sum_anom = qB_sum - qB(:,:,:,1);
qF_sum = cumsum(qF,4).*dt; qF_sum_anom = qF_sum - qF(:,:,:,1);
qAdv_sum = cumsum(qAdv,4).*dt; qAdv_sum_anom = qAdv_sum - qAdv(:,:,:,1);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures fq components %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /Users/Kayliec23/Desktop/plottingKit/

scrsz = get(0,'ScreenSize');
fontsize = 22;

framepos = [0 scrsz(4)/2 scrsz(3) scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(30);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');
t1 = 2; t2 = 4; t3 = 7;
gap = [0.03 0.04];
marg_h = [0.065 0.075];
marg_w = [0.075 0.1];
clim = [-1 1]*1e-10;
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
h1 = subtightplot(4,3,1,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(qt_sum_anom(:,:,:,t1))'); shading flat;
mymap = colormap('redblue');
colormap(h1,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
ylabel('z (m)','interpreter','latex')
xticklabels([])
%title('$q_t$','interpreter','latex')
title(['t = ',num2str(t1),' day'],'interpreter','latex','fontsize',28)
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
ax.FontName = 'latex';
txt={['(a)']};
h_annot=annotation('textbox',...
    [0.078  0.81 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
txt={['q_t']};
h_annot=annotation('textbox',...
    [0.28  0.68 1 0.11],...
    'String',txt,...
    'FontSize',28,...
    'FontName','latex',...
    'EdgeColor','none',...
    'LineWidth',1.5,...
    'BackgroundColor','none',...
    'FitBoxToText','on');

h2 = subtightplot(4,3,2,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(qt_sum_anom(:,:,:,t2))'); shading flat;
mymap = colormap('redblue');
colormap(h2,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
title(['t = ',num2str(t2),' day'],'interpreter','latex','fontsize',28)
xticklabels([])
yticklabels([])
%title('$q_t$','interpreter','latex')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(b)']};
h_annot=annotation('textbox',...
    [0.368  0.81 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');

h3 = subtightplot(4,3,3,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(qt_sum_anom(:,:,:,t3))'); shading flat;
mymap = colormap('redblue');
colormap(h3,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
title(['t = ',num2str(t3),' day'],'interpreter','latex','fontsize',28)
xticklabels([])
yticklabels([])
%title('$q_t$','interpreter','latex')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(c)']};
h_annot=annotation('textbox',...
    [0.658  0.81 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');

h4 = subtightplot(4,3,4,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(-qAdv_sum_anom_posPV(:,:,:,t1))'); shading flat;
colorbar('FontSize',20,'FontName','latex','linewidth',1.5,'position',[0.93    0.7350    0.0116    0.19000]);
mymap = colormap('redblue');
colormap(h4,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
%title(['t = ',num2str(t4),' day'],'interpreter','latex','fontsize',28)
xticklabels([])
ylabel('z (m)','interpreter','latex')
%title('$q_t$','interpreter','latex')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(d)']};
h_annot=annotation('textbox',...
    [0.078  0.58 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
txt={['q_{adv}']};
h_annot=annotation('textbox',...
    [0.28  0.46 1 0.11],...
    'String',txt,...
    'FontSize',28,...
    'FontName','latex',...
    'FontWeight','normal',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');

h5 = subtightplot(4,3,5,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(-qAdv_sum_anom_posPV(:,:,:,t2))'); shading flat;
mymap = colormap('redblue');
colormap(h5,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
xticklabels([])
yticklabels([])
%title('$\int q_t$','interpreter','latex')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
ax.FontName = 'latex';
txt={['(e)']};
h_annot=annotation('textbox',...
    [0.368  0.58 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');

h6 = subtightplot(4,3,6,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(-qAdv_sum_anom_posPV(:,:,:,t3))'); shading flat;
mymap = colormap('redblue');
colormap(h6,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
xticklabels([])
yticklabels([])
%title('$<q_t>$','interpreter','latex')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(f)']};
h_annot=annotation('textbox',...
    [0.658  0.58 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');

h7 = subtightplot(4,3,7,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(-qF_sum_anom_posPV(:,:,:,t1))'); shading flat;
mymap = colormap('redblue');
colormap(h7,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
xticklabels([])
ylabel('z (m)','interpreter','latex')
%title('$<q_t>$','interpreter','latex')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(g)']};
h_annot=annotation('textbox',...
    [0.078  0.36 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
txt={['q_F']};
h_annot=annotation('textbox',...
    [0.28  0.24 1 0.11],...
    'String',txt,...
    'FontSize',28,...
    'FontName','latex',...
    'FontWeight','normal',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');


h8 = subtightplot(4,3,8,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(-qF_sum_anom_posPV(:,:,:,t2))'); shading flat;
colorbar('FontSize',20,'FontName','latex','linewidth',1.5,'position',[0.93    0.074    0.0116    0.6204]);
mymap = colormap('redblue');
colormap(h8,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
xticklabels([])
yticklabels([])
%title('$<q_t>$','interpreter','latex')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(h)']};
h_annot=annotation('textbox',...
    [0.368  0.36 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');

h9 = subtightplot(4,3,9,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(-qF_sum_anom_posPV(:,:,:,t3))'); shading flat;
mymap = colormap('redblue');
colormap(h9,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
xticklabels([])
yticklabels([])
%title('$\int qF$','interpreter','latex')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
ax.FontName = 'latex';
txt={['(i)']};
h_annot=annotation('textbox',...
    [0.658  0.36 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');


h10 = subtightplot(4,3,10,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(-qB_sum_anom_posPV(:,:,:,t1))'); shading flat;
mymap = colormap('redblue');
colormap(h10,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
ylabel('z (m)','interpreter','latex')
xlabel('y (km)')
%title('$\int qF $','interpreter','latex')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(j)']};
h_annot=annotation('textbox',...
    [0.078  0.14 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');
txt={['q_B']};
h_annot=annotation('textbox',...
    [0.28  0.02 1 0.11],...
    'String',txt,...
    'FontSize',28,...
    'FontName','latex',...
    'EdgeColor','none',...
    'LineWidth',1.5,...
    'BackgroundColor','none',...
    'FitBoxToText','on');


h11 = subtightplot(4,3,11,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(-qB_sum_anom_posPV(:,:,:,t2))'); shading flat;
mymap = colormap('redblue');
colormap(h11,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
yticklabels([])
xlabel('y (km)')
%title('$\int qF $','interpreter','latex')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(k)']};
h_annot=annotation('textbox',...
    [0.368  0.14 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');


h12 = subtightplot(4,3,12,gap,marg_h,marg_w);
pcolor(YY/m1km,ZZ,squeeze(-qB_sum_anom_posPV(:,:,:,t3))'); shading flat;
mymap = colormap('redblue');
colormap(h12,mymap)
caxis(clim)
hold on;
set(gca,'fontsize',fontsize)
yticklabels([])
xlabel('y (km)')
%title('$\int qF $','interpreter','latex')
xlim(x_lim);
ylim(y_lim);
fill(x2, inBetween, c_ice,'EdgeColor','none');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
txt={['(l)']};
h_annot=annotation('textbox',...
    [0.658  0.14 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','latex',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate plume means %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% calculate d, distance from slope
Nd = 50;
d = zeros(Ny,Nd);
% Ns = 400;
Ns = 100:300;
for j = 1:Ny
    k = find(zz==floor(icetopo(j))+dz/2);
    if icetopo(j)-zz(k)>0
        d(j,:) = icetopo(j)-zz(k:k+(Nd-1));
    elseif icetopo(j)-zz(k)<0
        d(j,:) = icetopo(j)-zz(k+1:k+Nd);
    end
end
dd = mean(d(Ns,:),1); dd = dd-dd(1);
[D,NT] = meshgrid(dd,[1:Nt-1]);


qdot_mean = plume_mean(qt,Nd,Ns,icetopo,zz,dz);
qF_mean = plume_mean(qF,Nd,Ns,icetopo,zz,dz);
qB_mean = plume_mean(qB,Nd,Ns,icetopo,zz,dz);
qAdv_mean = plume_mean(qAdv,Nd,Ns,icetopo,zz,dz);
bdot_mean = plume_mean(bdot,Nd,Ns,icetopo,zz,dz);
q_mean = plume_mean(q,Nd,Ns,icetopo,zz,dz);
qt_sum_anom_mean = plume_mean(qt_sum_anom,Nd,Ns,icetopo,zz,dz);
qF_sum_anom_mean = plume_mean(qF_sum_anom,Nd,Ns,icetopo,zz,dz);
qB_sum_anom_mean = plume_mean(qB_sum_anom,Nd,Ns,icetopo,zz,dz);
qAdv_sum_anom_mean = plume_mean(qAdv_sum_anom,Nd,Ns,icetopo,zz,dz);
q_sum_anom_mean = plume_mean(q_sum_anom,Nd,Ns,icetopo,zz,dz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Enstrophy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Total, mean and perturbation enstrophies
E = plume_mean(0.5*q.^2,Nd,Ns,icetopo,zz,dz);
E_mean = 0.5*q_mean.^2;
E_pert = E - E_mean;

%%% Mean and perturbation enstrophy tendencies
Et_mean = diff(E_mean,1,2)./dt;
Et_pert = diff(E_pert,1,2)./dt;

%%% Total, mean and perturbation terms in the enstrophy budget

Eqdot = plume_mean(0.5*(q(:,:,:,1:end-1)+q(:,:,:,2:end)).*qdot,Nd,Ns,icetopo,zz,dz);
Eqdot_mean = 0.5*(q_mean(:,1:end-1)+q_mean(:,2:end)).*qdot_mean;
Eqdot_pert = Eqdot - Eqdot_mean;

EqB = plume_mean(0.5*(q(:,:,:,1:end-1)+q(:,:,:,2:end)).*qB,Nd,Ns,icetopo,zz,dz);
EqB_mean = 0.5*(q_mean(:,1:end-1)+q_mean(:,2:end)).*qB_mean;
EqB_pert = EqB - EqB_mean;

EqF = plume_mean(0.5*(q(:,:,:,1:end-1)+q(:,:,:,2:end)).*qF,Nd,Ns,icetopo,zz,dz);
EqF_mean = 0.5*(q_mean(:,1:end-1)+q_mean(:,2:end)).*qF_mean;
EqF_pert = EqF - EqF_mean;

EqAdv = plume_mean(0.5*(q(:,:,:,1:end-1)+q(:,:,:,2:end)).*qAdv,Nd,Ns,icetopo,zz,dz);
EqAdv_mean = 0.5*(q_mean(:,1:end-1)+q_mean(:,2:end)).*qAdv_mean;
EqAdv_pert = EqAdv - EqAdv_mean;

% Eterm1 = plume_mean(0.5*(q(:,:,:,1:end-1)+q(:,:,:,2:end)).*term1,Nd,Ns,icetopo,zz,dz);
% Eterm1_mean = 0.5*(q_mean(:,1:end-1)+q_mean(:,2:end)).*term1_mean;
% Eterm1_pert = Eterm1 - Eterm1_mean;
% 
% Eterm2 = plume_mean(0.5*(q(:,:,:,1:end-1)+q(:,:,:,2:end)).*term2,Nd,Ns,icetopo,zz,dz);
% Eterm2_mean = 0.5*(q_mean(:,1:end-1)+q_mean(:,2:end)).*term2_mean;
% Eterm2_pert = Eterm2 - Eterm2_mean;
% 
% Eterm3 = plume_mean(0.5*(q(:,:,:,1:end-1)+q(:,:,:,2:end)).*term3,Nd,Ns,icetopo,zz,dz);
% Eterm3_mean = 0.5*(q_mean(:,1:end-1)+q_mean(:,2:end)).*term3_mean;
% Eterm3_pert = Eterm3 - Eterm3_mean;
% 
% Eterm4 = plume_mean(0.5*(q(:,:,:,1:end-1)+q(:,:,:,2:end)).*term4,Nd,Ns,icetopo,zz,dz);
% Eterm4_mean = 0.5*(q_mean(:,1:end-1)+q_mean(:,2:end)).*term4_mean;
% Eterm4_pert = Eterm4 - Eterm4_mean;
% 
% Eterm5 = plume_mean(0.5*(q(:,:,:,1:end-1)+q(:,:,:,2:end)).*term5,Nd,Ns,icetopo,zz,dz);
% Eterm5_mean = 0.5*(q_mean(:,1:end-1)+q_mean(:,2:end)).*term5_mean;
% Eterm5_pert = Eterm5 - Eterm5_mean;

%%% Variance ratio budget terms
ENqdot = Eqdot_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eqdot_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2;
ENqB = EqB_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*EqB_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2;
ENqF = EqF_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*EqF_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2;
ENqAdv = EqAdv_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*EqAdv_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2;

% ENterm1 = Eterm1_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eterm1_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2;
% ENterm2 = Eterm2_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eterm2_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2;
% ENterm3 = Eterm3_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eterm3_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2;
% ENterm4 = Eterm4_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eterm4_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2;
% ENterm5 = Eterm5_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eterm5_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2;

%%

scrsz = get(0,'ScreenSize');
fontsize = 22;

framepos = [0 scrsz(4)/2 scrsz(3) scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(11);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');
gap = [0.05 0.02];
marg_h = [0.065 0.065];
marg_w = [0.04 0.03];

subtightplot(3,6,1,gap,marg_h,marg_w)
pcolor(YY',ZZ'-(YY'-4e4)*s,squeeze(q(1,:,:,120)));
shading flat;colorbar;caxis([-1 1]*3e-8);colormap redblue; set(gca,'YLim',[-20 0]);
title('snapshot of q','interpreter','latex','fontsize',20)

subtightplot(3,6,2,gap,marg_h,marg_w)
pcolor(NT,D,sqrt(E_pert(:,2:end))'-abs(q_mean(:,2:end))');
colorbar; shading flat; caxis([-1 1]*1e-8); colormap redblue;
xticklabels([]); hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$E^{\prime} - |\bar{q}|$','interpreter','latex','fontsize',20)

subtightplot(3,6,3,gap,marg_h,marg_w)
pcolor(NT,D,E_mean(:,2:end)');
colorbar;shading flat;caxis([-1 1]*1e-16);colormap redblue;
xticklabels([]); yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$\bar{E}$','interpreter','latex','fontsize',20)

subtightplot(3,6,4,gap,marg_h,marg_w)
pcolor(NT,D,E_pert(:,2:end)');
colorbar;shading flat;caxis([-1 1]*1e-16);colormap redblue;
xticklabels([]); yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$E^{\prime}$','interpreter','latex','fontsize',20)

subtightplot(3,6,5,gap,marg_h,marg_w)
pcolor(NT,D,(E_pert(:,2:end)'./E_mean(:,2:end)'));
colorbar;shading flat;caxis([-10 10]);colormap redblue;
yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$E^{\prime}/\bar{E}$','interpreter','latex','fontsize',20)

subtightplot(3,6,6,gap,marg_h,marg_w)
pcolor(NT,D,Et_mean');
colorbar;shading flat;caxis([-1 1]*1e-21);colormap redblue;
xticklabels([]); yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$\bar{E_t}$','interpreter','latex','fontsize',20)

subtightplot(3,6,7,gap,marg_h,marg_w);
pcolor(NT,D,Et_pert');colorbar;shading flat;caxis([-1 1]*1e-21);colormap redblue;
xticklabels([]); yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$E_t^{\prime}$','interpreter','latex','fontsize',20)

subtightplot(3,6,8,gap,marg_h,marg_w)
pcolor(NT,D,Et_pert'-Et_mean');
colorbar;shading flat;caxis([-1 1]*1e-21);colormap redblue;
xticklabels([]); hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$E_t^{\prime} - \bar{E_t}$','interpreter','latex','fontsize',20)

subtightplot(3,6,9,gap,marg_h,marg_w)
pcolor(NT,D,cumsum(EqB_mean',1)*dt);
colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
xticklabels([]); yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$\bar{E_{qB}}$','interpreter','latex','fontsize',20)

subtightplot(3,6,10,gap,marg_h,marg_w)
pcolor(NT,D,cumsum(EqB_pert',1)*dt);
colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
xticklabels([]); yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$E_{qB}^{\prime}$','interpreter','latex','fontsize',20)

subtightplot(3,6,11,gap,marg_h,marg_w)
pcolor(NT,D,cumsum(EqF_mean',1)*dt);
colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
xticklabels([]); yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$\bar{E_{qF}}$','interpreter','latex','fontsize',20)

subtightplot(3,6,12,gap,marg_h,marg_w)
pcolor(NT,D,cumsum(EqF_pert',1)*dt);
colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
xticklabels([]); yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$E_{qF}^{\prime}$','interpreter','latex','fontsize',20)

subtightplot(3,6,13,gap,marg_h,marg_w)
pcolor(NT,D,cumsum(EqAdv_mean',1)*dt);
colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
xticklabels([]); yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$\bar{E_{qAdv}}$','interpreter','latex','fontsize',20)

subtightplot(3,6,14,gap,marg_h,marg_w)
pcolor(NT,D,cumsum(EqAdv_pert',1)*dt);
colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$E_{qAdv}^{\prime}$','interpreter','latex','fontsize',20)

subtightplot(3,6,15,gap,marg_h,marg_w)
pcolor(NT,D,cumsum(ENqdot',1)*dt);
colorbar;shading flat;caxis([-10 10]);colormap redblue; title('qt');
yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$EN_{q_t}$','interpreter','latex','fontsize',20)

subtightplot(3,6,16,gap,marg_h,marg_w)
pcolor(NT,D,cumsum(ENqB',1)*dt);
colorbar;shading flat;caxis([-10 10]);colormap redblue; title('qB');
yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$EN_{qB}$','interpreter','latex','fontsize',20)

subtightplot(3,6,17,gap,marg_h,marg_w);
pcolor(NT,D,cumsum(ENqF',1)*dt);
colorbar;shading flat;caxis([-10 10]);colormap redblue; title('qF');
yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$EN_{qF}$','interpreter','latex','fontsize',20)

subtightplot(3,6,18,gap,marg_h,marg_w)
pcolor(NT,D,cumsum(ENqAdv',1)*dt);
colorbar;shading flat;caxis([-10 10]);colormap redblue; title('qAdv');
yticklabels([]);hold on;
plot([1:Nt],6*ones(length([1:Nt])),'linewidth',1,'color','k'); plot([1:Nt],10*ones(length([1:Nt])),'linewidth',1,'color','k');
title('$EN_{qAdv}$','interpreter','latex','fontsize',20)

% figure(8);pcolor(NT,D,cumsum(Eqdot_mean',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
% figure(9);pcolor(NT,D,cumsum(Eqdot_pert',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
% 
% figure(10);pcolor(NT,D,cumsum(Eterm1_mean',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
% figure(11);pcolor(NT,D,cumsum(Eterm1_pert',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
% 
% figure(12);pcolor(NT,D,cumsum(Eterm2_mean',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
% figure(13);pcolor(NT,D,cumsum(Eterm2_pert',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
% 
% figure(14);pcolor(NT,D,cumsum(Eterm3_mean',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
% figure(15);pcolor(NT,D,cumsum(Eterm3_pert',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
% 
% figure(16);pcolor(NT,D,cumsum(Eterm4_mean',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
% figure(17);pcolor(NT,D,cumsum(Eterm4_pert',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
% 
% figure(18);pcolor(NT,D,cumsum(Eterm5_mean',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
% figure(19);pcolor(NT,D,cumsum(Eterm5_pert',1)*dt);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;

% figure(28);pcolor(NT,D,cumsum(ENterm1',1)*dt);colorbar;shading flat;caxis([-10 10]);colormap redblue; title('term1');
% figure(29);pcolor(NT,D,cumsum(ENterm2',1)*dt);colorbar;shading flat;caxis([-10 10]);colormap redblue; title('term2');
% figure(30);pcolor(NT,D,cumsum(ENterm3',1)*dt);colorbar;shading flat;caxis([-10 10]);colormap redblue; title('term3');
% figure(31);pcolor(NT,D,cumsum(ENterm4',1)*dt);colorbar;shading flat;caxis([-10 10]);colormap redblue; title('term4');
% figure(32);pcolor(NT,D,cumsum(ENterm5',1)*dt);colorbar;shading flat;caxis([-10 10]);colormap redblue; title('term5');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scrsz = get(0,'ScreenSize');
fontsize = 26;
framepos = [0 scrsz(4)/2 scrsz(3) scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(30);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');

h1 = subplot(2,3,1);
pcolor(NT,D,Qt_mean'); colorbar; shading interp;
mymap = colormap('redblue');
colormap(h1,mymap)
caxis([-5 5]*1e-13)
hold on;
ylim([1 25])
set(gca,'fontsize',fontsize)
ylabel('d (m)')
xticklabels([])
title('$\Delta_t q$','interpreter','latex')

h2 = subplot(2,3,2);
pcolor(NT,D,qdot_mean'); colorbar; shading interp;
mymap = colormap('redblue');
colormap(h1,mymap)
caxis([-5 5]*1e-13)
hold on;
ylim([1 25])
set(gca,'fontsize',fontsize)
ylabel('d (m)')
xticklabels([])
title('$\dot{q}$','interpreter','latex')

h3 = subplot(2,3,3);
pcolor(NT,D,(qdot_mean - Qt_mean)'); colorbar; shading interp;
mymap = colormap('redblue');
colormap(h3,mymap)
caxis([-5 5]*1e-13)
hold on;
ylim([1 25])
set(gca,'fontsize',fontsize)
ylabel('d (m)')
xticklabels([])
title('$\dot{q} - \Delta_t q$','interpreter','latex')

h4 = subplot(2,3,4);
pcolor(NT,D,(term1_mean - qdot_mean)'); colorbar; shading interp;
mymap = colormap('redblue');
colormap(h4,mymap)
caxis([-5 5]*1e-13)
hold on;
ylim([1 25])
set(gca,'fontsize',fontsize)
ylabel('d (m)')
xticks([24,120,240,360]);
xticklabels({'1','5','10','15'})
xlabel('days')
title('$fD - \dot{q}$','interpreter','latex')

h5 = subplot(2,3,5);
pcolor(NT,D,(qvert_dot_mean - qdot_mean)'); colorbar; shading interp;
mymap = colormap('redblue');
colormap(h5,mymap)
caxis([-5 5]*1e-13)
hold on;
ylim([1 25])
set(gca,'fontsize',fontsize)
ylabel('d (m)')
xticks([24,120,240,360]);
xticklabels({'1','5','10','15'})
xlabel('days')
title('$\dot{q_{vert}} - \dot{q}$','interpreter','latex')

h6 = subplot(2,3,6);
pcolor(NT,D,(qbc_dot_mean - qdot_mean)'); colorbar; shading interp;
mymap = colormap('redblue');
colormap(h6,mymap)
caxis([-5 5]*1e-13)
hold on;
ylim([1 25])
set(gca,'fontsize',fontsize)
ylabel('d (m)')
xticks([24,120,240,360]);
xticklabels({'1','5','10','15'})
xlabel('days')
title('$\dot{q_{bc}} - \dot{q}$','interpreter','latex')


%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Kinetic Energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Uc = U; Vc = NaN(Nx,Ny,Nr,Nt); Wc = NaN(Nx,Ny,Nr,Nt);
% Vc(:,1:Ny-1,:,:) = (V(:,1:Ny-1,:,:)+V(:,2:Ny,:,:))/2; Vc(:,Ny,:,:) = 2*Vc(:,Ny-1,:,:) - Vc(:,Ny-2,:,:);
% Wc(:,:,2:Nr,:) = (W(:,:,2:Nr,:)+W(:,:,1:Nr-1,:))/2; Wc(:,:,1,:) = 2*Wc(:,:,2,:) - Wc(:,:,3,:);
% Uc = iceNaN(Uc,Nx,Ny,Nr,icetopo,zz); Vc = iceNaN(Vc,Nx,Ny,Nr,icetopo,zz); Wc = iceNaN(Wc,Nx,Ny,Nr,icetopo,zz);
% KE = Uc.^2 + Vc.^2 + Wc.^2 ;
%
% % advective flux terms
% vKE = Vc.*KE;
% wKE = Wc.*KE;
% vKEdy = diff(vKE,1,2)./dy;
% vKEdy = (vKEdy(:,:,2:Nr,:) + vKEdy(:,:,1:Nr-1,:))/2;
% vKEdy = interp_yz(vKEdy);
% wKEdz = (wKE(:,:,[1:Nr-1],1:Nt) - wKE(:,:,[2:Nr],1:Nt))./dz;
% wKEdz = (wKEdz(:,2:Ny,:,:) + wKEdz(:,1:Ny-1,:,:))/2;
% wKEdz = interp_yz(wKEdz);
%
% KEadvFlux = -(vKEdy + wKEdz);
%
% % pressure flux terms
% P = NaN(Nx,Ny,Nr,Nt); P(:,:,:,1:Nt-1) = PH + PNH(:,:,:,1:Nt-1); P(:,:,:,Nt) = 2*P(:,:,:,Nt-1) - P(:,:,:,Nt-2);
% vP = Vc.*P;
% wP = Wc.*P;
% vPdy = 1/rho0*diff(vP,1,2)./dy;
% vPdy = (vPdy(:,:,2:Nr,:) + vPdy(:,:,1:Nr-1,:))/2;
% vPdy = interp_yz(vPdy);
% wPdz = 1/rho0*(wP(:,:,[1:Nr-1],1:Nt) - wP(:,:,[2:Nr],1:Nt))./dz;
% wPdz = (wPdz(:,2:Ny,:,:) + wPdz(:,1:Ny-1,:,:))/2;
% wPdz = interp_yz(wPdz);
%
% KEpFlux = - (vPdy + wPdz);
%
% % viscous flux terms
% KEdy2(:,2:Ny-1,:,:) = (KE(:,3:Ny,:,:) - 2*KE(:,2:Ny-1,:,:) + KE(:,1:Ny-2,:,:))./(dy^2); % second derivative
% KEdy2(:,1,:,:) = 2*KEdy2(:,2,:,:) - KEdy2(:,3,:,:); KEdy2(:,Ny,:,:) = 2*KEdy2(:,Ny-1,:,:) - KEdy2(:,Ny-2,:,:); % adding boundaries
% KEdy2 = (KEdy2(:,2:Ny,:,:) + KEdy2(:,1:Ny-1,:,:))/2; KEdy2 = (KEdy2(:,:,2:Nr,:) + KEdy2(:,:,1:Nr-1,:))/2;
% KEdy2 = interp_yz(KEdy2);
% KEviscFlux_h = viscA4/(dy^2).*KEdy2;
% KEdz2(:,:,2:Nr-1,:) = (KE(:,:,1:Nr-2,:) - 2*KE(:,:,2:Nr-1,:,:) + KE(:,:,3:Nr,:))./(dz^2);
% KEdz2(:,:,1,:) = 2*KEdz2(:,:,2,:) - KEdz2(:,:,3,:); KEdz2(:,:,Nr,:) = 2*KEdz2(:,:,Nr-1,:) - KEdz2(:,:,Nr-2,:);
% KEdz2 = (KEdz2(:,2:Ny,:,:) + KEdz2(:,1:Ny-1,:,:))/2; KEdz2 = (KEdz2(:,:,2:Nr,:) + KEdz2(:,:,1:Nr-1,:))/2;
% KEdz2 = interp_yz(KEdz2);
% KEviscFlux_z = viscAr*KEdz2;
%
% KEviscFlux = - (KEviscFlux_z + KEviscFlux_h);
%
% % buoyancy source
% wb = Wc.*b_;
%
% % viscous dissipation terms
% dUdy = NaN(Nx,Ny,Nr,Nt); dUdy(:,1:Ny-1,:,:) = diff(U(:,:,:,1:Nt),1,2)./dy; % derivative wrt y, shifts to right face center
% dUdy(:,:,2:Nr,:) = (dUdy(:,:,[2:Nr],:)+dUdy(:,:,[1:Nr-1],:))/2; % average to top right corner
% dUdy(:,:,1,:) = 2*dUdy(:,:,2,:)-dUdy(:,:,3,:); dUdy(:,Ny,:,:) = 2*dUdy(:,Ny-1,:,:)-dUdy(:,Ny-2,:,:);
% dVdy = NaN(Nx,Ny,Nr,Nt); dVdy(:,1:Ny-1,:,:) = diff(Vc(:,:,:,1:Nt),1,2)./dy; % derivative wrt y, shifts to right face center
% dVdy(:,:,2:Nr,:) = (dVdy(:,:,[2:Nr],:)+dVdy(:,:,[1:Nr-1],:))/2; % average to top right corner
% dVdy(:,:,1,:) = 2*dVdy(:,:,2,:)-dVdy(:,:,3,:); dVdy(:,Ny,:,:) = 2*dVdy(:,Ny-1,:,:)-dVdy(:,Ny-2,:,:);
% dWdy = NaN(Nx,Ny,Nr,Nt); dWdy(:,1:Ny-1,:,:) = diff(Wc(:,:,:,1:Nt),1,2)./dy; % derivative wrt y, shifts to right face center
% dWdy(:,:,2:Nr,:) = (dWdy(:,:,[2:Nr],:)+dWdy(:,:,[1:Nr-1],:))/2; % average to top right corner
% dWdy(:,:,1,:) = 2*dWdy(:,:,2,:)-dWdy(:,:,3,:); dWdy(:,Ny,:,:) = 2*dWdy(:,Ny-1,:,:)-dWdy(:,Ny-2,:,:);
% KEdiss_h = viscA4/(dy^2)*(dUdy.^2 + dVdy.^2 + dWdy.^2);
%
% dUdz = NaN(Nx,Ny,Nr,Nt); dUdz(:,:,2:Nr,:) = (U(:,:,[1:Nr-1],1:Nt) - U(:,:,[2:Nr],1:Nt))./dz; % derivative wrt z, shifts to top face center
% dUdz(:,1:Ny-1,:,:) = (dUdz(:,[2:Ny],:,:)+dUdz(:,[1:Ny-1],:,:))/2; % average to top right corner
% dUdz(:,:,1,:) = 2*dUdz(:,:,2,:)-dUdz(:,:,3,:); dUdz(:,Ny,:,:) = 2*dUdz(:,Ny-1,:,:)-dUdz(:,Ny-2,:,:); % add boundaries
% dVdz = NaN(Nx,Ny,Nr,Nt); dVdz(:,:,2:Nr,:) = (Vc(:,:,[1:Nr-1],1:Nt) - Vc(:,:,[2:Nr],1:Nt))./dz; % derivative wrt z, shifts to top face center
% dVdz(:,1:Ny-1,:,:) = (dVdz(:,[2:Ny],:,:)+dVdz(:,[1:Ny-1],:,:))/2; % average to top right corner
% dVdz(:,:,1,:) = 2*dVdz(:,:,2,:)-dVdz(:,:,3,:); dVdz(:,Ny,:,:) = 2*dVdz(:,Ny-1,:,:)-dVdz(:,Ny-2,:,:); % add boundaries
% dWdz = NaN(Nx,Ny,Nr,Nt); dWdz(:,:,2:Nr,:) = (Wc(:,:,[1:Nr-1],1:Nt) - Wc(:,:,[2:Nr],1:Nt))./dz; % derivative wrt z, shifts to top face center
% dWdz(:,1:Ny-1,:,:) = (dWdz(:,[2:Ny],:,:)+dWdz(:,[1:Ny-1],:,:))/2; % average to top right corner
% dWdz(:,:,1,:) = 2*dWdz(:,:,2,:)-dWdz(:,:,3,:); dWdz(:,Ny,:,:) = 2*dWdz(:,Ny-1,:,:)-dWdz(:,Ny-2,:,:); % add boundaries
% KEdiss_z = viscAr*(dUdz.^2 + dVdz.^2 + dWdz.^2);
%
% KEdiss = KEdiss_h + KEdiss_z;
%
% % KE tendency
% KEt(:,:,:,1:Nt-1) = (KE(:,:,:,2:Nt) - KE(:,:,:,1:Nt-1))/dt;
% KEt(:,:,:,Nt) = 2*KEt(:,:,:,Nt-1) - KEt(:,:,:,Nt-2);
%
% % KEt_tot
% KEt_tot = KEadvFlux + KEpFlux + KEviscFlux + KEdiss;
%
% % residual
% res = KEt - KEt_tot;
%
%
% figure;
% hold on;
% plot(squeeze(nansum(nansum(KEadvFlux(:,:,:,360),1),3)));
% plot(squeeze(nansum(nansum(KEpFlux(:,:,:,360),1),3)));
% plot(squeeze(nansum(nansum(KEviscFlux(:,:,:,360),1),3)));
% plot(squeeze(nansum(nansum(wb(:,:,:,360),1),3)));
% plot(squeeze(nansum(nansum(KEdiss(:,:,:,360),1),3)));
% plot(squeeze(nansum(nansum(KEt(:,:,:,360),1),3)));
% plot(squeeze(nansum(nansum(KEt_tot(:,:,:,360),1),3)));
% legend('adv','p','visc','wb','diss','KEt','KEt_tot')


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

function Ahov = plume_mean(AA,Nd,Ns,icetopo,zz,dz)
Ny = size(AA,2);
Nt = size(AA,4);
AA_= zeros(Ny-1,Nd,Nt);
for j = 1:Ny-1
    k = find(zz==floor(icetopo(j))+(dz)/2);
    AA_(j,:,:) = AA(1,j,k:k+(Nd-1),:);
end
Ahov = squeeze(nanmean(AA_(Ns,:,:),1));
end

function Ahov = plume_var(AA,Nd,Ns,icetopo,zz,dz)
Ny = size(AA,2);
Nt = size(AA,4);
AA_= zeros(Ny-1,Nd,Nt);
for j = 1:Ny-1
    k = find(zz==floor(icetopo(j))+(dz)/2);
    AA_(j,:,:) = AA(1,j,k:k+(Nd-1),:);
end
% Ahov = squeeze(max(AA_(1:Ns,:,:),[],1,'omitnan'));
Ahov = squeeze(var(AA_(Ns,:,:),1,'omitnan'));
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
