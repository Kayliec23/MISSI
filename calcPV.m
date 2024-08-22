%function [Q,Qsym,Qcent,Qconv,Ri,Rii] = calcPV(T,S,U,V,W)

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
g = 9.8;
rho0 = 1027;
L = 334e3;
rhoi = 917;
t1min = 50;
t1hour = 60*t1min;
t1day = 24*t1hour;
dt = t1day;
[YY,ZZ] = meshgrid(yy,zz);

addpath /Users/Kayliec23/Desktop/analysis
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculating PV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear sigmatheta theta dWdy dVdz dUdz dUdy dRHOdy dRHOdz dbdy dbdz Q

% Calculate potential density (at grid center)
sigmatheta = zeros(Nx,Ny,Nr,Nt);
for n = 1:Nt
    SS = squeeze(S(1,:,:,n));
    TT = squeeze(T(1,:,:,n));
    theta(:,:) = densjmd95(SS,TT,0);
    sigmatheta(1,:,:,n) = theta -1000;
end


% calculate ZETAx
dWdy = NaN(Nx,Ny,Nr,Nt); dWdy(:,1:Ny-1,:,:) = diff(W(:,:,:,1:Nt),1,2)./dy; % derivative wrt y, moves to bottom right grid
dWdy(:,Ny,:,:) = 2*dWdy(:,Ny-1,:,:) - dWdy(:,Ny-2,:,:); % add boundary at the end
dVdz = NaN(Nx,Ny,Nr,Nt); dVdz(:,:,2:Nr,:) = diff(V(:,:,:,1:Nt),1,3)./dz; % derivative wrt z, moves to bottom left grid, shifted to bottom right
dVdz(:,:,1,:) = 2*dVdz(:,:,2,:) - dVdz(:,:,3,:); % add boundary at the beginning
ZETAx = dWdy - dVdz;

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
Qsym = dbdy.*ZETAy; % symmetric
Qconv = dbdz.*f0; % convective
Qcent = dbdz.*ZETAz; % centrifugal
Q = dbdy.*ZETAy + (f0 + ZETAz).*dbdz; % total

q = iceNaN(Q,Nx,Ny-1,Nr-1,icetopo,zz);
fq = f0*q;

q1 = Qconv+ Qcent; q1 = iceNaN(q1,Nx,Ny-1,Nr-1,icetopo,zz);
q2 = Qsym; q2 = iceNaN(q2,Nx,Ny-1,Nr-1,icetopo,zz);


% gradient Richardson number
Ri = f0.^2*(-g/rho0*dRHOdz)./((-g/rho0*dRHOdy).^2);
Ri(Ri==inf) = NaN; Ri(Ri==-inf) = NaN;
Ri = iceNaN(Ri,Nx,Ny,Nr,icetopo,zz);

Rii = squeeze(Ri);
for j = 1:Ny
    for k =1:Nr
        for n = 1:Nt
            if fq(1,j,k,n)>0
                Rii(j,k,n) = NaN;
            end
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PV budget %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time derivative
qt(:,:,:,1:Nt-1) = diff(q,1,4)./dt; qt(:,:,:,Nt) = 2*qt(:,:,:,Nt-1)-qt(:,:,:,Nt-2);
qt(:,:,:,2:Nt) = (qt(:,:,:,2:Nt)+qt(:,:,:,1:Nt-1))./2; qt(:,:,:,1) = 2*qt(:,:,:,2) - qt(:,:,:,3);

% advective term
vc(:,1:Ny-1,:,1:Nt) = (V(:,2:Ny,:,1:Nt) + V(:,1:Ny-1,:,1:Nt))/2; vc(:,Ny,:,:) = 2*vc(:,Ny-1,:,:) - vc(:,Ny-2,:,:);
qy(:,1:Ny-1,:,1:Nt) = diff(q,1,2)./dy; qy(:,Ny,:,:) = 2*qy(:,Ny-1,:,:) - qy(:,Ny-2,:,:);
qy(:,:,[2:Nr],:,:) = (qy(:,:,[2:Nr],:)+qy(:,:,[1:Nr-1],:))/2; qy(:,:,1,:) = 2*qy(:,:,2,:) - qy(:,:,3,:);

vqy = vc.*qy; 

wc(:,:,2:Nr,1:Nt) = (W(:,:,2:Nr,1:Nt) + W(:,:,1:Nr-1,1:Nt))/2; wc(:,:,1,:) = 2*wc(:,:,2,:) - wc(:,:,3,:);
qz(:,:,2:Nr,1:Nt) = (q(:,:,1:Nr-1,:)-q(:,:,2:Nr,:))./dz; qz(:,:,1,:) = 2*qz(:,:,2,:) - qz(:,:,3,:);
qz(:,1:Ny-1,:,:) = (qz(:,1:Ny-1,:,:)+qz(:,2:Ny,:,:))/2; qz(:,Ny,:,:) = 2*qz(:,Ny-1,:,:) - qz(:,Ny-2,:,:);

wqz = wc.*qz; 

qAdv = - (vqy + wqz) ; 


% viscous
Uyy(:,2:Ny-1,:,:) = 0.5*(U(:,3:Ny,:,1:Nt) - 2*U(:,2:Ny-1,:,1:Nt) + U(:,1:Ny-2,:,1:Nt))./(dy.^2);
Uyy(:,1,:,:) = 2*Uyy(:,2,:,:) - Uyy(:,3,:,:); Uyy(:,Ny,:,:) = 2*Uyy(:,Ny-1,:,:) - Uyy(:,Ny-2,:,:);
Fh = viscA4./(dy.^2).*Uyy;

Uzz(:,:,2:Nr-1,:) = 0.5*(U(:,:,1:Nr-2,1:Nt) - 2*U(:,:,2:Nr-1,1:Nt) + U(:,:,3:Nr,1:Nt))./(dz.^2);
Uzz(:,:,1,:) = 2*Uzz(:,:,2,:) - Uzz(:,:,3,:); Uzz(:,:,Nr,:) = 2*Uzz(:,:,Nr-1,:) - Uzz(:,:,Nr-2,:); 
Fz = viscAr.*Uzz;

F = Fz + Fh; 

dbdy_(:,1:Ny-1,:,:) = (dbdy(:,1:Ny-1,:,:) + dbdy(:,2:Ny,:,:))./2; dbdy_(:,Ny,:,:) = 2*dbdy_(:,Ny-1,:,:) - dbdy_(:,Ny-2,:,:);
dbdz_(:,:,2:Nr,:) = (dbdz(:,:,2:Nr,:)+dbdz(:,:,1:Nr-1,:))./2; dbdz_(:,:,1,:) = 2*dbdz_(:,:,2,:) - dbdz(:,:,3,:); 

qFj = dbdz_.*F;
qFk = -dbdy_.*F;
qFjy(:,1:Ny-1,:,:) = diff(qFj,1,2)./dy; qFjy(:,Ny,:,:) = 2*qFjy(:,Ny-1,:,:) - qFjy(:,Ny-2,:,:); 
qFjy(:,1:Ny-1,:,:) = (qFjy(:,2:Ny,:,:) + qFjy(:,1:Ny-1,:,:))./2; qFjy(:,Ny,:,:) = 2*qFjy(:,Ny-1,:,:) - qFjy(:,Ny-2,:,:); 
qFkz(:,:,2:Nr,:) = (qFk(:,:,1:Nr-1,:) - qFk(:,:,2:Nr,:))./dz; qFkz(:,:,1,:) = 2*qFkz(:,:,2,:) - qFkz(:,:,3,:);
qFkz(:,:,2:Nr,:) = (qFkz(:,:,2:Nr,:) + qFkz(:,:,1:Nr-1,:))./2; qFkz(:,:,1,:) = 2*qFkz(:,:,2,:) - qFkz(:,:,3,:);

qF = qFkz + qFjy; 


% diabatic
% Bdot = fwFlux.*L./rhoi./dy./dz;
% Bdot_dy(:,1:Ny-1,:,:) = diff(Bdot,1,1); Bdot_dy(:,Ny,:,:) = 2*Bdot_dy(:,Ny-1,:,:) - Bdot_dy(:,Ny-2,:,:); 
% Bdot_dz(:,:,2:Nr,:) = (Bdot(:,:,1:Nr-1,:) - Bdot(:,:,2:Nr,:))./dz; Bdot_dz(:,:,1,:) = 2*Bdot_dz(:,:,2,:) - Bdot_dz(:,:,3,:);
% 





dbdzz(:,:,2:Nr,:) = (dbdz(:,:,1:Nr-1,:) - dbdz(:,:,2:Nr,:))./dz; dbdzz(:,:,1,:) = 2*dbdzz(:,:,2,:) - dbdzz(:,:,3,:);
%dbdzzz(:,:,2:Nr,:) = (dbdzz(:,:,1:Nr-1,:) - dbdzz(:,:,2:Nr,:))./dz; dbdzzz(:,:,1,:) = 2*dbdzzz(:,:,2,:) - dbdzzz(:,:,3,:);

qBk = (f0-dUdy).*KPPdiffT(:,:,:,1:Nt).*dbdzz;

qB(:,:,2:Nr,:) = (qBk(:,:,1:Nr-1,:) - qBk(:,:,2:Nr,:))./dz; qB(:,:,1,:) = 2*qB(:,:,2,:) - qB(:,:,3,:);

% omega_a = f0-dUdy; 
% Bt_ = squeeze(1/rhoi/L*fwFlux);
% Bt = NaN(Ny,Nr,Nt-1);
% 
% for j =1:Ny
%     for k = 2:Nr
%         if zz(k-1)>icetopo(j) && zz(k)<icetopo(j)
%             Bt(j,k,:) = Bt_(j,:);
%         end
%     end
% end

%qB = squeeze(omega_a(:,:,:,1:Nt-1)).*Bt; 

% % qB_fComp_ = f0.*KPPdiffT(:,:,:,1:Nt).*dbdzz;
% % qB_fComp(:,:,2:Nr,:) = (qB_fComp_(:,:,1:Nr-1,:) - qB_fComp_(:,:,2:Nr,:))./dz; qB_fComp(:,:,1,:) = 2*qB_fComp(:,:,2,:) - qB_fComp(:,:,3,:);
% % qB_geoComp_ = -dUdy.*KPPdiffT(:,:,:,1:Nt).*dbdzz;
% % qB_geoComp(:,:,2:Nr,:) = (qB_geoComp_(:,:,1:Nr-1,:) - qB_geoComp_(:,:,2:Nr,:))./dz; qB_geoComp(:,:,1,:) = 2*qB_geoComp(:,:,2,:) - qB_geoComp(:,:,3,:);

qAdv_neg = negVal(qAdv);
qF_neg = negVal(qF);
qB_neg = negVal(qB);
% % qB_fComp_neg = negVal(qB_fComp);
% % qB_geoComp_neg = negVal(qB_geoComp);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting fq terms and Ri %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




[YY,ZZ] = meshgrid(yy,zz);
scrsz = get(0,'ScreenSize');
fontsize = 28;
framepos = [0 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(17);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');

hplot = 0.27;
wplot = 0.14;
t1 = 8; t2 = 20; t3 = 180;
use_linear_lims = 'true';
if (use_linear_lims)
    xlim_min = 5;
    xlim_max = 15;  
else
    xlim_min = 15; 
    xlim_max = 25;
end
ylim_min = -320;
ylim_max = -230;
curve1 = icetopo;
curve2 = zeros(1,length(curve1));
x1 = yy/m1km; %linspace(xlim_min, xlim_max,length(curve1));
x2 = [x1, fliplr(x1)];
inBetween = [curve1, fliplr(curve2)];
c_ice = [0.9 0.9 0.9];

%%% plotting fq
subplot('position',[0.08   0.6874    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(fq(1,:,:,t1))');
shading flat;
%colorbar('box','on','linewidth',2,'location','northoutside','position',[0.08+0.01    0.955    3*wplot-0.04   0.015]);
caxis([-5 5]*1e-13);
hold on;
%cmocean('balance','pivot',0);
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xticklabels([])
ylabel('z (m)')
set(gca,'Box','on','linewidth',2)
title('fQ')
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
txt={['t = 7 days']};
    h_annot=annotation('textbox',...
        [0.14  0.62 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');
txt={['(a)']};
    h_annot=annotation('textbox',...
        [0.085  0.835 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');
yticks([-290 -260 -230 -200])
yticklabels({'-290','-260','-230','-200'})
fill(x2, inBetween, c_ice,'EdgeColor','none');



subplot('position',[0.08    0.4034    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(fq(1,:,:,t2))');
shading flat;
%colorbar('box','on','linewidth',1.5);
caxis([-5 5]*1e-13);
hold on;
%cmocean('balance','pivot',0);
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xticklabels([])
ylabel('z (m)')
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
txt={['t = 30 days']};
    h_annot=annotation('textbox',...
        [0.135  0.336 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');
txt={['(e)']};
    h_annot=annotation('textbox',...
        [0.085  0.551 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');
yticks([-290 -260 -230 -200])
yticklabels({'-290','-260','-230','-200'});
fill(x2, inBetween, c_ice,'EdgeColor','none');


subplot('position',[0.08    0.1195    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(fq(1,:,:,t3))');
shading flat;
colorbar('box','on','linewidth',2,'location','northoutside','position',[0.08+0.01    0.0345    3*wplot-0.01   0.015]);
caxis([-5 5]*1e-13);
hold on;
%cmocean('balance','pivot',0);
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xlabel('y (km)');
ylabel('z (m)')
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
txt={['t = 180 days']};
    h_annot=annotation('textbox',...
        [0.125    0.0508 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');
txt={['(i)']};
    h_annot=annotation('textbox',...
        [0.085  0.2658 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');
xticks([9 12 15 18])
xticklabels({'9','12','15','18'})
yticks([-290 -260 -230 -200])
yticklabels({'-290','-260','-230','-200'});
fill(x2, inBetween, c_ice,'EdgeColor','none');



% %% finding q1 and q2 when fq<0 only
% q1 = Qconv+Qcent;
% q2 = Qsym;
% 
% for j = 1:Ny
%     for k = 1:Nr
%         for n = 1:Nt
%             if fq(1,j,k,n)>0
%                 q1(1,j,k,n) = NaN;
%                 q2(1,j,k,n) = NaN;
%             end
%         end
%     end
% end

%%% plotting q1
subplot('position',[0.235    0.6874    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*q1(1,:,:,t1))');
shading interp;
%colorbar('box','on','linewidth',2,'location','northoutside','position',[0.3556+0.01    0.955    wplot-0.02   0.015]);
caxis([-5 5]*1e-13);
hold on;
%cmocean('balance','pivot',0);
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xticklabels([])
yticklabels([]);
set(gca,'Box','on','linewidth',2)
title('fQ_{vert}')
xlim([xlim_min xlim_max]);
ylim([ylim_min ylim_max]);
fill(x2, inBetween, c_ice,'EdgeColor','none');
txt={['(b)']};
    h_annot=annotation('textbox',...
        [0.24  0.835 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');
    


subplot('position',[0.235    0.4034    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*q1(1,:,:,t2))');
shading interp;
%colorbar('box','on','linewidth',1.5);
caxis([-5 5]*1e-13);
hold on;
%cmocean('balance','pivot',0);
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xticklabels([])
yticklabels([]);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max]);
ylim([ylim_min ylim_max]);
fill(x2, inBetween, c_ice,'EdgeColor','none');
txt={['(f)']};
    h_annot=annotation('textbox',...
        [0.24  0.551 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');


subplot('position',[0.235    0.1195    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*q1(1,:,:,t3))');
shading interp;
%colorbar('box','on','linewidth',1.5);
caxis([-5 5]*1e-13);
hold on;
%cmocean('balance','pivot',0);
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xlabel('y (km)');
%xticklabels([])
yticklabels([]);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
xticks([9 12 15 18])
xticklabels({'9','12','15','18'});
fill(x2, inBetween, c_ice,'EdgeColor','none');
txt={['(j)']};
    h_annot=annotation('textbox',...
        [0.24  0.2658 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');


%%% plotting q2
subplot('position',[0.39    0.6874    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*q2(1,:,:,t1))');
shading interp;
%colorbar('box','on','linewidth',2,'location','northoutside','position',[0.4841+0.01    0.955    wplot-0.02   0.015]);
caxis([-5 5]*1e-13);
hold on;
%cmocean('balance','pivot',0);
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xticklabels([])
set(gca,'Box','on','linewidth',2)
title('fQ_{bc}')
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
yticklabels([]);
fill(x2, inBetween, c_ice,'EdgeColor','none');
txt={['(c)']};
    h_annot=annotation('textbox',...
        [0.395  0.835 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');
    


subplot('position',[0.39    0.4034    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*q2(1,:,:,t2))');
shading interp;
%colorbar('box','on','linewidth',1.5);
caxis([-5 5]*1e-13);
hold on;
%cmocean('balance','pivot',0);
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xticklabels([]);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
yticklabels([]);
fill(x2, inBetween, c_ice,'EdgeColor','none');
txt={['(g)']};
    h_annot=annotation('textbox',...
        [0.395  0.551 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');


subplot('position',[0.39    0.1195    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*q2(1,:,:,t3))');
shading interp;
%colorbar('box','on','linewidth',1.5);
caxis([-5 5]*1e-13);
hold on;
%cmocean('balance','pivot',0);
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xlabel('y (km)');
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
yticklabels([]);
xticks([9 12 15 18])
xticklabels({'9','12','15','18'});
fill(x2, inBetween, c_ice,'EdgeColor','none');
txt={['(k)']};
    h_annot=annotation('textbox',...
        [0.395  0.2658 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');


%%% plotting Ri
h = subplot('position',[0.546    0.6874    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(Rii(:,:,t1))');
shading flat;
%cb = colorbar('box','on','linewidth',2,'location','northoutside','position',[0.6131+0.01    0.955    wplot-0.02   0.015]);
hold on;
caxis([-0.33 1.33]);
mymap = cmocean('oxy');
colormap(h,mymap);
set(gca,'fontsize',22)
xticklabels([])
set(gca,'Box','on','linewidth',2)
title('Ri')
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
yticklabels([]);
fill(x2, inBetween, c_ice,'EdgeColor','none');
txt={['(d)']};
    h_annot=annotation('textbox',...
        [0.551  0.835 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');


h = subplot('position',[0.546    0.4034    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(Rii(:,:,t2))');
shading flat;
%colorbar('box','on','linewidth',1.5);
hold on;
caxis([-0.33 1.33]);
mymap = cmocean('oxy');
colormap(h,mymap);
set(gca,'fontsize',22)
xticklabels([])
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
yticklabels([]);
fill(x2, inBetween, c_ice,'EdgeColor','none');
txt={['(h)']};
    h_annot=annotation('textbox',...
        [0.551  0.551 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');

h = subplot('position',[0.546    0.1195   wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(Rii(:,:,t3))');
shading flat;
cb = colorbar('box','on','linewidth',2,'location','northoutside','position',[0.546+0.01    0.0345    wplot-0.01   0.015]);
hold on;
caxis([-0.33 1.33]);
cb.YTick = [0 1];
mymap = cmocean('oxy');
colormap(h,mymap);
set(gca,'fontsize',22)
xlabel('y (km)');
yticklabels([]);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
xticks([9 12 15 18])
xticklabels({'9','12','15','18'});
fill(x2, inBetween, c_ice,'EdgeColor','none');
txt={['(l)']};
h_annot=annotation('textbox',...
    [0.551  0.2658 1 0.11],...
    'String',txt,...
    'FontSize',20,...
    'FontName','Arial',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none',...
    'FitBoxToText','on');




% %%
% %%% calculate d, distance from slope
% Nd = 50;
% d = zeros(Ny,Nd);
% Ns = 400;
% for j = 1:Ny
%     k = find(zz==floor(icetopo(j))+dz/2);
%     if icetopo(j)-zz(k)>0
%         d(j,:) = icetopo(j)-zz(k:k+(Nd-1));
%     elseif icetopo(j)-zz(k)<0
%         d(j,:) = icetopo(j)-zz(k+1:k+Nd);
%     end
% end
% dd = mean(d(1:Ns,:),1); dd = dd-dd(1);
% [D,NT] = meshgrid(dd,[1:Nt]);
% 
% Qt_mean = plume_mean(Qt,Nd,Ns,icetopo,zz,dz); 
% q_mean = plume_mean(q,Nd,Ns,icetopo,zz,dz);
% 
% figure; pcolor(NT,D,q_mean'); colorbar; shading flat; caxis([-5 5]*1e-8);cmocean('balance','pivot',0)
% 
%% make movie
[YY,ZZ] = meshgrid(yy,zz);
scrsz = get(0,'ScreenSize');
fontsize = 26;
framepos = [0 scrsz(4)/2 scrsz(3) scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(10);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');
M = moviein(Nt);
vid=VideoWriter('dS04_fq_Ri','MPEG-4');
vid.FrameRate = 10;
open(vid)

for n = 1:Nt-1
    
    subplot(1,2,1)
    pcolor(YY/m1km,ZZ,squeeze(fq(:,:,:,n))');
    shading flat;
    colorbar;
    caxis([-5 5]*1e-15);
    %cmocean('balance','pivot',0)
    colormap('redblue')
    %xlim([18 25])
    %ylim([-260 -190])
    
    subplot(1,2,2)
    pcolor(YY/m1km,ZZ,squeeze(qt(:,:,:,n))');
    shading interp;
    colorbar;
    caxis([-5 5]*1e-16);
    colormap('redblue')
    %cmocean('balance','pivot',0)
    %xlim([18 25])
    %ylim([-260 -190])
    
    M(n) = getframe(gcf);
    writeVideo(vid,M(n))
end

close(vid)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting PV budget %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rbmap = b2r(0,5e-13);

[YY,ZZ] = meshgrid(yy,zz);
scrsz = get(0,'ScreenSize');
fontsize = 28;
framepos = [0 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(13);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');
xlim_min = 0; xlim_max = 40;
ylim_min = -550; ylim_max = 0;
curve1 = icetopo;
curve2 = zeros(1,length(curve1));
x1 = yy/m1km; %linspace(xlim_min, xlim_max,length(curve1));
x2 = [x1, fliplr(x1)];
inBetween = [curve1, fliplr(curve2)];
c_ice = [0.9 0.9 0.9];

t1 = 5; t2 = 60;

sb1 = subplot(4,2,1); hold on;
pcolor(YY/m1km,ZZ,squeeze(qt(:,:,:,t1))'); shading flat; 
%colorbar;
caxis([-1 1]*1e-13);
xlim([0 40]);
title('q_t')
%cmocean('balance','pivot',0);
colormap(sb1,'redblue');
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
fill(x2, inBetween, c_ice,'EdgeColor','none');
set(gca,'fontsize',fontsize)
sb2 = subplot(4,2,2); hold on;
pcolor(YY/m1km,ZZ,squeeze(qt(:,:,:,t2))'); shading flat; 
colorbar;
caxis([-1 1]*1e-13);
xlim([0 40]);
title('q_t')
colormap(sb2,'redblue');  %cmocean('balance','pivot',0);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
yticklabels([]);
fill(x2, inBetween, c_ice,'EdgeColor','none');
set(gca,'fontsize',fontsize)

sb3 = subplot(4,2,3); hold on;
pcolor(YY/m1km,ZZ,squeeze(qAdv_neg(:,:,:,t1))'); shading flat; 
%colorbar;
caxis([0 1]*1e-13);
xlim([0 40]);
title('$-\nabla \cdot (\mathbf{u} q)$','interpreter','latex')
colormap(sb3,rbmap); %cmocean('balance','pivot',0);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
fill(x2, inBetween, c_ice,'EdgeColor','none');
set(gca,'fontsize',fontsize)
sb4 = subplot(4,2,4); hold on;
pcolor(YY/m1km,ZZ,squeeze(qAdv_neg(:,:,:,t2))'); shading flat; 
colorbar;
caxis([0 1]*1e-13);
xlim([0 40]);
title('$-\nabla \cdot (\mathbf{u} q)$','interpreter','latex')
colormap(sb4,rbmap); % cmocean('balance','pivot',0);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
yticklabels([]);
fill(x2, inBetween, c_ice,'EdgeColor','none');
set(gca,'fontsize',fontsize)


sb5 = subplot(4,2,5); hold on;
pcolor(YY/m1km,ZZ,squeeze(qF_neg(:,:,:,t1))'); shading flat; 
%colorbar;
caxis([0 1]*1e-13);
xlim([0 40]);
title('$\nabla \cdot (\nabla b \times \mathbf{F})$','interpreter','latex')
colormap(sb5,rbmap); %cmocean('balance','pivot',0);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
fill(x2, inBetween, c_ice,'EdgeColor','none');
set(gca,'fontsize',fontsize)
sb6 = subplot(4,2,6); hold on;
pcolor(YY/m1km,ZZ,squeeze(qF_neg(:,:,:,t2))'); shading flat; 
%colorbar;
caxis([0 1]*1e-13);
xlim([0 40]);
title('$\nabla \cdot (\nabla b \times \mathbf{F})$','interpreter','latex')
colormap(sb6,rbmap); %cmocean('balance','pivot',0);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
yticklabels([]);
fill(x2, inBetween, c_ice,'EdgeColor','none');
set(gca,'fontsize',fontsize)

sb7 = subplot(4,2,7); hold on;
pcolor(YY/m1km,ZZ,squeeze(qB_neg(:,:,:,t1))'); shading flat; 
%colorbar;
caxis([0 1]*1e-13);
xlim([0 40]);
title('$\nabla \cdot (\omega_a \frac{Db}{Dt})$','interpreter','latex')
colormap(sb7,rbmap); %cmocean('balance','pivot',0);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
fill(x2, inBetween, c_ice,'EdgeColor','none');
set(gca,'fontsize',fontsize)
sb8 = subplot(4,2,8); hold on;
pcolor(YY/m1km,ZZ,squeeze(qB_neg(:,:,:,t2))'); shading flat; 
%colorbar;
caxis([0 1]*1e-13);
xlim([0 40]);
title('$\nabla \cdot (\omega_a \frac{Db}{Dt})$','interpreter','latex')
colormap(sb8,rbmap); %cmocean('balance','pivot',0);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
yticklabels([]);
fill(x2, inBetween, c_ice,'EdgeColor','none');
set(gca,'fontsize',fontsize)

% subplot(5,2,9); hold on;
% pcolor(YY/m1km,ZZ,squeeze(qB_fComp_neg(:,:,:,t1))'); shading flat; 
% colorbar;
% caxis([0 1]*1e-13);
% xlim([0 40]);
% title('$\nabla \cdot (\omega_a \frac{Db}{Dt})$','interpreter','latex')
% cmocean('balance','pivot',0);
% set(gca,'Box','on','linewidth',2)
% xlim([xlim_min xlim_max])
% ylim([ylim_min ylim_max])
% yticklabels([]);
% fill(x2, inBetween, c_ice,'EdgeColor','none');
% subplot(5,2,10); hold on;
% pcolor(YY/m1km,ZZ,squeeze(qB_geoComp_neg(:,:,:,t1))'); shading flat; 
% colorbar;
% caxis([0 1]*1e-13);
% xlim([0 40]);
% title('$\nabla \cdot (\omega_a \frac{Db}{Dt})$','interpreter','latex')
% cmocean('balance','pivot',0);
% set(gca,'Box','on','linewidth',2)
% xlim([xlim_min xlim_max])
% ylim([ylim_min ylim_max])
% yticklabels([]);
% fill(x2, inBetween, c_ice,'EdgeColor','none');

%colormap('redblue')














function Aneg = negVal(A)
    Ny = size(A,2);
    Nr = size(A,3);
    Nt = size(A,4);
    for j = 1:Ny
        for k = 1:Nr
            for n = 1:Nt
                if A(1,j,k,n)>0
                    Aneg(1,j,k,n) = A(1,j,k,n);
                else
                    Aneg(1,j,k,n) = NaN;
                end
            end
        end
    end
end


function Ahov = plume_mean(AA,Nd,Ns,icetopo,zz,dz)
Ny = size(AA,2);
Nt = size(AA,4);
AA_= zeros(Ny-1,Nd,Nt);
for j = 1:Ny-1
    k = find(zz==floor(icetopo(j))+(dz)/2);
    AA_(j,:,:) = AA(1,j,k:k+(Nd-1),:);
end
Ahov = squeeze(nanmean(AA_(1:Ns,:,:),1));
end

%end

