addpath /Users/Kayliec23/Desktop/analysis
addpath /Users/Kayliec23/Desktop/MISSI/

% Constants
Nx = length(delX);
Ny = length(delY);
Nr = length(delR);
Nt = size(T,4)-1;
%%% Domain dimensions
Lx = sum(delX);
Ly = sum(delY);
H = sum(delR);
%%% Gridpoint locations are at the centre of each grid cell
xx = cumsum((delX + [0 delX(1:Nx-1)])/2)-Lx/2;
yy = cumsum((delY + [0 delY(1:Ny-1)])/2);
zz = -cumsum((delR + [0 delR(1:Nr-1)])/2);
%%% grid size
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


%%% Calculate Q and take time derivative
[q, qvert, qbc, Rii, Rii_pos] = calculateQ(S,T,U,V,W,Nx,Ny,Nr,Nt,dy,dz,g,rho0,f0);
q = iceNaN(q,Nx,Ny,Nr,icetopo,zz);
Rii = iceNaN(Rii,Nx,Ny,Nr,icetopo,zz);
Rii_pos = iceNaN(Rii_pos,Nx,Ny,Nr,icetopo,zz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting fq terms and Ri %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




[YY,ZZ] = meshgrid(yy,zz);
scrsz = get(0,'ScreenSize');
fontsize = 28;
framepos = [0 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(13);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');

hplot = 0.27;
wplot = 0.14;
t1 = 8; t2 = 20; t3 = 60;
xlim_min = 6; %1;
xlim_max = 15; %Ny*dy(1)/m1km;
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
pcolor(YY/m1km,ZZ,squeeze(f0*q(1,:,:,t1))');
shading flat;
caxis([-5 5]*1e-13);
hold on;
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xticklabels([])
ylabel('z (m)')
set(gca,'Box','on','linewidth',2)
title('fq')
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
yticks([-290 -260 -230 -200])
yticklabels({'-290','-260','-230','-200'})
fill(x2, inBetween, c_ice,'EdgeColor','none');



subplot('position',[0.08    0.4034    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*q(1,:,:,t2))');
shading flat;
caxis([-5 5]*1e-13);
hold on;
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
yticks([-290 -260 -230 -200])
yticklabels({'-290','-260','-230','-200'});
fill(x2, inBetween, c_ice,'EdgeColor','none');


subplot('position',[0.08    0.1195    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*q(1,:,:,t3))');
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
xticks([9 12 15 18])
xticklabels({'9','12','15','18'})
yticks([-290 -260 -230 -200])
yticklabels({'-290','-260','-230','-200'});
fill(x2, inBetween, c_ice,'EdgeColor','none');

%%% plotting q1
subplot('position',[0.235    0.6874    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*qvert(1,:,:,t1))');
shading flat;
caxis([-5 5]*1e-13);
hold on;
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xticklabels([])
yticklabels([]);
set(gca,'Box','on','linewidth',2)
title('fq_{vert}')
xlim([xlim_min xlim_max]);
ylim([ylim_min ylim_max]);
fill(x2, inBetween, c_ice,'EdgeColor','none');


subplot('position',[0.235    0.4034    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*qvert(1,:,:,t2))');
shading flat;
caxis([-5 5]*1e-13);
hold on;
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xticklabels([])
yticklabels([]);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max]);
ylim([ylim_min ylim_max]);
fill(x2, inBetween, c_ice,'EdgeColor','none');


subplot('position',[0.235    0.1195    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*qvert(1,:,:,t3))');
shading flat;
caxis([-5 5]*1e-13);
hold on;
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xlabel('y (km)');
yticklabels([]);
set(gca,'Box','on','linewidth',2)
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
xticks([9 12 15 18])
xticklabels({'9','12','15','18'});
fill(x2, inBetween, c_ice,'EdgeColor','none');


%%% plotting q2
subplot('position',[0.39    0.6874    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*qbc(1,:,:,t1))');
shading flat;
caxis([-5 5]*1e-13);
hold on;
mymap = colormap('redblue');
colormap(mymap);
set(gca,'fontsize',22)
xticklabels([])
set(gca,'Box','on','linewidth',2)
title('fq_{bc}')
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
yticklabels([]);
fill(x2, inBetween, c_ice,'EdgeColor','none');


subplot('position',[0.39    0.4034    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*qbc(1,:,:,t2))');
shading flat;
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


subplot('position',[0.39    0.1195    wplot   hplot])
pcolor(YY/m1km,ZZ,squeeze(f0*qbc(1,:,:,t3))');
shading flat;
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


%%% plotting Ri
h = subplot('position',[0.546    0.6874    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(Rii_pos(:,:,:,t1))');
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


h = subplot('position',[0.546    0.4034    wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(Rii_pos(:,:,:,t2))');
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


h = subplot('position',[0.546    0.1195   wplot   hplot]);
pcolor(YY/m1km,ZZ,squeeze(Rii_pos(:,:,:,t3))');
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











function [Q,Qvert,Qbc,Ri, Rii_pos] = calculateQ(S,T,U,V,W,Nx,Ny,Nr,Nt,dy,dz,g,rho0,f0)
% Calculate potential density (at grid center)
sigmatheta = zeros(Nx,Ny,Nr,Nt);
for n = 1:Nt
    SS = squeeze(S(1,:,:,n));
    TT = squeeze(T(1,:,:,n));
    theta(:,:) = densjmd95(SS,TT,0);
    sigmatheta(1,:,:,n) = theta - 1000;
end


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
Qbc = dbdy.*ZETAy; % symmetric
Qvert = dbdz.*f0 + dbdz.*ZETAz; % convective
Q = dbdy.*ZETAy + (f0 + ZETAz).*dbdz; % total

% gradient Richardson number
%Ri = f0.^2*(-g/rho0*dRHOdz)./((-g/rho0*dRHOdy).^2);
Ri = (-g/rho0*dRHOdz)./(abs(dUdz).^2);
Ri(Ri==inf) = NaN; Ri(Ri==-inf) = NaN;

Rii_pos = Ri;
for j = 1:Ny
    for k =1:Nr
        for n = 1:Nt
            if f0*Q(1,j,k,n)>0 % Ri for positive fq
                Rii_pos(:,j,k,n) = NaN;
            end
        end
    end
end

end
