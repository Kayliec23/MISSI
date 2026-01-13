%%% anim_all
%set_paths_and_load_data;
addpath /Users/kalyliec23/Desktop/Research/MISSI
DZ_ = diff(zz,1,2); DZ = DZ_(1);
DY_ = diff(yy,1,2); DY = DY_(1);
[YY,ZZ] = meshgrid(yy,zz);
Nt = size(T,4);
dy = Ly/length(yy);
dx = Lx/length(xx);
dz = H/length(zz);g = 9.8;
rho0 = 1027;
m1km = 1000;
Omega = 2*pi*366/365/86400;
lat0 = -80; %%% Latitude at southern boundary
f0 = 2*Omega*sind(lat0); %%% Coriolis parameter

zoom = 1;
if zoom == 1
    xlims = [4 13]; %[10 20];
    ylims = [-400 -200]; %[-300 -70];
else
    xlims = [0 40];
    ylims = [-500 0];
end

%%
scrsz = get(0,'ScreenSize');
fontsize = 26;
framepos = [0 scrsz(4) scrsz(3)/2 scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(11);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');

txt='';
h_annot=annotation('textbox',...
    [0.76 0.04 0.3 0.11],...
    'String',txt,...
    'FontSize',32,...
    'FontName','Arial',...
    'FontWeight','Bold',...
    'EdgeColor','none',...
    'BackgroundColor','none');

%%% choose timestamps for figure
t1 = 3; t2 = 25; 

TT = iceNaN(T,Nx,Ny,Nr,icetopo,zz); TT = squeeze(TT);
SS = iceNaN(S,Nx,Ny,Nr,icetopo,zz); SS = squeeze(SS);  
UU = iceNaN(U,Nx,Ny,Nr,icetopo,zz); UU = squeeze(UU);  
VV = iceNaN(V,Nx,Ny,Nr,icetopo,zz); VV = squeeze(VV);  
WW = iceNaN(W,Nx,Ny,Nr,icetopo,zz); WW = squeeze(WW); 
visc = iceNaN(KPPvisc,Nx,Ny,Nr,icetopo,zz); visc = squeeze(visc); 

for n = [t1, t2]
    %%% density
    sigma(:,:,n) = densjmd95(SS(:,:,n),TT(:,:,n),0);
    %%% zonal shear
    dudz(:,1:Nr-1,n) = [UU(:,1:Nr-1,n)-UU(:,2:Nr,n)]./dz;
    dudz(:,Nr,n) = 2*dudz(:,Nr-1,n) - dudz(:,Nr-2,n);
    %%% PV, Ri and b
    [Q(:,:,n), Ri(:,:,n), b(:,:,n)] = calculateQ(S(1,:,:,n),T(1,:,:,n),U(1,:,:,n),V(1,:,:,n),...
        W(1,:,:,n),Nx,Ny,Nr,Nt,dy,dz,g,rho0,f0);
    %%% bz
    dbdz(:,1:Nr-1,n) = [b(:,1:Nr-1,n)-b(:,2:Nr,n)]./dz;
    dbdz(:,Nr,n) = 2*dbdz(:,Nr-1,n) - dbdz(:,Nr-2,n);
end
sigma = iceNaN(sigma,Nx,Ny,Nr,icetopo,zz);
q = squeeze(Q); 
fq = iceNaN(f0*q,Nx,Ny,Nr,icetopo,zz);
bz = iceNaN(dbdz,Nx,Ny,Nr,icetopo,zz);
uz = iceNaN(dudz,Nx,Ny,Nr,icetopo,zz);

% SUBPLOT 1
    ax1 = subplot('Position',[0.1600    0.8200    0.3455    0.160]);
    pcolor(ax1, YY/m1km,ZZ,fq(:,:,t1)');
    shading(ax1,'flat');
    caxis(ax1, [-5 5]*1e-13);
    mymap = colormap('redblue');
    colormap(ax1, mymap)
    hold on;
    contour(YY/m1km,ZZ,(bz(:,:,t1))',linspace(1e-5,2e-5,10),'k')
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    set(gca,'Box','on','linewidth',2)
    xticklabels([]);
    yticks([-380 -300 -220])
    xlim(xlims); ylim(ylims);
    ylabel('z (m)');
    hold off;


% SUBPLOT 6
    ax6 = subplot('position',[0.5250    0.8200    0.3455    0.160]) ;
    pcolor(ax6, YY/m1km,ZZ,fq(:,:,t2)');
    shading(ax6, 'flat');
    colorbar('box','on','linewidth',1.5,'Position',[0.88,0.83,0.016,0.13]);
    caxis(ax6, [-5 5]*1e-13);
    mymap = colormap('redblue');
    colormap(ax6, mymap)
    hold on;
    contour(YY/m1km,ZZ,(bz(:,:,t2))',linspace(1e-5,2e-5,10),'k')
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    yticklabels([]); xticklabels([]);
    set(gca,'Box','on','linewidth',2)
    xlim(xlims); ylim(ylims);
    hold off;


% SUBPLOT 2
    ax2 = subplot('Position',[0.1600    0.64    0.3455    0.160]);
    %%% salinity
    pcolor(ax2,YY/m1km,ZZ,Ri(:,:,t1)');
    shading( ax2, 'flat');
    caxis(ax2, [0 1]);
    hold on;
    contour(YY/m1km,ZZ,(bz(:,:,t1))',linspace(1e-5,2e-5,10),'k')
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    mymap = cmocean('-deep');
    colormap(ax2, mymap)
    set(gca,'fontsize',22)
    xticklabels([]);
    yticks([-380 -300 -220])
    set(gca,'Box','on','linewidth',2)
    xlim(xlims); ylim(ylims);
    ylabel('z (m)')
    hold off;


% SUBPLOT 3    
    ax3 = subplot('Position',[0.1600    0.46    0.3455    0.160]);
    hold on;
    pcolor(ax3, YY/m1km,ZZ,visc(:,:,t1)');
    shading(ax3, 'flat');
    caxis(ax3, [0 1]*1e-3);
    mymap = cmocean('matter');
    colormap(ax3, mymap)
    hold on;
    contour(YY/m1km,ZZ,(bz(:,:,t1))',linspace(1e-5,2e-5,10),'k')
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(ax3,'fontsize',22)
    set(ax3,'Box','on','linewidth',2)
    xlim(xlims); ylim(ylims);
    ylabel('z (m)');
    yticks([-380 -300 -220])
    xticklabels([]);
    hold off;

% SUBPLOT 4
    ax4 = subplot('position',[0.1600    0.28    0.3455    0.160]);
    pcolor(ax4, YY/m1km,ZZ,(bz(:,:,t1))');
    shading flat;
    caxis([-3 3]*1e-4)
    mymap = cmocean('balance','pivot',0);
    colormap(ax4, mymap)
    hold on;
    contour(YY/m1km,ZZ,(bz(:,:,t1))',linspace(1e-5,2e-5,10),'k')
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    set(gca,'Box','on','linewidth',2)
    xlim(xlims); ylim(ylims);
    xticklabels([]);
    yticks([-380 -300 -220])
    ylabel('z (m)')
    hold off;
  
% SUBPLOT 5
    ax5 = subplot('position',[0.1600    0.1    0.3455    0.160]);
    pcolor(ax5, YY/m1km,ZZ,(uz(:,:,t1))');
    shading flat;
    caxis([-4 4]*1e-2)
    mymap = cmocean('balance','pivot',0);
    colormap(ax5, mymap)
    hold on;
    contour(YY/m1km,ZZ,(bz(:,:,t1))',linspace(1e-5,2e-5,10),'k')
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    set(gca,'Box','on','linewidth',2)
    xlim(xlims); ylim(ylims);
    xlabel('y (km)'); ylabel('z (m)');
    yticks([-380 -300 -220])
    hold off;
 

% SUBPLOT 7
    ax7 = subplot('position',[0.5250    0.64    0.3455    0.160]);
    hold on;
    pcolor(ax7, YY/m1km,ZZ,Ri(:,:,t2)');
    shading flat;
    colorbar('box','on','linewidth',1.5,'Position',[0.88,0.65,0.016,0.13]);
    caxis([0 1]);
    mymap = cmocean('-deep');
    colormap(ax7, mymap)
    hold on;
    contour(YY/m1km,ZZ,(bz(:,:,t2))',linspace(1e-5,2e-5,10),'k')
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    set(gca,'Box','on','linewidth',2)
    xlim(xlims); ylim(ylims);
    xticklabels([]); yticklabels([]);
    hold off;


% SUBPLOT 8
    ax8 = subplot('position',[0.5250    0.46    0.3455    0.160]);
    pcolor(ax8, YY/m1km,ZZ,(visc(:,:,t2))');
    shading flat;
    colorbar('box','on','linewidth',1.5,'Position',[0.88,0.47,0.016,0.13]);
    mymap = cmocean('matter');
    colormap(ax8, mymap)
    hold on;
    contour(YY/m1km,ZZ,(bz(:,:,t2))',linspace(1e-5,2e-5,10),'k')
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    set(gca,'Box','on','linewidth',2)
    xlim(xlims); ylim(ylims);
    yticklabels([]); xticklabels([]);
    hold off;
 
% SUBPLOT 9
    ax9 = subplot('position',[0.5250    0.28    0.3455    0.160]);
    pcolor(ax9, YY/m1km,ZZ,(bz(:,:,t2))');
    shading flat;
    colorbar('box','on','linewidth',1.5,'Position',[0.88,0.29,0.016,0.13]);
    caxis([-3 3]*1e-4)
    mymap = cmocean('balance','pivot',0);
    colormap(ax9, mymap)
    hold on;
    contour(YY/m1km,ZZ,(bz(:,:,t2))',linspace(1e-5,2e-5,10),'k')
    contour(YY/m1km,ZZ,(bz(:,:,t2))')
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    set(gca,'Box','on','linewidth',2)
    xlim(xlims); ylim(ylims);
    yticklabels([]); xticklabels([]); 
    hold off;

 % SUBPLOT 10
    ax10 = subplot('position',[0.5250    0.1    0.3455    0.160]);
    pcolor(ax10, YY/m1km,ZZ,(uz(:,:,t2))');
    shading flat;
    colorbar('box','on','linewidth',1.5,'Position',[0.88,0.11,0.016,0.13]);
    caxis([-4 4]*1e-2)
    mymap = cmocean('balance','pivot',0);
    colormap(ax10, mymap)
    hold on;
    contour(YY/m1km,ZZ,(bz(:,:,t2))',linspace(1e-5,2e-5,10),'k')
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    set(gca,'Box','on','linewidth',2)
    xlim(xlims); ylim(ylims);
    yticklabels([]); 
    xlabel('y (km)')
    hold off;



function [Q,Ri,b] = calculateQ(S,T,U,V,W,Nx,Ny,Nr,Nt,dy,dz,g,rho0,f0)
% Calculate potential density (at grid center)
sigmatheta = zeros(Nx,Ny,Nr);
    SS = squeeze(S(1,:,:));
    TT = squeeze(T(1,:,:));
    theta(:,:) = densjmd95(SS,TT,0);
    sigmatheta(1:Nx,:,:) = theta - 1000;


% Calculate ZETAy
dUdz = NaN(Nx,Ny,Nr); dUdz(:,:,2:Nr) = (U(:,:,[1:Nr-1]) - U(:,:,[2:Nr]))./dz; % derivative wrt z, shifts to top face center
dUdz(:,1:Ny-1,:) = (dUdz(:,[2:Ny],:)+dUdz(:,[1:Ny-1],:))/2; % average to top right corner
dUdz(:,:,1) = 2*dUdz(:,:,2)-dUdz(:,:,3); dUdz(:,Ny,:) = 2*dUdz(:,Ny-1,:)-dUdz(:,Ny-2,:); % add boundaries
ZETAy = dUdz;

% Calcualte ZETAz
dUdy = NaN(Nx,Ny,Nr); dUdy(:,1:Ny-1,:) = diff(U(:,:,:),1,2)./dy; % derivative wrt y, shifts to right face center
dUdy(:,:,2:Nr) = (dUdy(:,:,[2:Nr])+dUdy(:,:,[1:Nr-1]))/2; % average to top right corner
dUdy(:,:,1) = 2*dUdy(:,:,2)-dUdy(:,:,3); dUdy(:,Ny,:) = 2*dUdy(:,Ny-1,:)-dUdy(:,Ny-2,:); % add boundaries
ZETAz = -dUdy;

% Density derivatives
dRHOdy = NaN(Nx,Ny,Nr); dRHOdy(:,1:Ny-1,:) = diff(sigmatheta,1,2)./dy; % derivative wrt y, shifts to right face center
dRHOdy(:,:,2:Nr) = (dRHOdy(:,:,[2:Nr])+dRHOdy(:,:,[1:Nr-1]))/2; % move to top right corner
dRHOdy(:,Ny,:) = 2*dRHOdy(:,Ny-1,:)-dRHOdy(:,Ny-2,:); dRHOdy(:,:,1) = 2*dRHOdy(:,:,2)-dRHOdy(:,:,3); % add boundaries
dRHOdz = NaN(Nx,Ny,Nr); dRHOdz(:,:,2:Nr) = (sigmatheta(:,:,[1:Nr-1]) - sigmatheta(:,:,[2:Nr]))./dz; % derivative wrt z, shifts to top face center
dRHOdz(:,1:Ny-1,:) = (dRHOdz(:,[2:Ny],:)+dRHOdz(:,[1:Ny-1],:))/2; % move to top right corner
dRHOdz(:,:,1) = 2*dRHOdz(:,:,2)-dRHOdz(:,:,3); dRHOdz(:,Ny,:) = 2*dRHOdz(:,Ny-1,:)-dRHOdz(:,Ny-2,:); % add boundaries

dbdy = -g/rho0.*dRHOdy;
dbdz = -g/rho0.*dRHOdz;
b = -g/rho0.*sigmatheta;

% Q and its symmetric, convecitve and centrifugal components
%Qbc = dbdy.*ZETAy; % symmetric
%Qvert = dbdz.*f0 + dbdz.*ZETAz; % convective
Q = dbdy.*ZETAy + (f0 + ZETAz).*dbdz; % total

% gradient Richardson number
%Ri = f0.^2*(-g/rho0*dRHOdz)./((-g/rho0*dRHOdy).^2);
Ri = (-g/rho0*dRHOdz)./(abs(dUdz).^2);
Ri(Ri==inf) = NaN; Ri(Ri==-inf) = NaN;

end




