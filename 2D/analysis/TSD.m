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
framepos = [0 scrsz(4) scrsz(3)/1.5 scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(10);
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

for n = [t1, t2]
    sigma(:,:,n) = densjmd95(SS(:,:,n),TT(:,:,n),0);
end
iceNaN(sigma,Nx,Ny,Nr,icetopo,zz);

% SUBPLOT 1
    subplot('Position',[0.1400    0.6800    0.3655    0.280])
    hold on;
    pcolor(YY/m1km,ZZ,TT(:,:,t1)');
    shading interp;
    caxis([-1.9 0.5]);
    cmocean('ice');
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    xticklabels([]);
    set(gca,'Box','on','linewidth',2)
    %title('T (^oC)')
    xlim(xlims); ylim(ylims);
    xticklabels([]);
    ylabel('z (m)')

% SUBPLOT 2
    subplot('Position',[0.1400    0.38    0.3655    0.280]);
    %%% salinity
    pcolor(YY/m1km,ZZ,SS(:,:,t1)');
    shading interp;
    %colorbar('box','on','linewidth',1.5);
    caxis([34 34.6]);
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    cmocean('-haline');
    set(gca,'fontsize',22)
    xticklabels([]);
    set(gca,'Box','on','linewidth',2)
    %title('S (psu)')
    xlim(xlims); ylim(ylims);
    ylabel('z (m)')
   
% SUBPLOT 3    
    subplot('Position',[0.1400    0.0792    0.3655    0.280])
    %%% Density
    pcolor(YY/m1km,ZZ,sigma(:,:,t1)');
    shading flat;
    %colorbar('box','on','linewidth',1.5);
    cmocean('dense');
    caxis([1027.3 1027.8]);
    %caxis([1026.7 1027.4]);
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    %xticklabels([]); yticklabels([]);
    set(gca,'Box','on','linewidth',2)
    %title('\sigma (kg m^{-3})')
    xlim(xlims); ylim(ylims);
    ylabel('z (m)'); xlabel('y (km)')

% SUBPLOT 4
    subplot('position',[0.5250    0.6800    0.3655    0.280]) 
    pcolor(YY/m1km,ZZ,(TT(:,:,t2))');
    shading interp;
    colorbar('box','on','linewidth',1.5,'Position',[0.9,0.68,0.016,0.28]);
    caxis([-1.9 0.5]);
    cmocean('ice');
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    set(gca,'Box','on','linewidth',2)
    xlim(xlims); ylim(ylims);
    yticklabels([]); xticklabels([]); 

% SUBPLOT 5
    subplot('position',[0.5250    0.38    0.3655    0.280]);
    hold on;
    pcolor(YY/m1km,ZZ,SS(:,:,t2)');
    shading interp;
    colorbar('box','on','linewidth',1.5,'Position',[0.9,0.38,0.016,0.28]);
    caxis([34 34.6]);
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    cmocean('-haline');
    set(gca,'fontsize',22)
    yticklabels([]); xticklabels([]); 
    set(gca,'Box','on','linewidth',2)
    xlim(xlims); ylim(ylims);


% SUBPLOT 6
    subplot('position',[0.5250    0.0792    0.3655    0.280])
    %%% vorticity
    pcolor(YY/m1km,ZZ,sigma(:,:,t2)');
    shading interp;
    colorbar('box','on','linewidth',1.5,'Position',[0.9,0.078,0.016,0.28]);
    cmocean('dense');
    caxis([1027.3 1027.8]);
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    set(gca,'fontsize',22)
    xlabel('y (km)');
    yticklabels([]);
    set(gca,'Box','on','linewidth',2)
    %title('\zeta /|f|')
    xlim(xlims); ylim(ylims);


    %txt = sprintf(['t = ',num2str(dumpIters(n)/t1day*deltaT),' days']);
    %set(h_annot,'String',txt)
    %drawnow
%     annotation('textbox',...
%     [0.78 0.08 0.3 0.15],...
%     'String',['t = ',num2str(n),'days'],...
%     'FontSize',24,...
%     'FontName','Arial',...
%     'FontWeight','Bold',...
%     'EdgeColor','none',...
%     'BackgroundColor','none');



