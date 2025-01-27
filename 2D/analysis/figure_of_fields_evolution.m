%%% Figure for the evolution of T,S,rho,U, psi and vort fields

m1km = 1000;
tstamp = [8,30,180];

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

sigmatheta = zeros(Nx,Ny,Nr,Nt);
for n = 1:Nt
    SS_ = squeeze(S(1,:,:,n));
    TT_ = squeeze(T(1,:,:,n));
    theta(:,:) = densjmd95(SS_,TT_,0);
    sigmatheta(1,:,:,n) = theta - 1000;
end


%%% Set figure position
scrsz = get(0,'ScreenSize');
fontsize = 26;
framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(10);
set(handle,'Position',framepos);


sp_col = length(tstamp);



for i=1:sp_col
    % plotting temperature
    %subplot(sp_col,3,3*i-2)
    TT = iceNaN(T,Nx,Ny,Nr,icetopo,zz); 
    Tmin = -1.9;
    Tmax = 0.5;
    levels = [Tmin :0.05: Tmax];
    sb = subplot('position',[0.15 0.68-(i-1)*(0.27+0.015)   0.25    0.27]);
    contourf(YY/m1km,ZZ,squeeze(TT(1,:,:,tstamp(i)))',levels,'Linecolor','none')
    shading interp;
    cmocean('ice')
    caxis([Tmin Tmax]) 
    sb.Box = 'on';
    sb.LineWidth = 1.5;
    set(gca,'fontsize',20)
    ylabel('z (m)')
    yticks([-450 -350 -250 -150 -50])
    yticklabels([-450 -350 -250 -150 -50]);
    if i ==1
        title ('T (^oC)')
    end
    if i == 3
        xlabel('y (km)')
    end
    if i ~=3
        cb = colorbar('horiz');
        cb.Position = [0.157 0.034 0.21 0.015];
        cb.Box = 'on';
        cb.LineWidth = 1.5;
        xticklabels([]);
    end
    xlim([1 40])
    ylim([-500 0])
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    
    % plotting salinity
    SS = iceNaN(S,Nx,Ny,Nr,icetopo,zz);
    Smin = 33.54;
    Smax = 34.6;
    levels = [Smin : 0.04 : Smax];
    sb = subplot('position',[0.43 0.68-(i-1)*(0.27+0.015)   0.25    0.27]);
    contourf(YY/m1km,ZZ,squeeze(SS(1,:,:,tstamp(i)))',levels,'Linecolor','none')
    shading interp;
    sb.Box = 'on';
    sb.LineWidth = 1.5;
    cmocean('-haline')
    caxis([Smin Smax]) 
    set(gca,'fontsize',20)
    yticklabels([]);
    if i == 1
        title('S (psu)')
    end
    if i ~= 3
        xticklabels([]);
    end
    if i == 3
        cb = colorbar('horiz');
        cb.Position = [0.45 0.034 0.21 0.015];
        cb.Box = 'on';
        cb.LineWidth = 1.5;
        xlabel('y (km)')
    end
    xlim([1 40])
    ylim([-500 0])
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    
    % plotting density
    sigma = iceNaN(sigmatheta,Nx,Ny,Nr,icetopo,zz);
    sig_min = 27;
    sig_max = 27.8;
    levels = [sig_min : 0.03 : sig_max];
    sb = subplot('position',[0.71 0.68-(i-1)*(0.27+0.015)   0.25    0.27]);
    contourf(YY/m1km,ZZ,squeeze(sigma(1,:,:,tstamp(i)))',levels,'Linecolor','none')
    shading interp;
    sb.Box = 'on';
    sb.LineWidth = 1.5;
    cmocean('dense')
    caxis([sig_min sig_max])
    set(gca,'fontsize',20)
    yticklabels([]);
    if i == 1
        title('\sigma (kg m^{-3})')
    end
    if i~= 3
        xticklabels([]);
    end 
    if i == 3
        %colorbar
        xlabel('y (km)')
        cb = colorbar('horiz');
        %cb.Location = 'southoutside';
        cb.Position = [0.73 0.034 0.21 0.015];
        cb.Box = 'on';
        cb.LineWidth = 1.5;
    end
    txt={['t = ',num2str(tstamp(i)),' days']};
    h_annot=annotation('textbox',...
        [0.8 0.61-(i-1)*(0.27+0.015) 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');
    xlim([1 40])
    ylim([-500 0])
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    
    
end


%% figure 2
addpath /Users/Kayliec23/Desktop/plottingKit

%%% Set figure position
scrsz = get(0,'ScreenSize');
fontsize = 26;
framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(11);
set(handle,'Position',framepos);


sp_col = length(tstamp);



for i=1:sp_col
    % plotting zonal velocity
    UU = iceNaN(U,Nx,Ny,Nr,icetopo,zz); 
    Umin = -0.2;
    Umax = 0.2;
    %levels = [Umin :0.001: Umax];
    sb = subplot('position',[0.15 0.68-(i-1)*(0.27+0.015)   0.25    0.27]);
    pcolor(YY/m1km,ZZ,squeeze(UU(1,:,:,tstamp(i)))') %,'Linecolor','none')
    shading flat;
    cmocean('balance','pivot',0)
    caxis([Umin Umax]) 
    sb.Box = 'on';
    sb.LineWidth = 1.5;
    set(gca,'fontsize',20)
    ylabel('z (m)')
    yticks([-450 -350 -250 -150 -50])
    yticklabels([-450 -350 -250 -150 -50]);
    if i ==1
        title ('U (ms^{-1})')
    end
    if i == 3
        xlabel('y (km)')
        cb = colorbar('horiz');
        cb.Position = [0.157 0.034 0.21 0.015];
        cb.Box = 'on';
        cb.LineWidth = 1.5;
    end
    if i ~=3
        xticklabels([]);
    end
    xlim([1 40])
    ylim([-500 0])
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    hold off;
    
    % plotting sreamfunction
    psi = zeros(Nx,Ny,Nr,Nt);
    for n = 1:Nt
        psi(:,:,:,n) = -nancumsum(squeeze(W(1,:,:,n)).*(delY'));
    end
    psi = iceNaN(psi,Nx,Ny,Nr,icetopo,zz);
    psi_min = -0.27;
    psi_max = 0.027;
    %levels = [psi_min : 0.001 : psi_max];
    sb = subplot('position',[0.43 0.68-(i-1)*(0.27+0.015)   0.25    0.27]);
    pcolor(YY/m1km,ZZ,squeeze(psi(1,:,:,tstamp(i)))') %,'Linecolor','none')
    shading flat;
    sb.Box = 'on';
    sb.LineWidth = 1.5;
    caxis([psi_min psi_max])
    cmocean('deep')
    set(gca,'fontsize',20)
    yticklabels([]);
    if i == 1
        title('\psi (m^2 s^{-1})')
    end
    if i ~= 3
        xticklabels([]);
    end
    if i == 3
        cb = colorbar('horiz');
        cb.Position = [0.45 0.034 0.21 0.015];
        cb.Box = 'on';
        cb.LineWidth = 1.5;
        xlabel('y (km)')
    end
    xlim([1 40])
    ylim([-500 0])
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    hold off;
    
    % plotting vorticity
    vort = ZETAz;
    vort = iceNaN(vort,Nx,Ny,Nr,icetopo,zz); 
    vort_min = -1.3;
    vort_max = 1.3;
    %levels = [vort_min : 0.001 : vort_max];
    sb = subplot('position',[0.71 0.68-(i-1)*(0.27+0.015)   0.25    0.27]);
    pcolor(YY/m1km,ZZ,squeeze(vort(1,:,:,tstamp(i))./abs(f0))') %,'Linecolor','none')
    shading flat;
    sb.Box = 'on';
    sb.LineWidth = 1.5;
    caxis([vort_min vort_max])
    cmocean('balance','pivot',0)
    set(gca,'fontsize',20)
    yticklabels([]);
    if i == 1
        title('\zeta/|f|')
    end
    if i~= 3
        xticklabels([]);
    end 
    if i == 3
        %colorbar
        xlabel('y (km)')
        cb = colorbar('horiz');
        %cb.Location = 'southoutside';
        cb.Position = [0.73 0.034 0.21 0.015];
        cb.Box = 'on';
        cb.LineWidth = 1.5;
    end
    txt={['t = ',num2str(tstamp(i)),' days']};
    h_annot=annotation('textbox',...
        [0.8 0.61-(i-1)*(0.27+0.015) 1 0.11],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'FontWeight','Bold',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','on');

    xlim([1 40])
    ylim([-500 0])
    hold on;
    plot(yy/m1km,icetopo,'k','linewidth',1.5)
    hold off;
    
    
end




