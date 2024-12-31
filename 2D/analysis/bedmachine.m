clear all;

addpath /Users/kalyliec23/Desktop/Research/MISSI/BedMachine
addpath /Users/kalyliec23/Desktop/Research/MISSI/BedMachine/github_repo-2
addpath /Users/kalyliec23/Desktop/Research/MISSI/BedMachine/github_repo-3
load('full_info.mat')




df = shelves.Pine_Island; %Thwaites, Pine Island, Crosson, Dotson
df(df>0) = NaN;
%df = smoothdata2(df,'movmean',5,'omitnan');
lon = squeeze(coords.Pine_Island(1,:,:));
lat = squeeze(coords.Pine_Island(2,:,:));
Nx = size(df,1);
Ny = size(df,2);
dx = 500;
dy = 500;
xx = linspace(0,dx*(Nx-1),Nx);
yy = linspace(0,dy*(Ny-1),Ny);
[XX YY] = meshgrid(xx,yy);
g = 9.8;
gamma = 2.5e-5;
ff = -1.4e-4;
Hpyc = 10;
m1km = 1000;

%%% first derivative - slope
s_x = zeros(Nx,Ny); s_y = zeros(Nx,Ny); s = zeros(Nx,Ny);
s_x(1:Nx-1,:) = (df(2:Nx,:)-df(1:Nx-1,:))./dx; %s_x(1:Nx-1,:) = (df(2:Nx,:)+df(1:Nx-1,:))/2; s_x(Nx,:) = 2*s_x(Nx-1,:) - s_x(Nx-2,:);
s_y(:,1:Ny-1) = (df(:,2:Ny)-df(:,1:Ny-1))./dy; %s_y(:,1:Ny-1) = (df(:,2:Ny)+df(:,1:Ny-1))/2; s_y(:,Ny) = 2*s_y(:,Ny-1) - s_y(:,Ny-2);
s = sqrt(s_x.^2 + s_y.^2);

% %%% second derivative - curvature
% s_xx = zeros(Nx,Ny); s_yy = zeros(Nx,Ny); s_xy = zeros(Nx,Ny); s_yx = zeros(Nx,Ny);
% s_xx(2:Nx-1,:) = (df(3:Nx,:)-2*df(2:Nx-1,:)+df(1:Nx-2,:)) ./ (dx.^2); %%% wrt x
% %s_xx(1,:) = 2*s_xx(2,:)-s_xx(3,:); s_xx(Nx,:) = 2*s_xx(Nx-1,:) - s_xx(Nx-2,:); % boundaries
% s_yy(:,2:Ny-1) = (df(:,3:Ny)-2*df(:,2:Ny-1)+df(:,1:Ny-2)) ./ (dy.^2); %%% wrt y
% %s_yy(:,1) = 2*s_yy(:,2)-s_yy(:,3); s_yy(:,Ny) = 2*s_yy(:,Ny-1) - s_yy(:,Ny-2); % boundaries
% s_xy(:,2:Ny) = (s_x(:,2:Ny)-s_x(:,1:Ny-1))./dy; %%% wrt x then y
% ds = sqrt(s_xx.^2+s_yy.^2);
% H = s_xx.*s_yy - s_xy.^2;


%dtheta_c = ff^2/g/gamma./(s.^2/Hml + ds);
dtheta_c = (ff^2*Hpyc) ./ (g*gamma*s.^2);


% ds_min = NaN(Nx,Ny);
% for i = 1:Nx
%     for j = 1:Ny
%         if H(i,j)>0
%             if s_xx(i,j) > 0
%                 ds_min(i,j) = dtheta_c(i,j);
%             end
%         end
%     end
% end



%%%%%% plots %%%%%%

% subplot(1,2,1)
% contourf(lat,lon,smoothdata2(df,'movmean',5,'omitnan'),[min(min(df)):100:0]);
% c = colorbar;
% cmocean('deep')
% set(gca,'fontsize',28)
% set(gca, 'XDir','reverse')
% title('Depth (m)')
% ax = gca;
% ax.Box = 'on';
% ax.LineWidth = 2;
%
dtheta_c_plot = dtheta_c;
dtheta_c_plot(dtheta_c_plot>8) = 10;


%subplot(3,2,5)
cn = contourf(lat,lon,smoothdata2(dtheta_c_plot,'movmean',5,'omitnan'),[0:1:8],'edgecolor','none');
%cn = contourf(lat,lon,dtheta_c_plot,[0:1:8]);
%pcolor(YY/m1km,XX/m1km,smoothdata2(dtheta_c,'gaussian',5,'omitnan')'); shading interp;
%c = colorbar;
cmocean('-curl','pivot',4)
set(gca,'fontsize',28)
set(gca, 'XDir','reverse')
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
title('Pine Island','fontsize',18)
xticklabels([])
yticklabels([])
sgtitle('\theta (^oC)','fontsize',36)



%%
%%% make plot for paper

addpath /Users/kalyliec23/Desktop/Research/MISSI/BedMachine
addpath /Users/kalyliec23/Desktop/Research/MISSI/BedMachine/github_repo-2
addpath /Users/kalyliec23/Desktop/Research/MISSI/BedMachine/github_repo-3
load('full_info.mat')


g = 9.8;
gamma = 2.5e-5;
ff = -1.4e-4;
Hpyc = 10;
m1km = 1000;

%%% masks
[lat_coast,lon_coast] = bedmachine_data('coast','geo');
lon_coast_plot = lon_coast;
[lat_gl,lon_gl] = bedmachine_data('gl','geo');
lon_gl_plot = lon_gl;



figure;

% Amundsen
%subplot('position',[0.1065    0.0333    0.1622    0.8821]);
ax1 = subplot('position',[0.1065    0.0955    0.1622    0.8308]);
%ax1 = axes;
hold on;
caxis([0 8]);
cm1 = cmocean('-curl','pivot',4);
colormap(ax1,cm1);


amundsen_shelves = {'Pine_Island','Crosson','Dotson','Thwaites'};
fn = fieldnames(shelves);
for i = 1:length(amundsen_shelves)
    df = shelves.(amundsen_shelves{i});
    df(df>0) = NaN;
    lon = squeeze(coords.(amundsen_shelves{i})(1,:,:));
    lon(lon==Inf) = NaN;
    lat = squeeze(coords.(amundsen_shelves{i})(2,:,:));
    lat(lat==Inf) = NaN;

    Nx = size(df,1);
    Ny = size(df,2);
    dx = 500;
    dy = 500;
    xx = linspace(0,dx*(Nx-1),Nx);
    yy = linspace(0,dy*(Ny-1),Ny);
    [XX YY] = meshgrid(xx,yy);

    s_x = zeros(Nx,Ny); s_y = zeros(Nx,Ny); s = zeros(Nx,Ny);
    s_x(1:Nx-1,:) = (df(2:Nx,:)-df(1:Nx-1,:))./dx; %s_x(1:Nx-1,:) = (df(2:Nx,:)+df(1:Nx-1,:))/2; s_x(Nx,:) = 2*s_x(Nx-1,:) - s_x(Nx-2,:);
    s_y(:,1:Ny-1) = (df(:,2:Ny)-df(:,1:Ny-1))./dy; %s_y(:,1:Ny-1) = (df(:,2:Ny)+df(:,1:Ny-1))/2; s_y(:,Ny) = 2*s_y(:,Ny-1) - s_y(:,Ny-2);
    s = sqrt(s_x.^2 + s_y.^2);

    dtheta_c = (ff^2*Hpyc) ./ (g*gamma*s.^2);
    dtheta_c_plot = dtheta_c;
    dtheta_c_plot(dtheta_c_plot>8) = 10;
    dtheta_smooth = smoothdata2(dtheta_c_plot,'movmean',5,'omitnan');

    pcolor(ax1,lat,lon,dtheta_smooth);
    shading interp;
    set(ax1,'fontsize',28)
    set(ax1, 'XDir','reverse')
    colormap(ax1,cm1);

    lat_min = -75.45; lat_max = -74.1;
    lon_min = -114.5; lon_max = -99;
    xlim([lat_min lat_max])
    ylim([lon_min lon_max])

    xlabel('Lat')
ylabel('Lon')


end
xticks(ax1,[-75,-74.1]);
xticklabels(ax1,{'-75','-74'});
yticks(ax1,[-110,-105,-100]);
yticklabels(ax1,{'-110','-105','-100'});

set(gca,'color',[0.95 0.95 0.95]);

% %%% create masks
[lati,loni] = meshgrid([lat_min:0.01:lat_max],[lon_min:0.1:lon_max]);
[Z_mask,Lat_mask,Lon_mask] = bedmachine_data('mask',lati,loni,'geo');
ocean_mask = Z_mask;
ocean_mask(ocean_mask~=0)=NaN;
%ocean_mask(ocean_mask==0)=10;

%%% plot mask
ax2 = axes('position',[0.1065    0.0955    0.1622    0.8308]);
hold on;
colormap(ax2,[0 0.3 1])
p = pcolor(ax2,Lat_mask,Lon_mask,ocean_mask);
shading interp;
alpha(p,0.3)
set(ax2, 'XDir','reverse')
%Turn off ax2 so you can see the first plot
ax2.Visible = 'off';
% set(ax2,'Color','None');
ax2.XTick = [];
ax2.YTick = [];

ylim([lon_min lon_max]);
xlim([lat_min lat_max]);

ax1.Box = 'on';
ax1.LineWidth = 2;
ax2.Box = 'on';
ax2.LineWidth = 2;


% Filchner-Ronne
%subplot('Position',[0.2809    0.4856    0.3184    0.4304])
ax1 = subplot('Position',[0.3497    0.5443    0.2973    0.3832]);
hold on;
caxis([0 8]);
cm1 = cmocean('-curl','pivot',4);
colormap(ax1,cm1);

Ronne_shelves = {'Ronne','Filchner'};
fn = fieldnames(shelves);

for i = 1:length(Ronne_shelves)
    df = shelves.(Ronne_shelves{i});
    df(df>0) = NaN;
    lon = squeeze(coords.(Ronne_shelves{i})(1,:,:));
    lon(lon==Inf) = NaN;
    lat = squeeze(coords.(Ronne_shelves{i})(2,:,:));
    lat(lat==Inf) = NaN;
    %%% masks
    [lat_coast,lon_coast] = bedmachine_data('coast','geo');
    lon_coast_plot = lon_coast;
    [lat_gl,lon_gl] = bedmachine_data('gl','geo');
    lon_gl_plot = lon_gl;

    %%% adjusting longitude for plot
    lon_plot = lon;
    lon_plot(lon_plot>=0) = lon_plot(lon_plot>=0)-360;
    lon_coast_plot(lon_coast_plot>=0) = lon_coast_plot(lon_coast_plot>=0)-360;
    lon_gl_plot(lon_gl_plot>=0) = lon_gl_plot(lon_gl_plot>=0)-360;
    Nx = size(df,1);
    Ny = size(df,2);
    dx = 500;
    dy = 500;
    xx = linspace(0,dx*(Nx-1),Nx);
    yy = linspace(0,dy*(Ny-1),Ny);
    [XX YY] = meshgrid(xx,yy);


    s_x = zeros(Nx,Ny); s_y = zeros(Nx,Ny); s = zeros(Nx,Ny);
    s_x(1:Nx-1,:) = (df(2:Nx,:)-df(1:Nx-1,:))./dx; %s_x(1:Nx-1,:) = (df(2:Nx,:)+df(1:Nx-1,:))/2; s_x(Nx,:) = 2*s_x(Nx-1,:) - s_x(Nx-2,:);
    s_y(:,1:Ny-1) = (df(:,2:Ny)-df(:,1:Ny-1))./dy; %s_y(:,1:Ny-1) = (df(:,2:Ny)+df(:,1:Ny-1))/2; s_y(:,Ny) = 2*s_y(:,Ny-1) - s_y(:,Ny-2);
    s = sqrt(s_x.^2 + s_y.^2);

    dtheta_c = (ff^2*Hpyc) ./ (g*gamma*s.^2);
    dtheta_c_plot = dtheta_c;
    dtheta_c_plot(dtheta_c_plot>8) = 10;
    dtheta_smooth = smoothdata2(dtheta_c_plot,'movmean',5,'omitnan');
    
    pcolor(ax1,lat,lon,dtheta_smooth); shading interp;
    hold on;
    set(ax1,'fontsize',28)
    set(ax1, 'XDir','reverse')
    cmocean('-curl','pivot',4);

    lat_min = -83.8; lat_max = -74.5;
    lon_min = -84; lon_max = -28;
    ylim([lon_min lon_max]);
    xlim([lat_min lat_max]);

    xlabel('Lat')
ylabel('Lon')

end
xticks(ax1,[-82,-79,-76]);
xticklabels(ax1,{'-82','-79','-76'});
yticks(ax1,[-80,-60,-40]);
yticklabels(ax1,{'-80','-60','-40'});

colorbar('box','on','linewidth',2,'location','northoutside','position',[0.11 0.95 0.675 0.02])
set(gca,'color',[0.95 0.95 0.95]);

% %%% create masks
[lati,loni] = meshgrid([lat_min:0.01:lat_max],[lon_min:0.1:lon_max]);
[Z_mask,Lat_mask,Lon_mask] = bedmachine_data('mask',lati,loni,'geo');
ocean_mask = Z_mask;
ocean_mask(ocean_mask~=0)=NaN;
%ocean_mask(ocean_mask==0)=10;

%%% plot mask
ax2 = axes('position',[0.3497    0.5443    0.2973    0.3832]);
hold on;
colormap(ax2,[0 0.3 1])
p = pcolor(ax2,Lat_mask,Lon_mask,ocean_mask);
shading interp;
alpha(p,0.3)
set(ax2, 'XDir','reverse')
%Turn off ax2 so you can see the first plot
ax2.Visible = 'off';
% set(ax2,'Color','None');
ax2.XTick = [];
ax2.YTick = [];

ylim([lon_min lon_max]);
xlim([lat_min lat_max]);

ax1.Box = 'on';
ax1.LineWidth = 2;
ax2.Box = 'on';
ax2.LineWidth = 2;

% Larsen
ax1 = subplot('Position',[0.7253    0.0932    0.1811    0.8331]);
hold on;
caxis([0 8]);
cm1 = cmocean('-curl','pivot',4);
colormap(ax1,cm1);


Larsen_shelves = {'LarsenA','LarsenB','LarsenC','LarsenD','LarsenD_1','LarsenE','LarsenF','LarsenG'};

for i = 1:length(Larsen_shelves)
    df = shelves.(Larsen_shelves{i});
    df(df>0) = NaN;
    lon = squeeze(coords.(Larsen_shelves{i})(1,:,:));
    lon(lon==Inf) = NaN;
    lat = squeeze(coords.(Larsen_shelves{i})(2,:,:));
    lat(lat==Inf) = NaN;

    Nx = size(df,1);
    Ny = size(df,2);
    dx = 500;
    dy = 500;
    xx = linspace(0,dx*(Nx-1),Nx);
    yy = linspace(0,dy*(Ny-1),Ny);
    [XX YY] = meshgrid(xx,yy);

    s_x = zeros(Nx,Ny); s_y = zeros(Nx,Ny); s = zeros(Nx,Ny);
    s_x(1:Nx-1,:) = (df(2:Nx,:)-df(1:Nx-1,:))./dx; %s_x(1:Nx-1,:) = (df(2:Nx,:)+df(1:Nx-1,:))/2; s_x(Nx,:) = 2*s_x(Nx-1,:) - s_x(Nx-2,:);
    s_y(:,1:Ny-1) = (df(:,2:Ny)-df(:,1:Ny-1))./dy; %s_y(:,1:Ny-1) = (df(:,2:Ny)+df(:,1:Ny-1))/2; s_y(:,Ny) = 2*s_y(:,Ny-1) - s_y(:,Ny-2);
    s = sqrt(s_x.^2 + s_y.^2);

    dtheta_c = (ff^2*Hpyc) ./ (g*gamma*s.^2);
    dtheta_c_plot = dtheta_c;
    dtheta_c_plot(dtheta_c_plot>8) = 10;
    dtheta_smooth = smoothdata2(dtheta_c_plot,'movmean',5,'omitnan');
    
    pcolor(ax1,lon,lat,dtheta_smooth); shading interp;
    hold on;
    cmocean('-curl','pivot',4);
    set(gca,'fontsize',28)

    lat_min = -75.5; lat_max = -64.2;
lon_min = -66; lon_max = -57;
ylim([lat_min lat_max])
xlim([lon_min lon_max])
xlabel('Lat')
ylabel('Lon')
    
end
xticks(ax1,[-65,-60]);
xticklabels(ax1,{'-65','-60'});
yticks(ax1,[-75,-70,-65]);
yticklabels(ax1,{'-75','-70','-65'});
set(gca,'color',[0.95 0.95 0.95]);

% %%% add masks
[lati,loni] = meshgrid([lat_min:0.01:lat_max],[lon_min:0.1:lon_max]);
[Z_mask,Lat_mask,Lon_mask] = bedmachine_data('mask',lati,loni,'geo');
ocean_mask = Z_mask;
ocean_mask(ocean_mask~=0)=NaN;


%%% plot mask
ax2 = axes('position',[0.7253    0.0932    0.1811    0.8331]);
hold on;
colormap(ax2,[0 0.3 1])
p = pcolor(ax2,Lon_mask,Lat_mask,ocean_mask); shading interp;
shading interp;
alpha(p,0.3)
%Turn off ax2 so you can see the first plot
ax2.Visible = 'off';
% set(ax2,'Color','None');
ax2.XTick = [];
ax2.YTick = [];

ylim([lat_min lat_max])
xlim([lon_min lon_max])

ax1.Box = 'on';
ax1.LineWidth = 2;
ax2.Box = 'on';
ax2.LineWidth = 2;



% Ross
ax1 = subplot('Position',[0.3483    0.0932    0.2986    0.3613]);
hold on;
caxis([0 8]);
cm1 = cmocean('-curl','pivot',4);
colormap(ax1,cm1);

Ross_shelves = {'Ross_West','Ross_East'};
fn = fieldnames(shelves);

for i = 1:length(Ross_shelves)
    df = shelves.(Ross_shelves{i});
    df(df>0) = NaN;
    lon = squeeze(coords.(Ross_shelves{i})(1,:,:));
    lon(lon==Inf) = NaN;
    lat = squeeze(coords.(Ross_shelves{i})(2,:,:));
    lat(lat==Inf) = NaN;
    %%% masks
    [lat_coast,lon_coast] = bedmachine_data('coast','geo');
    lon_coast_plot = lon_coast;
    [lat_gl,lon_gl] = bedmachine_data('gl','geo');
    lon_gl_plot = lon_gl;

    %%% adjusting longitude for plot
    lon_plot = lon;
    lon_plot(lon_plot>=0) = lon_plot(lon_plot>=0)-360;
    lon_coast_plot(lon_coast_plot>=0) = lon_coast_plot(lon_coast_plot>=0)-360;
    lon_gl_plot(lon_gl_plot>=0) = lon_gl_plot(lon_gl_plot>=0)-360;

    Nx = size(df,1);
    Ny = size(df,2);
    dx = 500;
    dy = 500;
    xx = linspace(0,dx*(Nx-1),Nx);
    yy = linspace(0,dy*(Ny-1),Ny);
    [XX YY] = meshgrid(xx,yy);

    s_x = zeros(Nx,Ny); s_y = zeros(Nx,Ny); s = zeros(Nx,Ny);
    s_x(1:Nx-1,:) = (df(2:Nx,:)-df(1:Nx-1,:))./dx; %s_x(1:Nx-1,:) = (df(2:Nx,:)+df(1:Nx-1,:))/2; s_x(Nx,:) = 2*s_x(Nx-1,:) - s_x(Nx-2,:);
    s_y(:,1:Ny-1) = (df(:,2:Ny)-df(:,1:Ny-1))./dy; %s_y(:,1:Ny-1) = (df(:,2:Ny)+df(:,1:Ny-1))/2; s_y(:,Ny) = 2*s_y(:,Ny-1) - s_y(:,Ny-2);
    s = sqrt(s_x.^2 + s_y.^2);

    dtheta_c = (ff^2*Hpyc) ./ (g*gamma*s.^2);
    dtheta_c_plot = dtheta_c;
    dtheta_c_plot(dtheta_c_plot>8) = 10;
    dtheta_smooth = smoothdata2(dtheta_c_plot,'movmean',5,'omitnan');
    
    pcolor(ax1,lon_plot,lat,dtheta_smooth); shading interp;
    hold on;
    cmocean('-curl','pivot',4);
    set(gca,'fontsize',28)
    set(gca, 'XDir','reverse')
    set(gca, 'YDir','reverse')


    lat_min = -86; lat_max = -77.4;
    lon_min = -202; lon_max = -146;
    xlim([lon_min lon_max]);
    ylim([lat_min lat_max]);

    xlabel('Lat');
    ylabel('Lon');


end
xticks(ax1,[-200, -180, -160]);
xticklabels(ax1,{'-200','-180','-160'});
yticks(ax1,[-86, -82, -78]);
yticklabels(ax1,{'-86','-82','-78'});
set(gca,'color',[0.95 0.95 0.95]);

% %%% add masks
[lati,loni] = meshgrid([lat_min:0.01:lat_max],[-179:0.1:lon_max]);
[Z_mask,Lat_mask,Lon_mask] = bedmachine_data('mask',lati,loni,'geo');
ocean_mask = Z_mask;
ocean_mask(ocean_mask~=0)=NaN;
ocean_mask(ocean_mask==0)=10;


%%% plot mask
ax2 = axes('position',[0.3483    0.0932    0.2986    0.3613]);
hold on;
colormap(ax2,[0 0.3 1])
p = pcolor(ax2,Lon_mask,Lat_mask,ocean_mask); shading interp;
shading interp;
alpha(p,0.3)
set(ax2, 'XDir','reverse')
set(ax2, 'YDir','reverse')
%Turn off ax2 so you can see the first plot
ax2.Visible = 'off';
% set(ax2,'Color','None');
ax2.XTick = [];
ax2.YTick = [];

ylim([lat_min lat_max]);
xlim([lon_min lon_max]);

% %%% add masks 2
[lati,loni] = meshgrid([lat_min:0.01:lat_max],[160:0.1:180]);
[Z_mask,Lat_mask,Lon_mask] = bedmachine_data('mask',lati,loni,'geo');
ocean_mask = Z_mask;
ocean_mask(ocean_mask~=0)=NaN;
ocean_mask(ocean_mask==0)=10;


%%% plot mask 2
ax3 = axes('position',[0.52801    0.0932    0.2986    0.3613]);
hold on;
colormap(ax3,[0 0.3 1])
p = pcolor(ax3,Lon_mask,Lat_mask,ocean_mask); shading interp;
alpha(p,0.3)
set(ax3, 'XDir','reverse')
set(ax3, 'YDir','reverse')
%Turn off ax2 so you can see the first plot
ax3.Visible = 'off';
% set(ax2,'Color','None');
ax3.XTick = [];
ax3.YTick = [];

ylim([lat_min lat_max]);
xlim([125 180]);

ax1.Box = 'on';
ax1.LineWidth = 2;
ax2.Box = 'off';
ax2.LineWidth = 2;
ax2.Box = 'off';
ax2.LineWidth = 2;




% % % % % % Shackleton
% % % % % subplot('Position',[0.7823    0.6352    0.1537    0.2808])
% % % % % hold on;
% % % % %
% % % % % df = shelves.Shackleton;
% % % % % df(df>0) = NaN;
% % % % % lon = squeeze(coords.Shackleton(1,:,:));
% % % % % lon(lon==Inf) = NaN;
% % % % % lon_plot = lon;
% % % % % lon_plot(lon_plot>=0) = lon_plot(lon_plot>=0)-360;
% % % % % lat = squeeze(coords.Shackleton(2,:,:));
% % % % % lat(lat==Inf) = NaN;
% % % % %
% % % % % %%% adjusting longitude for plot
% % % % % lon_plot = lon;
% % % % % lon_plot(lon_plot>=0) = lon_plot(lon_plot>=0)-360;
% % % % % lon_coast_plot(lon_coast_plot>=0) = lon_coast_plot(lon_coast_plot>=0)-360;
% % % % % lon_gl_plot(lon_gl_plot>=0) = lon_gl_plot(lon_gl_plot>=0)-360;
% % % % %
% % % % % Nx = size(df,1);
% % % % % Ny = size(df,2);
% % % % % dx = 500;
% % % % % dy = 500;
% % % % % xx = linspace(0,dx*(Nx-1),Nx);
% % % % % yy = linspace(0,dy*(Ny-1),Ny);
% % % % % [XX YY] = meshgrid(xx,yy);
% % % % %
% % % % % s_x = zeros(Nx,Ny); s_y = zeros(Nx,Ny); s = zeros(Nx,Ny);
% % % % % s_x(1:Nx-1,:) = (df(2:Nx,:)-df(1:Nx-1,:))./dx; %s_x(1:Nx-1,:) = (df(2:Nx,:)+df(1:Nx-1,:))/2; s_x(Nx,:) = 2*s_x(Nx-1,:) - s_x(Nx-2,:);
% % % % % s_y(:,1:Ny-1) = (df(:,2:Ny)-df(:,1:Ny-1))./dy; %s_y(:,1:Ny-1) = (df(:,2:Ny)+df(:,1:Ny-1))/2; s_y(:,Ny) = 2*s_y(:,Ny-1) - s_y(:,Ny-2);
% % % % % s = sqrt(s_x.^2 + s_y.^2);
% % % % %
% % % % % dtheta_c = (ff^2*Hpyc) ./ (g*gamma*s.^2);
% % % % % dtheta_c_plot = dtheta_c;
% % % % % dtheta_c_plot(dtheta_c_plot>8) = 10;
% % % % %
% % % % % cn = contourf(lon_plot,lat,smoothdata2(dtheta_c_plot,'movmean',5,'omitnan'),[0:1:8.5],'edgecolor','none');
% % % % % cmocean('-curl','pivot',4)
% % % % % set(gca,'fontsize',28)
% % % % % set(gca, 'XDir','reverse')
% % % % % set(gca, 'YDir','reverse')
% % % % % ax = gca;
% % % % % ax.Box = 'on';
% % % % % ax.LineWidth = 2;
% % % % % xticklabels([])
% % % % % yticklabels([])
% % % % % set(gca,'XTick',[])
% % % % % set(gca,'YTick',[])
% % % % %
% % % % %
% % % % % %%% add masks
% % % % % % ocean_coords = [10757:11659];
% % % % % % y_ocean = lat_coast(ocean_coords);
% % % % % % x_ocean = lon_coast_plot(ocean_coords);
% % % % % %
% % % % % % set(gca,'color',[0.95 0.95 0.95]);
% % % % % % area(x_ocean,y_ocean,'FaceColor',[0 0.3 1],'EdgeColor','none','FaceAlpha',0.15);
% % % % %
% % % % %
% % % % % lat_min = -64.95; lat_max = -66.8;
% % % % % lon_min = -265.2; lon_max = -259.25;
% % % % % ylim([-66.8 -64.95])
% % % % % xlim([-265.2 -259.25])
% % % % %
% % % % % % %%% add masks
% % % % %
% % % % % [lati,loni] = meshgrid([lat_min:0.01:lat_max],[lon_min:0.1:lon_max]);
% % % % % [Z_mask,Lat_mask,Lon_mask] = bedmachine_data('mask',lati,loni,'geo');
% % % % % ocean_mask = Z_mask;
% % % % % ocean_mask(ocean_mask~=0)=NaN;
% % % % % ocean_mask(ocean_mask==0)=10;
% % % % % pcolor(Lat_mask,Lon_mask,ocean_mask); shading interp;
% % % % %


% % Amery
% subplot('Position',[0.7823    0.3329    0.1537    0.2808])
% hold on;
% Amery_shelves = {'Amery'};
% 
% for i = 1:length(Amery_shelves)
%     df = shelves.(Amery_shelves{i});
%     df(df>0) = NaN;
%     lon = squeeze(coords.(Amery_shelves{i})(1,:,:));
%     lon(lon==Inf) = NaN;
%     lon_plot = lon;
%     lon_plot(lon_plot>=0) = lon_plot(lon_plot>=0)-360;
%     lat = squeeze(coords.(Amery_shelves{i})(2,:,:));
%     lat(lat==Inf) = NaN;
% 
%     %%% adjusting longitude for plot
%     lon_plot = lon;
%     lon_plot(lon_plot>=0) = lon_plot(lon_plot>=0)-360;
%     lon_coast_plot(lon_coast_plot>=0) = lon_coast_plot(lon_coast_plot>=0)-360;
%     lon_gl_plot(lon_gl_plot>=0) = lon_gl_plot(lon_gl_plot>=0)-360;
% 
%     Nx = size(df,1);
%     Ny = size(df,2);
%     dx = 500;
%     dy = 500;
%     xx = linspace(0,dx*(Nx-1),Nx);
%     yy = linspace(0,dy*(Ny-1),Ny);
%     [XX YY] = meshgrid(xx,yy);
% 
%     s_x = zeros(Nx,Ny); s_y = zeros(Nx,Ny); s = zeros(Nx,Ny);
%     s_x(1:Nx-1,:) = (df(2:Nx,:)-df(1:Nx-1,:))./dx; %s_x(1:Nx-1,:) = (df(2:Nx,:)+df(1:Nx-1,:))/2; s_x(Nx,:) = 2*s_x(Nx-1,:) - s_x(Nx-2,:);
%     s_y(:,1:Ny-1) = (df(:,2:Ny)-df(:,1:Ny-1))./dy; %s_y(:,1:Ny-1) = (df(:,2:Ny)+df(:,1:Ny-1))/2; s_y(:,Ny) = 2*s_y(:,Ny-1) - s_y(:,Ny-2);
%     s = sqrt(s_x.^2 + s_y.^2);
% 
%     dtheta_c = (ff^2*Hpyc) ./ (g*gamma*s.^2);
%     dtheta_c_plot = dtheta_c;
%     dtheta_c_plot(dtheta_c_plot>8) = 10;
% 
%     cn = contourf(lat,lon,smoothdata2(dtheta_c_plot,'movmean',5,'omitnan'),[0:1:8.5],'edgecolor','none');
% end
% 
% set(gca,'color',[0.95 0.95 0.95]);
% %%% add masks
% lat_min = -73.4; lat_max = -68.4;
% lon_min = 66.75; lon_max = 74.3;
% ylim([66.7 74.3])
% xlim([-73.4 -68.4])
% 
% % %%% add masks
% freezeColors;
% [lati,loni] = meshgrid([lat_min:0.01:lat_max],[lon_min:0.1:lon_max]);
% [Z_mask,Lat_mask,Lon_mask] = bedmachine_data('mask',lati,loni,'geo');
% ocean_mask = Z_mask;
% ocean_mask(ocean_mask~=0)=NaN;
% ocean_mask(ocean_mask==0)=100;
% p = pcolor(Lat_mask,Lon_mask,ocean_mask);shading flat;
% 
% 
% caxis([0 8]);
% cm = cmocean('-curl','pivot',4);
% cm = [cm; 0 0.3 1];
% ax = gca;
% ax.Box = 'on';
% ax.LineWidth = 2;
% set(gca,'fontsize',28)
