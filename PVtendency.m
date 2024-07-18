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
[YY,ZZ] = meshgrid(yy,zz);

addpath /Users/Kayliec23/Desktop/analysis

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% calculate b and bdot terms %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% define b and its components
% find alpha and beta
SS = squeeze(nanmean(S(1,:,:,:),4)); 
TT = squeeze(nanmean(T(1,:,:,:),4));
dPT = 1e-3;
dSP = 1e-3;
alpha = - (densjmd95(SS,TT+0.5*dPT,0) - densjmd95(SS,TT-0.5*dPT,0))./densjmd95(SS,TT,0)/dPT; 
beta = (densjmd95(SS+0.5*dSP,TT,0) - densjmd95(SS-0.5*dSP,TT,0))./densjmd95(SS,TT,0)/dSP;


%%% calculate b
for n = 1:Nt
    b(1,:,:,n) = g.*(alpha.*squeeze(T(1,:,:,n)) - beta.*squeeze(S(1,:,:,n))); %- ...
        %g.*(alpha.*squeeze(T(1,:,:,1)) - beta.*squeeze(S(1,:,:,1)));
end
b = (b(:,:,:,1:Nt-1)+b(:,:,:,2:Nt))./2;
b = iceNaN(b,Nx,Ny,Nr,icetopo,zz);
% calculate bdot_tot
[bdot bdotdot] = buoyancy(TOTSTEND./t1day,TOTTTEND./t1day,alpha,beta,icetopo,zz);
bdot = iceNaN(bdot,Nx,Ny,Nr,icetopo,zz);
% % calculate rho_diff
% [bDIFFx bDIFFx_dot] = buoyancy(DFxE_SLT,DFxE_TH,alpha,beta,icetopo,zz);
% [bDIFFy bDIFFy_dot] = buoyancy(DFyE_SLT,DFyE_TH,alpha,beta,icetopo,zz);
% [bDIFFrR bDIFFrE_dot] = buoyancy(DFrE_SLT,DFrE_TH,alpha,beta,icetopo,zz);
% [bDIFFrI bDIFFrI_dot] = buoyancy(DFrI_SLT,DFrI_TH,alpha,beta,icetopo,zz);
% 
% % calculate rho_adv
% [bADVx bADVx_dot] = buoyancy(ADVx_SLT,ADVx_TH,alpha,beta,icetopo,zz);
% [bADVy bADVy_dot] = buoyancy(ADVy_SLT,ADVy_TH,alpha,beta,icetopo,zz);
% [bADVr bADVr_dot] = buoyancy(ADVr_SLT,ADVr_TH,alpha,beta,icetopo,zz);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% calculate qdot mean terms %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% define u and its components
u = (U(:,:,:,1:Nt-1)+U(:,:,:,2:Nt))/2;
% U_adv = ADVx_Um + ADVy_Um + ADVrE_Um;
% U_visc = VISCx_Um + VISCy_Um + VISrE_Um + VISrI_Um;
% Udot_visc = (U_visc(:,:,:,2:Nt)-U_visc(:,:,:,1:Nt-1)) ./ t1hour;

%%% define udot and its components
udot = TOTUTEND./t1day;
udot = iceNaN(udot,Nx,Ny,Nr,icetopo,zz);
% Udot_Diss = (Um_Diss(:,:,:,2:Nt)-Um_Diss(:,:,:,1:Nt-1)) ./ t1hour;
% Udot_ImplD = (Um_Impl(:,:,:,2:Nt)-Um_ImplD(:,:,:,1:Nt-1)) ./ t1hour;
% U_diff = U_Diss + U_ImplD;

component = 'tot';

switch component
    case 'tot'
        udot_comp = udot;
        bdot_comp = bdot;
    case 'adv'
        udot_comp = Um_Advec;
        bdot_comp = bdot_adv;
    case 'diff'
        udot_comp = udot_diff;
        bdot_comp = bdot_diff;
    case 'visc'
        udot_comp = Udot_visc;
        bdot_comp = bdot_visc;
end


%%% term1 
term1 = (diff(f0*bdot_comp,1,3))./dz;    % diff in z
term1 = (term1(:,1:Ny-1,:,:) + term1(:,2:Ny,:,:)) / 2;    % average locally in y
term1 = interp_yz(term1);    % interpolate in z and y

%%% term2 
A = diff(udot_comp,1,2); % diff in y 
A = (A(:,:,1:Nr-1,:) + A(:,:,2:Nr,:)) / 2; % averaging locally in z
A = interp_yz(A);    % interpolate in z and y
B = diff(b,1,3);   % diff in z
B = (B(:,1:Ny-1,:,:) + B(:,2:Ny,:,:)) / 2;    % average locally in y
B = interp_yz(B);    % interpolate in z and y
term2 = -(A.*B)./dy./dz;

%%% term3
A = diff(u,1,2); % diff in y 
A = (A(:,:,1:Nr-1,:) + A(:,:,2:Nr,:)) / 2; % averaging locally in z
A = interp_yz(A);    % interpolate in z and y
B = diff(bdot_comp,1,3);    % diff in z 
B = (B(:,1:Ny-1,:,:) + B(:,2:Ny,:,:)) / 2;    % average locally in y
B = interp_yz(B);    % interpolate in z and y
term3 = -(A.*B)./dy./dz;

%%% term4 
A = diff(udot_comp,1,3); % diff in z 
A = (A(:,1:Ny-1,:,:) + A(:,2:Ny,:,:)) / 2;    % average locally in y
A = interp_yz(A);    % interpolate in z and y
B = diff(b,1,2); % diff in y 
B = (B(:,:,1:Nr-1,:) + B(:,:,2:Nr,:)) / 2;    % average locally in z
B = interp_yz(B);    % interpolate in z and y
term4 = (A.*B)./dy./dz;

%%% term5 
A = diff(u,1,3);    % diff in z
A = (A(:,1:Ny-1,:,:) + A(:,2:Ny,:,:)) / 2;    % average locally in y
A = interp_yz(A);    % interpolate in z and y
B = diff(bdot_comp,1,2);    % diff in y 
B = (B(:,:,1:Nr-1,:) + B(:,:,2:Nr,:)) / 2;    % average locally in z
B = interp_yz(B);    % interpolate in z and y
term5 = (A.*B)./dy./dz;

%%% qdot
qdot = term1 + term2 + term3 + term4 + term5; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate plume means %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% calculate d, distance from slope
Nd = 50;
d = zeros(Ny,Nd);
Ns = 400;
for j = 1:Ny
    k = find(zz==floor(icetopo(j))+dz/2);
    if icetopo(j)-zz(k)>0
        d(j,:) = icetopo(j)-zz(k:k+(Nd-1));
    elseif icetopo(j)-zz(k)<0
        d(j,:) = icetopo(j)-zz(k+1:k+Nd);
    end
end
dd = mean(d(1:Ns,:),1); dd = dd-dd(1);
[D,NT] = meshgrid(dd,[1:Nt-1]);
%%
qdot_mean = plume_mean(qdot,Nd,Ns,icetopo,zz,dz);
term1_mean = plume_mean(term1,Nd,Ns,icetopo,zz,dz);
term2_mean = plume_mean(term2,Nd,Ns,icetopo,zz,dz);
term3_mean = plume_mean(term3,Nd,Ns,icetopo,zz,dz);
qvert_dot_mean = term3_mean+term2_mean;
term4_mean = plume_mean(term4,Nd,Ns,icetopo,zz,dz);
term5_mean = plume_mean(term5,Nd,Ns,icetopo,zz,dz);
qbc_dot_mean = term4_mean+term5_mean;
bdot_mean = plume_mean(bdot,Nd,Ns,icetopo,zz,dz);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scrsz = get(0,'ScreenSize');
fontsize = 26;
framepos = [0 scrsz(4)/2 scrsz(3) scrsz(4)];
plotloc = [0.17 0.18 0.62 0.7];
handle = figure(17);
set(handle,'Position',framepos);
clf;
set(gcf,'Color','w');

h1 = subplot(2,3,1);
pcolor(NT,D,q_t_mean'); colorbar; shading interp;
mymap = colormap('redblue');
colormap(h1,mymap)
caxis([-5 5]*1e-13)
hold on;
ylim([0 25])
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
ylim([0 25])
set(gca,'fontsize',fontsize)
ylabel('d (m)')
xticklabels([])
title('$\dot{q}$','interpreter','latex')

h3 = subplot(2,3,3);
pcolor(NT,D,(qdot_mean - q_t_mean)'); colorbar; shading interp;
mymap = colormap('redblue');
colormap(h3,mymap)
caxis([-5 5]*1e-13)
hold on;
ylim([0 25])
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
ylim([0 25])
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
ylim([0 25])
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
ylim([0 25])
set(gca,'fontsize',fontsize)
ylabel('d (m)')
xticks([24,120,240,360]);
xticklabels({'1','5','10','15'})
xlabel('days')
title('$\dot{q_{bc}} - \dot{q}$','interpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b bdot] = buoyancy(S,T,alpha,beta,icetopo,zz)
    t1min = 60;
    t1hour = 60*t1min;
    rho0 = 1027;
    g = 9.81;
    Nx = size(S,1);
    Ny = size(S,2);
    Nr = size(S,3);
    Nt = size(S,4);
    for n = 1:Nt
        b(1,:,:,n) = g.*(alpha.*squeeze(T(1,:,:,n))-beta.*squeeze(S(1,:,:,n))); %- ...
        %g.*(alpha.*squeeze(T(1,:,:,1)) - beta.*squeeze(S(1,:,:,1)));
    end
    b = iceNaN(b,Nx,Ny,Nr,icetopo,zz);
    bdot = (b(:,:,:,2:Nt)-b(:,:,:,1:Nt-1)) ./ t1hour;
end

function BB = interp_yz(AA)
    Nx = size(AA,1);
    Ny = size(AA,2);
    Nr = size(AA,3);
    Nt = size(AA,4);
    BB = NaN(Nx,Ny+1,Nr+1,Nt);
    BB(:,2:Ny+1,2:Nr+1,:) = AA;
    BB(:,:,1,:) = 2*BB(:,:,2,:) - BB(:,:,3,:);   % interpolate in z for boundary point
    BB(:,1,:,:) = 2*BB(:,2,:,:) - BB(:,3,:,:);    % interpolate in y for boundary point
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

