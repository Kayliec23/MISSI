%%% Calculating PV tendency using MITgcm model diagnostics 

%%% need u,b and udot, bdot (dot = tendency).
%%% Um_Diss - U momentum tendency from Dissipation (Explicit part)
%%% Um_Impl - U momentum tendency from Dissipation (Implicit part)
%%% Um_Advec - U momentum tendency from Advection terms
%%% Um_Cori - U momentum tendency from Coriolis term (this term is part of Um_Advec in current model setting) 
%%% Um_dPhiX - U momentum tendency from Pressure/Potential grad (in 2D this term is zero)
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
dy = diff(yy,1);
dx = diff(xx,1);
dz = diff(zz,1);
Omega = 2*pi*366/365/86400;
lat0 = -80; %%% Latitude at southern boundary
f0 = 2*Omega*sind(lat0); %%% Coriolis parameter
m1km = 1000;
g = 9.81;
rho0 = 1027;
t1min = 50;
t1hour = 60*t1min;
t1day = 24*t1hour;
[YY,ZZ] = meshgrid(yy,zz);

addpath /Users/Kayliec23/Desktop/analysis

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% define b and u and their components %%%%%%%%%%%%%%%%%%%%%%%%
 
%%% define b and its components
% find alpha and beta
S = squeeze(nanmean(SS(1,:,:,:),4));
T = squeeze(nanmean(TT(1,:,:,1),4));
g = 9.81;
dPT = 1e-3;
dSP = 1e-3;
alpha = - (densjmd95(S,T+0.5*dPT,0) - densjmd95(S,T-0.5*dPT,0))./densjmd95(S,T,0)/dPT;
beta = (densjmd95(S+0.5*dSP,T,0) - densjmd95(S-0.5*dSP,T,0))./densjmd95(S,T,0)/dSP;
b = g.*(alpha.*T-beta.*S);
bdot = (b(:,:,:,2:Nt)-b(:,:,:,1:Nt-1)) ./ t1hour; 

% % calculate rho_diff
bDiff = g.*(alpha.*THdiffX-beta.*SdiffX);
bDiff_dot = (bDiff(:,:,:,2:Nt)-bDiff(:,:,:,1:Nt-1)) ./ t1hour; 
% 
% % calculate rho_adv
bAdv = g.*(alpha.*THdiffX-beta.*SdiffX);
bAdv_dot = (bAdv(:,:,:,2:Nt)-bAdv(:,:,:,1:Nt-1)) ./ t1hour; 

% % calculate rho_adv
btend = g.*(alpha.*TOTTTEND-beta.*TOTSTEND);

%%% define u and its components
u = U;

%%% define udot
udot = TOTUTEND/86400;


%%%%%%%%%%%%%%%%%%%%%%%%%% calculate qdot mean terms %%%%%%%%%%%%%%%%%%%%%%%%
 
%%% term1 
A = nanmean(f0*bdot,2,'omitnan');    % mean in y
term1 = (diff(A,1,3))./dz;    % diff in z
term1= interp_center_y(term1,Ny);    % interpolate and center in z direction --> mean in z

%%% term2 
A = diff(udot,1,2); % diff in y 
A = nanmean(interp_center_y(A,Ny),3);    % interpolate and center in y direction --> mean in z
B = diff(b,1,3); % diff in z 
B = nanmean(interp_center_z(B,Nr),2);    % interpolate and center in z direction --> mean in y
term2 = -(A.*B)./dy./dz;

%%% term3
A = diff(u,1,2); % diff in y 
A = nanmean(interp_center_y(A,Ny),3);    % interpolate and center in z direction --> mean in z
B = diff(bdot,1,3);    % diff in z 
B = nanmean(interp_center_z(B,Nr),2);    % interpolate and center in z direction --> mean in y
term3 = -(A.*B)./dy./dz;

%%% term4 
A = diff(udot,1,3); % diff in z 
A = nanmean(interp_center_z(A,Nr),2);    % interpolate and center in z direction --> mean in y
B = diff(b,1,2);   % diff in y 
B = nanmean(interp_center_y(B,Ny),3);    % interpolate and center in z direction --> mean in z
term4 = (A.*B)./dy.dz;

%%% term5 
A = diff(u,1,3);    % diff in z
A = nanmean(interp_center_z(A,Nr),2);    % interpolate and center in z direction --> mean in y
B = diff(bdot,1,2);    % diff in y 
B = nanmean(interp_center_y(B,Ny),3);    % interpolate and center in z direction --> mean in z
term5 = (A.*B)./dy./dz;

%%% qdot
qdot = term1 + term2 + term3 + term4 + term5; 

%%%%%%%%%%%%%%%%%%%%%%%%%% calculate qdot perturbations terms %%%%%%%%%%%%%%%%%%%%%%%%
qmean = 









function BB = interp_center_z(AA,Nr)
  AA(:,:,Nr,:) = 2*AA(:,:,Nr-1,:) - AA(:,:,Nr-2,:);    % interpolate in z for boundary point
  AA = ( AA(:,1:Nr-1,:,:) + AA(:,2:Nr,:,:) ) / 2;    % recenter
  BB = AA;
end

function BB = interp_center_y(AA,Ny)
  AA(:,Ny,:,:) = 2*AA(:,Ny-1,:,:)-AA(:,Ny-2,:,:);    % interpolate in y for boundary point
  AA = ( AA(:,1:Ny-1,:,:) + AA(:,2:Ny,:,:) ) / 2;    % recenter
  BB = AA;
end







