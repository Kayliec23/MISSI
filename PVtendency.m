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

%%% define b
% calculate rho
rho = NaN(Nx,Ny,Nr,Nt);
for n = 1:Nt
  SS = squeeze(S(1,:,:,n));
  TT = squeeze(T(1,:,:,n));
  rho(Nx,:,:,n) = densjmd95(SS,TT,0);
end

b = -g/rho0*rho;

%%% define bdot
% calculate rho_dot from TOTSTEND and TOTTTEND
rho_dot = NaN(Nx,Ny,Nr,Nt);
for n = 1:Nt
  SS = squeeze(TOTSEND(1,:,:,n));
  TT = squeeze(TOTTEND(1,:,:,n));
  rho_dot(Nx,:,:,n) = densjmd95(SS,TT,0);
end
bdot = -g/rho0*rho_dot;

% calculate rho_diff
rho_diff = NaN(Nx,Ny,Nr,Nt);
for n = 1:Nt
  SS = squeeze(DFxE_SLT(1,:,:,n));
  TT = squeeze(DFxE_TH(1,:,:,n));
  rho_diff(Nx,:,:,n) = densjmd95(SS,TT,0);
end
bdot_diff = -g/rho0*rho_diff;

% calculate rho_adv
rho_adv = NaN(Nx,Ny,Nr,Nt);
for n = 1:Nt
  SS = squeeze(ADVx_SLT(1,:,:,n));
  TT = squeeze(ADVx_TH(1,:,:,n));
  rho_adv(Nx,:,:,n) = densjmd95(SS,TT,0);
end
bdot_adv = -g/rho0*rho_adv;

%%% define u
u = U;

%%% define udot
udot = TOTUTEND:
udot_diss = Um_diss;
udot_adv = Um_Advec;
udot_cori = Um_cori;
udot_p = UmdPhiX;

%%% term1 
A = mean(f0*bdot,2,'omitnan');    % mean in y
term1 = (diff(A,1,3))./dz;    % diff in z
term1= interp_center_y(term1);    % interpolate and center in z direction

%%% term2 
A = mean(diff(udot,1,2),3); % diff in y --> mean in z
A = interp_center_y(A);    % interpolate and center in z direction
B = mean(diff(b,1,3),2); % diff in z --> mean in y
B = interp_center_z(B);    % interpolate and center in z direction
term2 = -(A.*B)./dy./dz;

%%% term3
A = mean(diff(u,1,2),3); % diff in y --> mean in z
A = interp_center_y(A);    % interpolate and center in z direction
B = mean(diff(bdot,1,3),2);    % diff in z --> mean in y
B = interp_center_z(B);    % interpolate and center in z direction
term3 = -(A.*B)./dy./dz;

%%% term4 
A = mean(diff(udot,1,3),2); % diff in z --> mean in y
A = interp_center_z(A);    % interpolate and center in z direction
B = mean(diff(b,1,2),3);   % diff in y --> mean in z
B = interp_center_y(B);    % interpolate and center in z direction
term4 = (A.*B)./dy.dz;

%%% term5 
A = mean(diff(u,1,3),2);    % diff in z --> mean in y
A = interp_center_z(A);    % interpolate and center in z direction
B = mean(diff(bdot,1,2),3);    % diff in y --> mean in z
B = interp_center_y(B);    % interpolate and center in z direction
term5 = (A.*B)./dy./dz;

%%% qdot
qdot = term1 + term2 + term3 + term4 + term5; 


function BB = interp_center_z(AA)
  AA(:,:,Nr,:) = 2*AA(:,:,Nr-1,:) - AA(:,:,Nr-2,:);    % interpolate in z for boundary point
  AA = ( AA(:,1:Nr-1,:,:) + AA(:,2:Nr,:,:) ) / 2;    % recenter
  BB = AA;
end

function BB = interp_center_y(AA)
  AA(:,Ny,:,:) = 2*AA(:,Ny-1,:,:)-AA(:,Ny-2,:,:);    % interpolate in y for boundary point
  AA = ( AA(:,1:Ny-1,:,:) + AA(:,2:Ny,:,:) ) / 2;    % recenter
  BB = AA;
end








