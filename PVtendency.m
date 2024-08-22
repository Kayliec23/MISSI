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

% addpath /Users/Kayliec23/Desktop/analysis

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% calculate b and bdot terms %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% define b and its components
% find alpha and beta
SS = iceNaN(S,Nx,Ny,Nr,icetopo,zz); SS = squeeze(nanmean(SS(1,:,:,:),4));
TT = iceNaN(T,Nx,Ny,Nr,icetopo,zz); TT = squeeze(nanmean(TT(1,:,:,:),4));
dPT = 1e-3;
dSP = 1e-3;
alpha = - (densjmd95(SS,TT+0.5*dPT,0) - densjmd95(SS,TT-0.5*dPT,0))./densjmd95(SS,TT,0)/dPT;
beta = (densjmd95(SS+0.5*dSP,TT,0) - densjmd95(SS-0.5*dSP,TT,0))./densjmd95(SS,TT,0)/dSP;


%%% calculate b
for n = 1:Nt
    b_(1,:,:,n) = g.*(alpha.*squeeze(T(1,:,:,n)) - beta.*squeeze(S(1,:,:,n))); % - ...
    %g.*(alpha.*squeeze(T(1,:,:,1)) - beta.*squeeze(S(1,:,:,1)));
end
%b = (b(:,:,:,1:Nt-1)+b(:,:,:,2:Nt))./2;
b_ = iceNaN(b_,Nx,Ny,Nr,icetopo,zz);
b = (b_(:,:,:,1:Nt-1)+b_(:,:,:,2:Nt))/2;
bt = diff(b,1,4)./t1hour;
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
%     case 'adv'
%         udot_comp = Um_Advec;
%         bdot_comp = bdot_adv;
%     case 'diff'
%         udot_comp = udot_diff;
%         bdot_comp = bdot_diff;
%     case 'visc'
%         udot_comp = Udot_visc;
%         bdot_comp = bdot_visc;
end


%%% term1
term1 = f0*(bdot_comp(:,:,[1:Nr-1],:) - bdot_comp(:,:,[2:Nr],:))./dz;    % diff in z
term1 = (term1(:,1:Ny-1,:,:) + term1(:,2:Ny,:,:)) / 2;    % average locally in y
term1 = interp_yz(term1);    % interpolate in z and y

%%% term2
A = diff(udot_comp,1,2); % diff in y
A = (A(:,:,1:Nr-1,:) + A(:,:,2:Nr,:)) / 2; % averaging locally in z
A = interp_yz(A);    % interpolate in z and y
B = (b(:,:,[1:Nr-1],:) - b(:,:,[2:Nr],:));   % diff in z
B = (B(:,1:Ny-1,:,:) + B(:,2:Ny,:,:)) / 2;    % average locally in y
B = interp_yz(B);    % interpolate in z and y
term2 = -(A.*B)./dy./dz;

%%% term3
A = diff(u,1,2); % diff in y
A = (A(:,:,1:Nr-1,:) + A(:,:,2:Nr,:)) / 2; % averaging locally in z
A = interp_yz(A);    % interpolate in z and y
B = (bdot_comp(:,:,[1:Nr-1],:) - bdot_comp(:,:,[2:Nr],:));    % diff in z
B = (B(:,1:Ny-1,:,:) + B(:,2:Ny,:,:)) / 2;    % average locally in y
B = interp_yz(B);    % interpolate in z and y
term3 = -(A.*B)./dy./dz;

%%% term4
A = (udot_comp(:,:,[1:Nr-1],:) - udot_comp(:,:,[2:Nr],:)); % diff in z
A = (A(:,1:Ny-1,:,:) + A(:,2:Ny,:,:)) / 2;    % average locally in y
A = interp_yz(A);    % interpolate in z and y
B = diff(b,1,2); % diff in y
B = (B(:,:,1:Nr-1,:) + B(:,:,2:Nr,:)) / 2;    % average locally in z
B = interp_yz(B);    % interpolate in z and y
term4 = (A.*B)./dy./dz;

%%% term5
A = (u(:,:,[1:Nr-1],:) - u(:,:,[2:Nr],:));    % diff in z
A = (A(:,1:Ny-1,:,:) + A(:,2:Ny,:,:)) / 2;    % average locally in y
A = interp_yz(A);    % interpolate in z and y
B = diff(bdot_comp,1,2);    % diff in y
B = (B(:,:,1:Nr-1,:) + B(:,:,2:Nr,:)) / 2;    % average locally in z
B = interp_yz(B);    % interpolate in z and y
term5 = (A.*B)./dy./dz;

%%% qdot
qdot = term1 + term2 + term3 + term4 + term5;

%%% Calculate Q and take time derivative
Q = calculateQ(S,T,U,V,W,Nx,Ny,Nr,Nt,dy,dz,g,rho0,f0);
Q = iceNaN(Q,Nx,Ny,Nr,icetopo,zz);
Qt = diff(Q,1,4)./t1hour;

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


qdot_mean = plume_mean(qdot,Nd,Ns,icetopo,zz,dz);
term1_mean = plume_mean(term1,Nd,Ns,icetopo,zz,dz);
term2_mean = plume_mean(term2,Nd,Ns,icetopo,zz,dz);
term3_mean = plume_mean(term3,Nd,Ns,icetopo,zz,dz);
qvert_dot_mean = term3_mean+term2_mean;
term4_mean = plume_mean(term4,Nd,Ns,icetopo,zz,dz);
term5_mean = plume_mean(term5,Nd,Ns,icetopo,zz,dz);
qbc_dot_mean = term4_mean+term5_mean;
bdot_mean = plume_mean(bdot,Nd,Ns,icetopo,zz,dz);
Qt_mean = plume_mean(Qt,Nd,Ns,icetopo,zz,dz);
Q_mean = plume_mean(Q,Nd,Ns,icetopo,zz,dz);

% %%
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
% KEt(:,:,:,1:Nt-1) = (KE(:,:,:,2:Nt) - KE(:,:,:,1:Nt-1))/t1hour;
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Enstrophy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Total, mean and perturbation enstrophies
E = plume_mean(0.5*Q.^2,Nd,Ns,icetopo,zz,dz);
E_mean = 0.5*Q_mean.^2; 
E_pert = E - E_mean;

%%% Mean and perturbation enstrophy tendencies
Et_mean = diff(E_mean,1,2)./t1hour;
Et_pert = diff(E_pert,1,2)./t1hour;

%%% Total, mean and perturbation terms in the enstrophy budget

Eqdot = plume_mean(0.5*(Q(:,:,:,1:end-1)+Q(:,:,:,2:end)).*qdot,Nd,Ns,icetopo,zz,dz);
Eqdot_mean = 0.5*(Q_mean(:,1:end-1)+Q_mean(:,2:end)).*qdot_mean;
Eqdot_pert = Eqdot - Eqdot_mean;

Eterm1 = plume_mean(0.5*(Q(:,:,:,1:end-1)+Q(:,:,:,2:end)).*term1,Nd,Ns,icetopo,zz,dz);
Eterm1_mean = 0.5*(Q_mean(:,1:end-1)+Q_mean(:,2:end)).*term1_mean;
Eterm1_pert = Eterm1 - Eterm1_mean;

Eterm2 = plume_mean(0.5*(Q(:,:,:,1:end-1)+Q(:,:,:,2:end)).*term2,Nd,Ns,icetopo,zz,dz);
Eterm2_mean = 0.5*(Q_mean(:,1:end-1)+Q_mean(:,2:end)).*term2_mean;
Eterm2_pert = Eterm2 - Eterm2_mean;

Eterm3 = plume_mean(0.5*(Q(:,:,:,1:end-1)+Q(:,:,:,2:end)).*term3,Nd,Ns,icetopo,zz,dz);
Eterm3_mean = 0.5*(Q_mean(:,1:end-1)+Q_mean(:,2:end)).*term3_mean;
Eterm3_pert = Eterm3 - Eterm3_mean;

Eterm4 = plume_mean(0.5*(Q(:,:,:,1:end-1)+Q(:,:,:,2:end)).*term4,Nd,Ns,icetopo,zz,dz);
Eterm4_mean = 0.5*(Q_mean(:,1:end-1)+Q_mean(:,2:end)).*term4_mean;
Eterm4_pert = Eterm4 - Eterm4_mean;

Eterm5 = plume_mean(0.5*(Q(:,:,:,1:end-1)+Q(:,:,:,2:end)).*term5,Nd,Ns,icetopo,zz,dz);
Eterm5_mean = 0.5*(Q_mean(:,1:end-1)+Q_mean(:,2:end)).*term5_mean;
Eterm5_pert = Eterm5 - Eterm5_mean;
% 
% %%
% bdotdz(:,:,2:Nr,:) = (bdot_comp(:,:,[1:Nr-1],:) - bdot_comp(:,:,[2:Nr],:))./dz;
% bdotdz(:,:,1,:) = 2*bdotdz(:,:,2,:) - bdotdz(:,:,3,:);
% wbdotdz = 0.5*(W(:,:,:,1:Nt-1)+W(:,:,:,2:Nt)).*bdotdz;
% wbdotdz_mean = plume_mean(wbdotdz,Nd,Ns,icetopo,zz,dz);
% bdotdy(:,1:Ny-1,:,:) = diff(bdot_comp,1,2)./dy;
% bdotdy(:,Ny,:,:) = 2*bdotdy(:,Ny-1,:,:) - bdotdy(:,Ny-2,:,:);
% vbdotdy = 0.5*(V(:,:,:,1:Nt-1)+V(:,:,:,2:Nt)).*bdotdy;
% vbdotdy_mean = plume_mean(vbdotdy,Nd,Ns,icetopo,zz,dz);
% %%
% dbdz(:,:,2:Nr,:) = (b(:,:,[1:Nr-1],:) - b(:,:,[2:Nr],:))./dz;
% dbdz(:,:,1,:) = 2*dbdz(:,:,2,:) - dbdz(:,:,3,:);
% wdbdz = 0.5*(W(:,:,:,1:Nt-1)+W(:,:,:,2:Nt)).*dbdz;
% wdbdz_mean = plume_mean(wdbdz,Nd,Ns,icetopo,zz,dz);
% dbdy(:,1:Ny-1,:,:) = diff(b,1,2)./dy;
% dbdy(:,Ny,:,:) = 2*dbdy(:,Ny-1,:,:) - dbdy(:,Ny-2,:,:);
% vdbdy = 0.5*(V(:,:,:,1:Nt-1)+V(:,:,:,2:Nt)).*dbdy;
% vdbdy_mean = plume_mean(vdbdy,Nd,Ns,icetopo,zz,dz);
%%
%%% Quick plot of Q snapshot
figure(1);pcolor(YY',ZZ'-(YY'-4e4)*s,squeeze(Q(1,:,:,120)));shading flat;colorbar;caxis([-1 1]*3e-8);colormap redblue; set(gca,'YLim',[-20 0]);

figure(2);pcolor(NT,D,sqrt(E_pert(:,2:end))'-abs(Q_mean(:,2:end))');colorbar;shading flat;caxis([-1 1]*1e-7);colormap redblue;

figure(3);pcolor(NT,D,E_mean(:,2:end)');colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
figure(4);pcolor(NT,D,E_pert(:,2:end)');colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;

figure(5);pcolor(NT,D,Et_mean');colorbar;shading flat;caxis([-1 1]*1e-21);colormap redblue;
figure(6);pcolor(NT,D,Et_pert');colorbar;shading flat;caxis([-1 1]*1e-21);colormap redblue;
figure(7);pcolor(NT,D,Et_pert'-Et_mean');colorbar;shading flat;caxis([-1 1]*1e-21);colormap redblue;

figure(8);pcolor(NT,D,cumsum(Eqdot_mean',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
figure(9);pcolor(NT,D,cumsum(Eqdot_pert',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;

figure(10);pcolor(NT,D,cumsum(Eterm1_mean',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
figure(11);pcolor(NT,D,cumsum(Eterm1_pert',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;

figure(12);pcolor(NT,D,cumsum(Eterm2_mean',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
figure(13);pcolor(NT,D,cumsum(Eterm2_pert',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;

figure(14);pcolor(NT,D,cumsum(Eterm3_mean',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
figure(15);pcolor(NT,D,cumsum(Eterm3_pert',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;

figure(16);pcolor(NT,D,cumsum(Eterm4_mean',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
figure(17);pcolor(NT,D,cumsum(Eterm4_pert',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;

figure(18);pcolor(NT,D,cumsum(Eterm5_mean',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;
figure(19);pcolor(NT,D,cumsum(Eterm5_pert',1)*t1hour);colorbar;shading flat;caxis([-1 1]*1e-15);colormap redblue;

%%% Variance ratio budget terms
ENqdot = Eqdot_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eqdot_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2; 
ENterm1 = Eterm1_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eterm1_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2; 
ENterm2 = Eterm2_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eterm2_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2; 
ENterm3 = Eterm3_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eterm3_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2; 
ENterm4 = Eterm4_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eterm4_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2; 
ENterm5 = Eterm5_pert./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))) - (0.5*(E_pert(:,1:end-1)+E_pert(:,2:end))).*Eterm5_mean./(0.5*(E_mean(:,1:end-1)+E_mean(:,2:end))).^2; 
figure(20);pcolor(NT,D,(E_pert(:,2:end)'./E_mean(:,2:end)'));colorbar;shading flat;caxis([-10 10]);colormap redblue;
figure(21);pcolor(NT,D,cumsum(ENqdot',1)*t1hour);colorbar;shading flat;caxis([-10 10]);colormap redblue;
figure(22);pcolor(NT,D,cumsum(ENterm1',1)*t1hour);colorbar;shading flat;caxis([-10 10]);colormap redblue;
figure(23);pcolor(NT,D,cumsum(ENterm2',1)*t1hour);colorbar;shading flat;caxis([-10 10]);colormap redblue;
figure(24);pcolor(NT,D,cumsum(ENterm3',1)*t1hour);colorbar;shading flat;caxis([-10 10]);colormap redblue;
figure(25);pcolor(NT,D,cumsum(ENterm4',1)*t1hour);colorbar;shading flat;caxis([-10 10]);colormap redblue;
figure(26);pcolor(NT,D,cumsum(ENterm5',1)*t1hour);colorbar;shading flat;caxis([-10 10]);colormap redblue;

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
