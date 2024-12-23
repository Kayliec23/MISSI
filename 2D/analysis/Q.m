addpath /Users/Kayliec23/Desktop/analysis
addpath /Users/Kayliec23/Desktop/MISSI/

% Constants
Nx = length(delX);
Ny = length(delY);
Nr = length(delR);
%Nt = size(T,4)-1;
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
q = calculateQ(S,T,U,V,W,Nx,Ny,Nr,Nt,dy,dz,g,rho0,f0);
q = iceNaN(q,Nx,Ny,Nr,icetopo,zz);

T = iceNaN(T,Nx,Ny,Nr,icetopo,zz);

figure; 
subplot(2,1,1)
pcolor(YY/m1km,ZZ,squeeze(q(:,:,:,Nt))'); shading flat; colorbar;
colormap redblue;
caxis([-1 1]*max(abs(q(:)))*1e-5)
hold on; plot(10*ones(length(zz)),zz);

subplot(2,1,2)
pcolor(YY/m1km,ZZ,squeeze(T(:,:,:,Nt))'); shading flat; colorbar;
cmocean('ice');
caxis([min(min(min(T(:,:,:,1)))) max(max(max(T(:,:,:,1))))])
hold on; plot(10*ones(length(zz)),zz);














function Q = calculateQ(S,T,U,V,W,Nx,Ny,Nr,Nt,dy,dz,g,rho0,f0)
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
%Qsym = dbdy.*ZETAy; % symmetric
%Qconv = dbdz.*f0; % convective
%Qcent = dbdz.*ZETAz; % centrifugal
Q = dbdy.*ZETAy + (f0 + ZETAz).*dbdz; % total

end
