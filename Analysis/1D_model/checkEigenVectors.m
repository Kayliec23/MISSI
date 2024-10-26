
%%
%%% Checks eigenvector magnitudes for a specific instability problem
%%%
%%% Parameters
lamY = 10.^3;
lamZ = 20;
N2 = 2e-4; %%% Buoyancy frequency
l = 2*pi./lamY; %%% Wavenumber vector
m = 2*pi./lamZ;
f = -1e-4;
s = 1e-2;
Ri0 = f^2./(s^2*N2);
[K0,KRi] = calc_kappa (Ri0);
%%% Solution vector [u,v,w,b,phi]
A = [(K0-2*KRi*Ri0)*m^2,    -f-s^2*N2/f,    s*N2/f,     (s/f)*KRi*Ri0*m^2,      0;
     f,                     m^2*K0,         0,          0,                      1i*(l-s*m); ...
     0,                     0,              0,          -1,                     1i*m; ...
     -2*(f/s)*KRi*Ri0*m^2,  -s*N2,          N2,         (K0+KRi*Ri0)*m^2,       0;
     0,                     l-s*m,          m,          0,                      0];
B = [1,     0,    0,    0,    0; ...
     0,     1,    0,    0,    0; ...
     0,     0,    0,    0,    0; ...
     0      0,    0,    1,    0;
     0      0,    0,    0,    0];
[V,D] = eig(A,B);
eigvals = diag(D);
eigvals(isinf(eigvals)) = [];
gr = max(real(-eigvals));
idx = find(real(-eigvals)==gr);
eigvec = V(:,idx);
yy = [-5000:10:5000];
zz = [-50:0.1:50];
[ZZ,YY] = meshgrid(zz,yy);
bb_ref = N2*(ZZ-s*YY);
bb = real(eigvec(4)*exp(1i*(l*YY+m*(ZZ-s*YY))));
vv = real(eigvec(2)*exp(1i*(l*YY+m*(ZZ-s*YY))));
uu = real(eigvec(1)*exp(1i*(l*YY+m*(ZZ-s*YY))));
psi = real(1i/m * eigvec(2)*exp(1i*(l*YY+m*(ZZ-s*YY))));

Bz = N2;
Uz = s.*Bz/f;
Uy = -s.*Uz; 
q = (f-Uy)*bz + (f-Uy-uy).*Bz + uz.*By - Uz,*by + Uz.*By

figure(1);
pcolor(YY,ZZ,bb);
shading interp;
colorbar;
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref,10,'EdgeColor','k');
contour(YY,ZZ,psi,10,'EdgeColor',[.3 .3 .3]);
hold off;
caxis([-1 1]*max(bb(:))*2);
figure(2);
pcolor(YY,ZZ,vv);
shading interp;
colorbar;
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref,10,'EdgeColor','k');
contour(YY,ZZ,psi,10,'EdgeColor',[.3 .3 .3]);
hold off;
caxis([-1 1]*max(vv(:))*2);
figure(3);
pcolor(YY,ZZ,bb);
shading interp;
colorbar;
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb/10,10,'EdgeColor','k');
contour(YY,ZZ,psi,10,'EdgeColor',[.3 .3 .3]);
hold off;
% caxis([-1 1]*max(bb(:))*2);