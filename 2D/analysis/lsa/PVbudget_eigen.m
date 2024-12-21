%%
%%% Checks eigenvector magnitudes for a specific instability problem
%%% 
addpath /Users/kalyliec23/Desktop/MISSI/Analysis/DiffusiveInstability-main


%%% check less stable modes for directions

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
%%% Arbitrary amplitude rescaling
eigvec = eigvec / 50;
Uz0 = s*N2/f;
yy = [-5000:10:5000];
zz = [-50:0.1:50];
[ZZ,YY] = meshgrid(zz,yy);
bb = real(eigvec(4)*exp(1i*(l*YY+m*(ZZ-s*YY))));
uu = real(eigvec(1)*exp(1i*(l*YY+m*(ZZ-s*YY))));
vv = real(eigvec(2)*exp(1i*(l*YY+m*(ZZ-s*YY))));
ww = real(eigvec(3)*exp(1i*(l*YY+m*(ZZ-s*YY))));
psi = real(1i/m * eigvec(2)*exp(1i*(l*YY+m*(ZZ-s*YY))));
bb_ref = N2*(ZZ-s*YY);

uu_cmplx = eigvec(1)*exp(1i*(l*YY+m*(ZZ-s*YY)));
vv_cmplx = eigvec(2)*exp(1i*(l*YY+m*(ZZ-s*YY)));
ww_cmplx = eigvec(3)*exp(1i*(l*YY+m*(ZZ-s*YY)));
bb_cmplx = eigvec(4)*exp(1i*(l*YY+m*(ZZ-s*YY)));

Ri = Ri0 + real(1i*m*bb_cmplx / Uz0^2 - 2*N2/Uz0^3 * 1i*m*uu_cmplx);
Uy0 = - s*Uz0;
M2 = - s*N2;
q0 = (f - Uy0)*N2 + Uz0*M2;
q = q0 + real((f-Uy0)*1i*m*bb_cmplx - 1i*l*uu_cmplx*N2 + 1i*m*M2 + Uz0*1i*l*bb_cmplx);
M = uu-f*YY; 


% - du/dt = (K0-2*KRi*Ri0)*m^2 *u   + (-f-s^2*N2/f) * v   + s*N2/f * w   +   (s/f)*KRi*Ri0*m^2 * b
% - db/dt = -2*(f/s)*KRi*Ri0*m^2 * u  +  -s*N2 * v    +    N2 * w     +    (K0+KRi*Ri0)*m^2 * b
% dq/dt = (f-Uy0)*1i*m * db/dt - 1i*l*N2 * du/dt + 1i*m*M2 * du/dt +  Uz0*1i*l * db/dt
dudt_visc = - (K0-2*KRi*Ri0)*m^2 * uu_cmplx      -   (s/f)*KRi*Ri0*m^2 * bb_cmplx;
dudt_adv = - (-f-s^2*N2/f) * vv_cmplx   - s*N2/f * ww_cmplx;
dbdt_diff = 2*(f/s)*KRi*Ri0*m^2 * uu_cmplx     -   (K0+KRi*Ri0)*m^2 * bb_cmplx;
dbdt_adv = s*N2 * vv_cmplx    -   N2 * ww_cmplx;
dqdt_visc = real(- 1i*(l-s*m)*N2 * dudt_visc + 1i*m*M2 * dudt_visc);
dqdt_diff = real((f-Uy0)*1i*m * dbdt_diff +  Uz0*1i*(l-s*m) * dbdt_diff);
dqdt_adv = real((f-Uy0)*1i*m * dbdt_adv - 1i*(l-s*m)*N2 * dudt_adv + 1i*m*M2 * dudt_adv +  Uz0*1i*(l-s*m) * dbdt_adv);
dqdt_tot = dqdt_adv + dqdt_diff + dqdt_visc;

%%
%%%%%%%%%%%%% figure of PV budget  %%%%%%%%%%%%%
fontsize = 20;
gap = [0.05 0.02];
marg_h = [0.09 0.1];
marg_w = [0.14 0.1];

figure;
subtightplot(3,1,1,gap,marg_h,marg_w)
pcolor(YY,ZZ,dqdt_tot);
shading interp;
colorbar('location','northoutside');
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,psi,10,'EdgeColor',[.3 .3 .3]);
hold off;
xticklabels([])
caxis([-1 1]*max(abs(dqdt_diff(:)))*2);
set(gca,'fontsize',fontsize)
ylabel('z (m)')
title('Total PV tendency','interpreter','latex');

subtightplot(3,1,2,gap,marg_h,marg_w)
pcolor(YY,ZZ,dqdt_visc);
shading interp;
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,psi,10,'EdgeColor',[.3 .3 .3]);
hold off;
xticklabels([])
set(gca,'fontsize',fontsize)
ylabel('z (m)')
caxis([-1 1]*max(abs(dqdt_diff(:)))*2);
title('Viscous PV tendency','interpreter','latex');

subtightplot(3,1,3,gap,marg_h,marg_w)
pcolor(YY,ZZ,dqdt_diff);
shading interp;
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,psi,10,'EdgeColor',[.3 .3 .3]);
hold off;
set(gca,'fontsize',fontsize)
ylabel('z (m)')
xlabel('y (m)')
caxis([-1 1]*max(abs(dqdt_diff(:)))*2);
title('Diffusive PV tendency','interpreter','latex');
%%
%%%%%%%%%%%%% figure of b,u,v,w,q,Ri  %%%%%%%%%%%%%
gap = [0.06 0.07];
marg_h = [0.1 0.05];
marg_w = [0.07 0.08];
fontsize = 28;

figure;
subtightplot(2,3,4,gap,marg_h,marg_w)
pcolor(YY,ZZ,bb);
shading interp;
colorbar('location','eastoutside','linewidth',1.5,'TickLabelInterpreter','latex');
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,Ri,10,'EdgeColor',[.3 .3 .3]);
hold off;
caxis([-1 1]*max(bb(:))*2);
ylabel('z (m)','Interpreter','latex')
xlabel('y (km)','Interpreter','latex')
set(gca,'fontsize',fontsize,'linewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
title('(d) b','interpreter','latex')

subtightplot(2,3,2,gap,marg_h,marg_w)
pcolor(YY,ZZ,vv);
shading interp;
colorbar('location','eastoutside','linewidth',1.5,'TickLabelInterpreter','latex');
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,Ri,10,'EdgeColor',[.3 .3 .3]);
hold off;
xticklabels([])
yticklabels([])
caxis([-1 1]*max(vv(:))*2);
set(gca,'fontsize',fontsize,'linewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
title('(b) v','interpreter','latex')

subtightplot(2,3,1,gap,marg_h,marg_w)
pcolor(YY,ZZ,uu);
shading interp;
colorbar('location','eastoutside','linewidth',1.5,'TickLabelInterpreter','latex');
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,Ri,10,'EdgeColor',[.3 .3 .3]);
hold off;
xticklabels([])
ylabel('z (m)','Interpreter','latex')
caxis([-1 1]*max(uu(:))*2);
set(gca,'fontsize',fontsize,'linewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
title('(a) u','interpreter','latex')

subtightplot(2,3,5,gap,marg_h,marg_w);
pcolor(YY,ZZ,q);
shading interp;
colorbar('location','eastoutside','linewidth',1.5,'TickLabelInterpreter','latex');
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,Ri,10,'EdgeColor',[.3 .3 .3]);
hold off;
xlabel('y (km)','Interpreter','latex')
yticklabels([])
caxis([-1 1]*max(abs(q(:)))*1.25);
set(gca,'fontsize',fontsize,'linewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
title('(e) q','interpreter','latex')

subtightplot(2,3,3,gap,marg_h,marg_w)
pcolor(YY,ZZ,ww);
shading interp;
colorbar('location','eastoutside','LineWidth',1.5,'TickLabelInterpreter','latex');
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,Ri,10,'EdgeColor',[.3 .3 .3]);
hold off;
xticklabels([])
yticklabels([])
caxis([-1 1]*max(ww(:))*2);
set(gca,'fontsize',fontsize,'linewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
title('(c) w','interpreter','latex')

sb = subtightplot(2,3,6,gap,marg_h,marg_w)
pcolor(YY,ZZ,Ri);
shading interp;
colorbar('location','eastoutside','linewidth',1.5,'TickLabelInterpreter','latex');
colormap(sb,'pink');
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,psi,10,'EdgeColor',[.3 .3 .3]);
hold off;
xlabel('y (km)','Interpreter','latex')
yticklabels([])
set(gca,'fontsize',fontsize,'linewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
title('(f) Ri','interpreter','latex')
caxis([0 1])



%%% meridoional mom balance
figure;

subplot(2,2,1)
pcolor(YY,ZZ,vv);
shading interp;
colorbar;
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,psi,10,'EdgeColor',[.3 .3 .3]);
hold off;
title('v')

subplot(2,2,2)
dvdt_cor = real(-f*eigvec(1)*exp(1i*(l*YY+m*(ZZ-s*YY))));
pcolor(YY,ZZ,dvdt_cor);
shading interp;
colorbar;
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,psi,10,'EdgeColor',[.3 .3 .3]);
hold off;
caxis([-1 1]*max(abs(dvdt_cor(:)))*2);
title('V tendency due to Coriolis');

subplot(2,2,3)
dvdt_visc = real(-m^2*K0*eigvec(2)*exp(1i*(l*YY+m*(ZZ-s*YY))));
pcolor(YY,ZZ,dvdt_visc);
shading interp;
colorbar;
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,psi,10,'EdgeColor',[.3 .3 .3]);
hold off;
caxis([-1 1]*max(abs(dvdt_visc(:)))*2);
title('V tendency due to viscosity');

subplot(2,2,4)
dvdt_pgf = real(-1i*(l-s*m)*eigvec(5)*exp(1i*(l*YY+m*(ZZ-s*YY))));
pcolor(YY,ZZ,dvdt_pgf);
shading interp;
colorbar;
colormap redblue;
hold on;
contour(YY,ZZ,bb_ref+bb,10,'EdgeColor','k');
contour(YY,ZZ,psi,10,'EdgeColor',[.3 .3 .3]);
hold off;
caxis([-1 1]*max(abs(dvdt_pgf(:)))*2);
title('V tendency due to pressure gradient');
