%%% Finding an emperical relationship between g' and Tforcing

addpath /Users/Kayliec23/Desktop/MISSI/experiments/Tforcing/T0
addpath /Users/Kayliec23/Desktop/MISSI/experiments/Tforcing/T1
addpath /Users/Kayliec23/Desktop/MISSI/experiments/Tforcing/T1.5
addpath /Users/Kayliec23/Desktop/MISSI/experiments/Tforcing/T2
addpath /Users/Kayliec23/Desktop/MISSI/experiments/Tforcing/Tneg0.5
addpath /Users/Kayliec23/Desktop/MISSI/experiments/Tforcing/Tneg1
addpath /Users/Kayliec23/Desktop/MISSI/experiments/ref2D
addpath /Users/Kayliec23/Desktop/plottingKit/rgb

%load('icetopo_T0')
tt = [3 7 10 30 60];
g = 9.8;
rho0 = 1027;
lat0 = -80; %%% Latitude at southern boundary
Omega = 2*pi*366/365/86400;
ff = 2*Omega*sind(lat0);


%%% Instability theory
theta_f = -1.9;
theta_a = linspace(theta_f,2,500);
dTheta = theta_a - theta_f;
g = 9.81;
beta = 7.8e-4; % haline coefficient
alpha = 2.5e-5; % thermal expansion coefficient
alpha_theta = 1.41e-5;
cp = 3974; % heat capacity
Lf = 3.34e5; %latent heat of fusion
Sref = 34.4; % reference salinity set to minimum salinity
gamma = beta*cp*Sref/Lf - alpha;
gprime_gade = g*dTheta*gamma;
Hi = 500;

L = linspace(0,400*m1km,500);
Hpyc = 10;
c = 3/8;
c2 = 1; %0.75;
%dTheta_SItheory = 2*ff^2*(L*c).^2/(Hi*g*gamma)/(1+Hi/Hml);
%dTheta_CItheory = 2*ff^2*(L*c).^2/(Hi*g*gamma);
dTheta_SItheory = c2*ff^2*(L*c).^2*Hpyc/(Hi^2*g*gamma);

figure;
%plot(L/m1km,theta_a,'--')
%plot(L/m1km,dTheta_theory,'--','color',[0.7 0.7 0.7],'linewidth',1.5)
hold on;
plot(L/m1km,dTheta_SItheory,'linewidth',3,'color',[rgb('MediumSpringGreen')])



%%% plot stable/unstable simluations

%%% a few examples so the legend turns out well
scatter(40,2.4,360,'g','o','linewidth',1) % Cd = 2.5e-2 exampl for legend
scatter(40,2.4,760,'m','o','linewidth',1) % dS = 0.4 example for legend
scatter(400,2.4,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4) % stable example for legend
scatter(200,2.4,450,'d','filled','MarkerFaceColor',[0.5000 0 0.5000]) % unstable example for legend


%%% Other configurations
% drag
scatter(40,2.4,460,'g','o','linewidth',1) % Cd = 7.5e-3
scatter(40,2.4,560,'g','o','linewidth',1) % Cd = 7.5e-4
scatter(40,2.4,660,'g','o','linewidth',1) % Cd = 2.5e-4
% dS
scatter(40,2.4,860,'m','o','linewidth',1) % dS = 0.6
scatter(40,2.4,960,'m','o','linewidth',1) % dS = 0.8
scatter(40,2.4,1060,'m','o','linewidth',1) % dS = 1


%%% stable runs - reference configuration
%scatter(400,4.4,250,'X','MarkerFaceColor',[1 0.41 0.16],'linewidth',2)
scatter(320,2.4,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
%scatter(320,4.4,250,'X','MarkerFaceColor',[1 0.41 0.16],'linewidth',2)
scatter(200,0.5,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4) % need to check for later time
scatter(120,0.75,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4) % need to check for later time
scatter(120,0.5,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(120,0.25,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(80,0.25,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(80,0.5,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
%scatter(120,0.25,250,'X','MarkerFaceColor',[1 0.41 0.16],'linewidth',4)

%%% unstable runs - reference configuration
scatter(320,4.5,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(200,4.5,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(120,1,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(120,2.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(120,4.5,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(80,2.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(80,4.5,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(80,6.3,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(80,2.15,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(40,2.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(40,1.9,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(40,3.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(40,2.9,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(40,1.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(40,0.9,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(10,2.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(20,2.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(10,0.5,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(10,0.25,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])

%plot(L/m1km,dTheta_CItheory,'linewidth',3,'color',[rgb('Orange')])
%a1 = patch([L/m1km fliplr(L/m1km)], [dTheta_SItheory fliplr(dTheta_CItheory)], [rgb('MediumSpringGreen')],'Edgecolor','none');
a1 = patch([L/m1km fliplr(L/m1km)], [dTheta_SItheory fliplr(L)], [rgb('MediumSpringGreen')],'Edgecolor','none');
a1.FaceAlpha = 0.1;
%a2 = patch([L/m1km fliplr(L/m1km)], [dTheta_CItheory fliplr(L)], [rgb('Orange')],'Edgecolor','none');
%a2.FaceAlpha = 0.2;
a3 = patch([L/m1km fliplr(L/m1km)], [zeros(1,length(L)) fliplr(dTheta_SItheory)], [rgb('Blue')],'Edgecolor','none');
a3.FaceAlpha = 0.15;
fg = gca;
fg.Box = 'on';
fg.LineWidth = 4; 
xlabel('$L (km)$','interpreter','latex')
ylabel('$\Delta\theta$','interpreter','latex')
set(gca,'fontsize',38,'fontname','Times')
ylim([0 8])

[h,icons] = legend('theory','drag','\DeltaS','Stable','Unstable');



% %%
% figure;
% for idx = 1:length(tt)
%     for d = [2 3 8 20]
%         t_id = tt(idx);
%         y_loc = 250;
%         ice_dist = d;
%         min_zloc = round(abs(icetopo(y_loc)-ice_dist));
%         max_zloc = 250;
%         
%         load('Diag_T0');
%         for n = 1:Nt
%             SS = squeeze(S(1,:,:,n));
%             TT = squeeze(T(1,:,:,n));
%             rho(:,:,n) = densjmd95(SS,TT,0);
%         end
%         dRho = rho(y_loc,max_zloc,t_id) - rho(y_loc,min_zloc,t_id);
%         gprime0 = g/rho0*dRho;
%         
%         load('Diag_T1');
%         for n = 1:Nt
%             SS = squeeze(S(1,:,:,n));
%             TT = squeeze(T(1,:,:,n));
%             rho(:,:,n) = densjmd95(SS,TT,0);
%         end
%         dRho = rho(y_loc,max_zloc,t_id) - rho(y_loc,min_zloc,t_id);
%         gprime1 = g/rho0*dRho;
%         
%         load('Diag_T15');
%         for n = 1:Nt
%             SS = squeeze(S(1,:,:,n));
%             TT = squeeze(T(1,:,:,n));
%             rho(:,:,n) = densjmd95(SS,TT,0);
%         end
%         dRho = rho(y_loc,max_zloc,t_id) - rho(y_loc,min_zloc,t_id);
%         gprime15 = g/rho0*dRho;
%         
%         load('Diag_Tneg05');
%         for n = 1:Nt
%             SS = squeeze(S(1,:,:,n));
%             TT = squeeze(T(1,:,:,n));
%             rho(:,:,n) = densjmd95(SS,TT,0);
%         end
%         dRho = rho(y_loc,max_zloc,t_id) - rho(y_loc,min_zloc,t_id);
%         gprime_neg05 = g/rho0*dRho;
%         
%         load('Diag_Tneg1');
%         for n = 1:Nt
%             SS = squeeze(S(1,:,:,n));
%             TT = squeeze(T(1,:,:,n));
%             rho(:,:,n) = densjmd95(SS,TT,0);
%         end
%         dRho = rho(y_loc,max_zloc,t_id) - rho(y_loc,min_zloc,t_id);
%         gprime_neg1 = g/rho0*dRho;
%         
%         load('diag_ref2D');
%         for n = 1:Nt
%             SS = squeeze(S(1,:,:,n));
%             TT = squeeze(T(1,:,:,n));
%             rho(:,:,n) = densjmd95(SS,TT,0);
%         end
%         dRho = rho(y_loc,max_zloc,t_id) - rho(y_loc,min_zloc,t_id);
%         gprime05 = g/rho0*dRho;
%         
%         x = [-1 -0.5 0 0.5 1 1.5];
%         y = [gprime_neg1 gprime_neg05 gprime0 gprime05 gprime1 gprime15];
%         
%         subplot(1,length(tt),idx)
%         scatter(x,y,'filled')
%         hold on;
%         
%     end
%     
%     %%% Gade line
%     theta_f = -1.9;
%     theta_a = linspace(theta_f,2,500);
%     dTheta = theta_a - theta_f;
%     g = 9.81;
%     beta = 7.8e-4; % haline coefficient
%     alpha = 2.5e-5; % thermal expansion coefficient
%     alpha_theta = 1.41e-5;
%     cp = 3974; % heat capacity
%     Lf = 3.34e5; %latent heat of fusion
%     Sref = 34.4; % reference salinity set to minimum salinity
%     gamma = beta*cp*Sref/Lf - alpha;
%     gprime_gade = g*dTheta*gamma;
%     
%     Hi = 500;
%     c = 3/8;
%     L = linspace(0,400*m1km,500);
%     gprime_theory = ff^2.*L.^2*c^2/Hi;
%     
%     plot(theta_a,gprime_gade,'linewidth',2,'color',[0.7 0.7 0.7])
%     plot(theta_a,0.65*gprime_gade,'--','linewidth',1.5,'color',[0.7 0.7 0.7])
%     %plot(theta_a,gprime_theory)
% 
%     title(['t = ',num2str(tt(idx)),'day'])
%     set(gca,'fontsize',20)
%     ylabel('g^{\prime} (ms^{-2})')
%     xlabel('T_{forcing} (^oC)')
%     
% end
% 
% legend('z = 2m', 'z = 3m', 'z = 8m', 'z = 20m', 'theory', 'theory, C = 0.65')

