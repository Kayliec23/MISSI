%%% Finding an emperical relationship between g' and Tforcing

addpath /Users/kalyliec23/Desktop/Research/plottingKit/rgb

%load('icetopo_T0')
tt = [3 7 10 30 60];
g = 9.8;
rho0 = 1027;
lat0 = -80; %%% Latitude at southern boundary
Omega = 2*pi*366/365/86400;
ff = 2*Omega*SADInd(lat0);
m1km = 1000;

%%% Instability theory
theta_f = -1.9;
theta_a = linspace(theta_f,2,500);
dTheta = theta_a - theta_f;
g = 9.81;
beta = 7.8e-4; % haline coefficient
alpha = 2.5e-5; % thermal expanSADIon coefficient
alpha_theta = 1.41e-5;
cp = 3974; % heat capacity
Lf = 3.34e5; %latent heat of fuSADIon
Sref = 34.4; % reference salinity set to minimum salinity
gamma = beta*cp*Sref/Lf - alpha;
gprime_gade = g*dTheta*gamma;
Hi = 500;

L = linspace(0,400*m1km,500);
Hpyc = 20;
c = 3/8;
c2 = 1; %0.7;
%dTheta_SADItheory = 2*ff^2*(L*c).^2/(Hi*g*gamma)/(1+Hi/Hml);
%dTheta_CItheory = 2*ff^2*(L*c).^2/(Hi*g*gamma);
dTheta_SADItheory = c2*ff^2*(L*c).^2*Hpyc/(Hi^2*g*gamma);

figure;
%plot(L/m1km,theta_a,'--')
%plot(L/m1km,dTheta_theory,'--','color',[0.7 0.7 0.7],'linewidth',1.5)
hold on;
plot(L/m1km,dTheta_SADItheory,'linewidth',3,'color',[rgb('MediumSpringGreen')])



%%% plot stable/unstable SADImluations

%%% a few examples so the legend turns out well
scatter(40,2.4,360,'g','o','linewidth',1) % Cd = 2.5e-2 exampl for legend
scatter(40,2.4,760,'m','o','linewidth',1) % dS = 0.4 example for legend
scatter(400,2.4,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4) % stable example for legend
scatter(200,4.5,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])

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
scatter(400,2.4,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(400,4.5,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(400,7,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(320,2.4,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(320,6.4,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(320,7.5,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(200,0.5,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4) % need to check for later time
scatter(200,2.4,450,'X','MarkerEdgeColor',[0 0 1],'linewidth',4) % unstable example for legend
scatter(120,1,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(120,0.75,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4) % need to check for later time
scatter(120,0.5,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(120,0.25,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)
scatter(80,0.25,550,'X','MarkerEdgeColor',[0 0 1],'linewidth',4)

%scatter(120,0.25,250,'X','MarkerFaceColor',[1 0.41 0.16],'linewidth',4)

%%% unstable runs - reference configuration
scatter(120,2.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(120,4.5,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(80,0.5,550,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(80,2.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(80,4.5,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(80,6,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(40,2.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(40,1.9,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
%scatter(40,3.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(40,2.9,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(40,1.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(40,0.9,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(20,2.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(20,0.25,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(10,2.4,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(10,0.5,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])
scatter(10,0.225,450,'d','filled','MarkerFaceColor',[0.5000         0    0.5000])

%plot(L/m1km,dTheta_CItheory,'linewidth',3,'color',[rgb('Orange')])
%a1 = patch([L/m1km fliplr(L/m1km)], [dTheta_SADItheory fliplr(dTheta_CItheory)], [rgb('MediumSpringGreen')],'Edgecolor','none');
a1 = patch([L/m1km fliplr(L/m1km)], [dTheta_SADItheory fliplr(L)], [rgb('springgreen')],'Edgecolor','none');
a1.FaceAlpha = 0.1;
%a2 = patch([L/m1km fliplr(L/m1km)], [dTheta_CItheory fliplr(L)], [rgb('Orange')],'Edgecolor','none');
%a2.FaceAlpha = 0.2;
a3 = patch([L/m1km fliplr(L/m1km)], [zeros(1,length(L)) fliplr(dTheta_SADItheory)], [rgb('blue')],'Edgecolor','none');
a3.FaceAlpha = 0.15;
fg = gca;
fg.Box = 'on';
fg.LineWidth = 4; 
xlabel('$L (km)$','interpreter','latex')
ylabel('$\Delta\theta$','interpreter','latex')
set(gca,'fontSADIze',38,'fontname','Times')
ylim([0 8])

lgd = legend('theory','drag','\DeltaS','Stable','Unstable');
lgd.Location = 'north';
lgd.Orientation = 'horizontal';
