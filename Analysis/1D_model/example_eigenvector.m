%%%
%%% Example use of OILS to extract eigenvector corresponding to
%%% fastest-growing mode.
%%%

%%% Initialize input parameter structure and set up bathymetry and
%%% along-slope velocity profiles
example_parameters;

%%% Have OILS calculate fastest growth rates
params.incrop_no_flux = false;
params.return_growth = true;
params.return_eivecs = true;
params.adjust_eta = true;
  
%%% Supply wavenumber as an input parameter to OILS
lam = 20*m1km;
l = 2*pi/lam; 
params.l = l;
  
%%% Compute fastest-growing modes
result = oils(params);
  
%%% Extract eigenvector
V = result.V;
D = result.D;
D = D./1i;
D = (D - real(D))*1i;
[m,n] = find(D==max(max(D)));
vec = V(:,n)./1i; 

%%% Extract grids
xx1_h = result.xx1_h;
xx2_h = result.xx2_h;
xx1_u = result.xx1_u;
xx2_u = result.xx2_u;
M = result.M;
N = result.N;

%%% Plot vbar
figure(1);
plot(result.xx1_h/m1km,result.v1bar);
hold on;
plot(result.xx2_h/m1km,result.v2bar);
hold off;
xlabel('$x$ (km)','interpreter','latex');
ylabel('$\overline{v}$ (m/s)','interpreter','latex');
handle = legend('$\overline{v}_1$','$\overline{v}_2$');
set(handle,'interpreter','latex');
set(gca,'FontSize',14);

%%% Plot etabar
figure(2);
plot(result.xx1_h/m1km,result.eta1bar);
hold on;
plot(result.xx2_h/m1km,result.eta2bar);
plot(result.xx1_h/m1km,result.etab);
hold off;
xlabel('$x$ (km)','interpreter','latex');
ylabel('$z$ (m)','interpreter','latex');
handle = legend('$\overline{\eta}_1$','$\overline{\eta}_2$','$\eta_b$');
set(handle,'interpreter','latex');
set(gca,'FontSize',14);

%%% Plot u1^
figure(3);
plot(xx1_u/m1km,abs(vec(1:M-1)));
xlabel('$x$ (km)','interpreter','latex');
ylabel('$|\hat{u}_1|$','interpreter','latex');
set(handle,'interpreter','latex','Location','NorthWest');
set(gca,'FontSize',14);
set(gca,'xlim',[50 120]);

%%% Plot u2^
figure(4);
plot(xx2_u/m1km,abs(vec(3*M-1+1:3*M-1+N-1)));
xlabel('$x$ (km)','interpreter','latex');
ylabel('$|\hat{u}_2|$','interpreter','latex');
set(handle,'interpreter','latex','Location','NorthWest');
set(gca,'FontSize',14);
set(gca,'xlim',[50 120]);

%%% Plot v1^
figure(5);
plot(xx1_h/m1km,abs(vec(M:2*M-1)));
xlabel('$x$ (km)','interpreter','latex');
ylabel('$|\hat{v}_1|$','interpreter','latex');
set(handle,'interpreter','latex','Location','NorthWest');
set(gca,'FontSize',14);
set(gca,'xlim',[50 120]);

%%% Plot v2^
figure(6);
plot(xx2_h/m1km,abs(vec(3*M-1+N:3*M-1+2*N-1)));
xlabel('$x$ (km)','interpreter','latex');
ylabel('$|\hat{v}_2|$','interpreter','latex');
set(handle,'interpreter','latex','Location','NorthWest');
set(gca,'FontSize',14);
set(gca,'xlim',[50 120]);

%%% Plot eta2^
figure(7);
plot(xx2_h/m1km,abs(vec(3*M-1+2*N:3*M-1+3*N-1)));
xlabel('$x$ (km)','interpreter','latex');
ylabel('$|\hat{\eta}_2|$','interpreter','latex');
set(handle,'interpreter','latex','Location','NorthWest');
set(gca,'FontSize',14);
set(gca,'xlim',[50 120]);
