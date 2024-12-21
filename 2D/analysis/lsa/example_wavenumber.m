%%%
%%% Example to determine max growth rates over a range of
%%% wavenumbers.
%%%

%%% Initialize input parameters
example_parameters;

%%% Define range of wavenumbers and vectors to store results
lam_max = 1*m1km;
Nlam = 100;
Nmodes = 10;
lam = zeros(1,Nlam);
ll = zeros(1,Nlam);
maxgrowth = zeros(Nmodes,Nlam);

%%% Loop over wavenumbers
for n=1:1:Nlam

  disp(n)

  lam(n) = lam_max / n;
  l = 2*pi/lam(n);
  ll(n) = l;
  
  %%% Supply wavenumber as an input parameter to OILS
  params.l = l;
  
  %%% Compute fastest-growing modes
  result = diffusive_instability(params);

  %%% Remove NaNs
  growth = result.growth;
  growth(isnan(growth)) = [];
  
  %%% Extract fastest growth rates
  maxgrowth(:,n) = growth(1:Nmodes);  
  
end

%%% Plot growth rates as a function of wavenumber
figure(10);
semilogy(2*pi./ll(1:end)/m1km,maxgrowth(1,1:end),'o-');
xlabel('Wavelength (km)','interpreter','latex');
ylabel('Growth rate (s$^{-1}$)','interpreter','latex');
set(gca,'FontSize',14);
% set(gca,'ylim',[0 3]*1e-6);
