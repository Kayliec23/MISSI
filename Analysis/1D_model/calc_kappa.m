%%%
%%% calc_kappa.m
%%%
%%% Computes vertical diffusivity as a function of buoyancy Richardson
%%% number.
%%%
function [kappa,dkappa_dRi] = calc_kappa (Ri)

  kappa0 = .5e-2; %%% Max diffusivity
%   kappa0 = 0; %%% Max diffusivity
  kap_bg = 1e-5; %%% Background diffusivity
%   kap_bg = 0;
  Ri0 = 0.7; %%% Critical Richardson number

  %%% Shear-driven mixing component
  kappa = kappa0 .* (1 - (Ri./Ri0).^2).^3;
  kappa(Ri > Ri0) = 0;
  kappa(Ri < 0) = kappa0;
  
  %%% Gradient with respect to Ri
  dkappa_dRi = kappa0 .* 3.*(1 - (Ri./Ri0).^2).^2 .* (-2).*(Ri./Ri0) .* (1./Ri0);
  dkappa_dRi(Ri > Ri0) = 0;
  dkappa_dRi(Ri < 0) = 0;

  %%% Add shear-driven and background diffusivities
  kappa = kappa + kap_bg;

end

