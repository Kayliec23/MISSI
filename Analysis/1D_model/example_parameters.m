%%%
%%% Example input parameter setup.
%%%

%%% Initialize input parameter structure
params = struct;

%%% Scalar parameters
m1km = 1000;
Lz = 100;
M = 200; %%% Sufficient resolution to test examples
f = -1.4e-4;
s = 1e-2;

%%% Add parameters to input parameter structure
params.Lz = Lz;
params.f = f;
params.M = M;
params.s = s;

%%% Grids on h and u points
dz = params.Lz/params.M;
zz = -Lz/2+0.5*dz:dz:Lz/2-0.5*dz;

%%% Exponential pycnocline
bz0 = 1e-6; %%% Ambient stratification
bzmax = 6e-4; %%% Max stratification
% bzmax = 0;
Hpyc = 20; %%% Pycnocline e-folding scale
params.Bz = (bz0 + bzmax * exp(-(zz./Hpyc).^2));