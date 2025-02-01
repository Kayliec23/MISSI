%%i%
%%% setParams.m
%%%
%%% Sets basic MITgcm parameters plus parameters for included packages, and
%%% writes out the appropriate input files.,
%%%
function nTimeSteps = setParams (inputpath,codepath,listterm,Nx,Ny,Nr)  

  %%% Load EOS utilities
  addpath /data3/kayliec23/MITgcm_IP/GSW;
  addpath /data3/kayliec23/MITgcm_IP/GSW/html;
  addpath /data3/kayliec23/MITgcm_IP/GSW/library;
  addpath /data3/kayliec23/MITgcm_IP/GSW/pdf;
  addpath /data3/kayliec23/MITgcm_IP/utils/matlab
  
  %%% set EOS type (linear, jmd95z, mjdwf)
  eos = 'JMD95Z';
    
  
  %%%%%%%%%%%%%%%%%%
  %%%%% SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%      
  
  %%% If set true, plots of prescribed variables will be shown
  showplots = true;      
  fignum = 1;
  
  %%% Data format parameters
  ieee='b';
  prec='real*8';
  realdigits = 8;
  realfmt=['%.',num2str(realdigits),'e'];
  
  %%% Get parameter type definitions
  paramTypes;
  
  %%% To store parameter names and values
  parm01 = parmlist;
  parm02 = parmlist;
  parm03 = parmlist;
  parm04 = parmlist;
  parm05 = parmlist;
  PARM={parm01,parm02,parm03,parm04,parm05}; 
  
  %%% Seconds in one hour
  t1min = 60;
  %%% Seconds in one hour
  t1hour = 60*t1min;
  %%% hours in one day
  t1day = 24*t1hour;
  %%% Seconds in 1 year
  t1year = 365*t1day;  
  %%% Metres in one kilometre
  m1km = 1000; 
  %%% Pascal in 1 dbar
  Pa1dbar = 1e4;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SIMULATION CONTROL PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  nonHydrostatic = true; %%% Whether to run in nonhydrostatic mode
  use_seaIce = false; 
  use_shelfIce = true; 
  use_RBCS = false;
  use_OBCS = true; 
  use_3D = true; %%% Whether to run a 3D vs 2D simulation
  use_KPP =true;
  use_GGL90 = false; % ----> added. For now set to false ; this is a new addition to isomip exepriment, maybe worth checking out?
  %useOrlanskiNorth = false;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% FIXED PARAMETER VALUES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  
  simTime = 30*t1day; %%% Simulation time  
  nIter0 = 0; %%% Initial iteration 
  if (use_3D)
    Lx = 40*m1km; %%% Domain size in x 
  else
    Lx = 1*m1km; %%% Domain size in x 
  end
  Ly = 180*m1km; %%% Domain size in y  ----> depends on defined variable grid resolution. Better way to do this? 
  H = 700; %%% Domain size in z 
  g = 9.81; %%% Gravity
  Omega = 2*pi*366/365/86400;  
  lat0 = -80; %%% Latitude at southern boundary
  f0 = 2*Omega*sind(lat0); %%% Coriolis parameter      
  rho0 = 1027;
  cp = 3974.0; %%% specific heat capacity for the ocean (J/kg/K) ------------- ADDED
  
  %%% Diffusion parameters
  viscAh = 0; %%% 
  viscA4 = 30000; %%% Biharmonic viscosity ------> does this number even make sense? Also, should it be viscA4W since it'gs NH?
  viscAhGrid = 0; %%% Grid-dependent viscosity
%   viscA4Grid = 0.2; %%% Grid-dependent biharmonic viscosity    
%   viscC4smag = 0; %%% Smagorinsky biharmonic viscosity
%   diffK4Tgrid = 0.2; %%% Grid-dependent biharmonic diffusivity
  viscA4Grid = 0; %%% Grid-dependent biharmonic viscosity -------> Should I be using this instead of viscA4? I don't quite understand the difference between them, I thought biharmonic viscosity is grid dependent    
  viscC4smag = 0; %%% Smagorinsky biharmonic viscosity
  diffK4Tgrid = 0; %%% Grid-dependent biharmonic diffusivity
  viscAr = 1.67e-4; %%% Vertical viscosity
  diffKhT = 0; %%% Horizontal temp diffusion
  diffKrT = 0; %%% Vertical temp diffusion
  diffKhS = 0; %%% Horizontal salt diffusion
  diffKrS = 0; %%% Vertical salt diffusion
  
  %%% Parameters related to periodic forcing
  periodicExternalForcing = false;
  externForcingCycle = 2*simTime;  
  if (~periodicExternalForcing)
    externForcingPeriod = externForcingCycle;
    nForcingPeriods = 1;
  else
    externForcingPeriod = 0.01*t1day;
    nForcingPeriods = externForcingCycle/externForcingPeriod;
  end
  
  %%% Partial cell capabilities ------------------- ADDED
  hFacMin = 0.2;
  hFacMinDr = 0.1;
  hFacSup = 2.0;
  hFacInf = 0.02;
  
  
  %%% PARM01
  %%% momentum scheme
  parm01.addParm('vectorInvariantMomentum',true,PARM_BOOL);
  %%% viscosity  
  parm01.addParm('viscAr',viscAr,PARM_REAL);
  parm01.addParm('viscA4',viscA4,PARM_REAL);
  parm01.addParm('viscAh',viscAh,PARM_REAL);
  parm01.addParm('viscA4Grid',viscA4Grid,PARM_REAL);
  parm01.addParm('viscAhGrid',viscAhGrid,PARM_REAL);
  parm01.addParm('viscA4GridMax',0.5,PARM_REAL);
  parm01.addParm('viscAhGridMax',1,PARM_REAL);
  parm01.addParm('useAreaViscLength',false,PARM_BOOL);
  parm01.addParm('useFullLeith',false,PARM_BOOL);
  parm01.addParm('viscC4leith',0,PARM_REAL);
  parm01.addParm('viscC4leithD',0,PARM_REAL);      
  parm01.addParm('viscC2leith',0,PARM_REAL);
  parm01.addParm('viscC2leithD',0,PARM_REAL); 
  parm01.addParm('viscC2smag',0,PARM_REAL); 
  parm01.addParm('viscC4smag',viscC4smag,PARM_REAL); 
  
  %%% diffusivity
  parm01.addParm('diffKrT',diffKrT,PARM_REAL);
  parm01.addParm('diffKhT',diffKhT,PARM_REAL); 
  parm01.addParm('diffKrS',diffKrS,PARM_REAL);
  parm01.addParm('diffKhS',diffKhS,PARM_REAL); 
  %%% advection schemes
  parm01.addParm('tempAdvScheme',77,PARM_INT);
  parm01.addParm('saltAdvScheme',77,PARM_INT);
  parm01.addParm('multiDimAdvection',true,PARM_BOOL);
  parm01.addParm('tempStepping',true,PARM_BOOL);
  parm01.addParm('saltStepping',true,PARM_BOOL);
  parm01.addParm('staggerTimeStep',true,PARM_BOOL);
  parm01.addParm('convertFW2Salt',33.4,PARM_REAL); %%%%%%%%%%%%%% was 0
  %%% equation of state
  parm01.addParm('eosType',eos,PARM_STR);   
  %%% boundary conditions
  parm01.addParm('no_slip_sides',false,PARM_BOOL);
  parm01.addParm('no_slip_bottom',false,PARM_BOOL);
  parm01.addParm('bottomDragLinear',0,PARM_REAL);
  parm01.addParm('bottomDragQuadratic',2.5e-3,PARM_REAL); %%%%%%%%%%%%%%%% was 0
  parm01.addParm('implicitFreeSurface',true,PARM_BOOL); %%%%%%%%%%%%%%%%%% ADDED
  parm01.addParm('rigidLid',false,PARM_BOOL); %%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDED
  %%% physical parameters
  parm01.addParm('f0',f0,PARM_REAL);
  parm01.addParm('gravity',g,PARM_REAL);
  parm01.addParm('rhonil',rho0,PARM_REAL);
  parm01.addParm('rhoConst',rho0,PARM_REAL);
  parm01.addParm('HeatCapacity_cp',cp,PARM_REAL); %%%%%%%%%%%% ADDED
  %%% partial cell  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDED
  parm01.addParm('hFacMin',hFacMin,PARM_REAL);  
  parm01.addParm('hFacMinDr',hFacMinDr,PARM_REAL);  
  parm01.addParm('hFacSup',hFacSup,PARM_REAL);
  parm01.addParm('hFacInf',hFacInf,PARM_REAL);
 
  %%% implicit diffusion and convective adjustment  
  parm01.addParm('ivdc_kappa',0,PARM_REAL);
  parm01.addParm('implicitDiffusion',true,PARM_BOOL);
  parm01.addParm('implicitViscosity',true,PARM_BOOL);
  parm01.addParm('nonHydrostatic',nonHydrostatic,PARM_BOOL);
  %%% exact volume conservation
  parm01.addParm('exactConserv',true,PARM_BOOL);
  %%% C-V scheme for Coriolis term
  parm01.addParm('useCDscheme',false,PARM_BOOL);  

  %%% file IO stuff
  parm01.addParm('readBinaryPrec',64,PARM_INT);
  parm01.addParm('useSingleCpuIO',true,PARM_BOOL);  
  parm01.addParm('debugLevel',5,PARM_INT); %%%%%%%%%%%%%%%%%%%% was -1

  %%%%% Vertical advection
%   parm01.addParm('momImplVertAdv',true,PARM_BOOL);
%   parm01.addParm('tempImplVertAdv',true,PARM_BOOL);
%   parm01.addParm('saltImplVertAdv',true,PARM_BOOL);


  %%% Wet-point method at boundaries - may improve boundary stability
  parm01.addParm('useJamartWetPoints',true,PARM_BOOL); %%%%%%%%%%%%%%%% was true
  parm01.addParm('useJamartMomAdv',true,PARM_BOOL); %%%%%%%%%%%%%%%%%%% was true
  parm01.addParm('useRealFreshWaterFlux',false,PARM_BOOL)

  %%% PARM02
  parm02.addParm('useSRCGSolver',false,PARM_BOOL); %%%%%%%%%%%%%%% was true  
  parm02.addParm('cg2dMaxIters',1000,PARM_INT);  
  parm02.addParm('cg2dTargetResidual',1e-12,PARM_REAL); %%%%%%%%%%% was 1e-12
  %parm02.addParm('cg3dMaxIters',400,PARM_INT); %%%%%%%%%%%%%%%%%%%% was 300 
  %parm02.addParm('cg3dTargetResidual',1e-13,PARM_REAL); %%%%%%%%%%% was 1e-7  
 
  %%% PARM03
  %parm03.addParm('alph_AB',1/2,PARM_REAL);
  %parm03.addParm('beta_AB',5/12,PARM_REAL);
  %parm03.addParm('momDissip_In_AB',true,PARM_BOOL); %%%%%%%% was false (not specified in isomip, but default is true)
  
  %%%%%%% Added to make SEAICE work
  %parm03.addParm('momForcingOutAB',1,PARM_INT);
  %parm03.addParm('tracForcingOutAB',1,PARM_INT);


  %%%%%%%
  parm03.addParm('nIter0',nIter0,PARM_INT);
  parm03.addParm('abEps',0.1,PARM_REAL);
  parm03.addParm('chkptFreq',1*t1day,PARM_REAL);
  parm03.addParm('pChkptFreq',5*t1day,PARM_REAL);
  parm03.addParm('taveFreq',0,PARM_REAL);
  parm03.addParm('dumpFreq',1*t1day,PARM_REAL);
  parm03.addParm('monitorFreq',1*t1day,PARM_REAL);
  parm03.addParm('cAdjFreq',0,PARM_REAL);
  parm03.addParm('dumpInitAndLast',true,PARM_BOOL);
  parm03.addParm('pickupStrictlyMatch',false,PARM_BOOL);
  if (periodicExternalForcing)
    parm03.addParm('periodicExternalForcing',periodicExternalForcing,PARM_BOOL);
    parm03.addParm('externForcingPeriod',externForcingPeriod,PARM_REAL);
    parm03.addParm('externForcingCycle',externForcingCycle,PARM_REAL);
  end
  
  %%% PARM04
  parm04.addParm('usingCartesianGrid',true,PARM_BOOL);
  parm04.addParm('usingSphericalPolarGrid',false,PARM_BOOL);    
 
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% GRID SPACING %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%    
    

  %%%% Uniform meridional grid   
  %dy = (Ly/Ny)*ones(1,Ny);
  %yy = cumsum(dy);
  %yy = yy-mean(yy);
  
  %%%% Zonal grid
  %if (use_3D)
  %  dx = Lx/Nx*ones(1,Nx);  
  %  xx = cumsum(dx);
  %  xx = xx-mean(xx);
  %else
  %  dx = dy(1);
  %  xx = 0;
  %end
  
  
  
  %%% Variable meridional grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDED
  resMin = 100; % minimum resolution
  resSponge = 1*m1km; % resolution in sponge
  Lmin = 45*m1km; % distance over which resMin applies
  Lsponge = 20*m1km; % distance over which sponge applies %%% should define this variable earlier, since it is used for applying the sponge
  nSponge = Lsponge/resSponge;
  Lmax = Ly-Lmin-Lsponge; % distance over which resMax applies (see below)
  nResMax = Ny - Lmin/resMin-nSponge; % number of grid point for maximum resolution
  resMax = 1000; %Lmax/nResMax; % maximum resolution
  nWidthRes = 30; % number of grid points to transition between resMin and resMax
  A = resSponge; % area A
  %B = Lsponge/resSponge - Lmin/resMin;
  B = resSponge - resMin; % area B
%   nswitch_min = Lmin/resMin; % index change from minRes to larger resolution
%   nswitch_max = Ny - nSponge; % index change to sponge resolution (creating an approximately constant resolution within the sponge)
%   for i= 1:nswitch_max
%     dy(i) = round(resMin + 0.5*B*(tanh((i-nswitch_min)/nWidthRes)+1)); % don't have to round, but was causing some issues in the past
%     %dy(i) = round(Lmin - 0.5*B*(tanh((i-nswitch_min)/nWidthRes)+1));
%   end
%   for i = nswitch_max:Ny
%     dy(i) = resSponge; 
%   end

  
  %%% Better code for grid-spacing. 
  dy_vec = [resMin resMin  resMax resMax];
  yy_vec = [0 Lmin Ly-Lsponge Ly];
  NN_vec = cumsum([1 410 220 9]);
  dy = interp1(NN_vec,dy_vec,1:NN_vec(end),'pchip');
  yy_s = [0 cumsum(dy)];
  yy = 0.5*(yy_s(1:end-1)+yy_s(2:end));
  
  
  %%% Zonal grid
  if (use_3D)
    dx = Lx/Nx*ones(1,Nx);  
    xx_s = [0 cumsum(dx)];
    xx = 0.5*(xx_s(1:end-1)+xx_s(2:end));
    %%% Plotting mesh
    [Y,X] = meshgrid(yy,xx);
  else
    dx = min(dy); %
    xx_s = [0 cumsum(dx)];
    xx = 0.5*(xx_s(1:end-1)+xx_s(2:end));
  end
 
 
 
  
  %%% Grid spacing increases with depth, but spacings exactly sum to H
  zidx = 1:Nr;
  dz = H/Nr*ones(1,Nr);
  zz = -cumsum((dz+[0 dz(1:end-1)])/2);
  
  %%% Heights of cell vertical faces
  zzf = -cumsum([0 dz]);

  %%% Store grid spacings
  parm04.addParm('delX',dx,PARM_REALS);
  parm04.addParm('delY',dy,PARM_REALS);
  parm04.addParm('delR',dz,PARM_REALS);          
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% BATHYMETRY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%
   
  %%% Flat bottom  
  h = -H*ones(Nx,Ny);
  h(:,1) = 0;
%    if (use_3D)
%      h(:,end) = 0;
%    end
  
  
  %%% Save as a parameter
  writeDataset(h,fullfile(inputpath,'bathyFile.bin'),ieee,prec);
  parm05.addParm('bathyFile','bathyFile.bin',PARM_STR);
 
  
   
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% ICE SHELF DRAFT %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Linear ice topography
  Y = [0:dy:dy*Ny-dy];
  yc = Y+dy/2;
  Z = [0:-dz:-dz*Nr+dz];
  Hmin = -500; % deepest point of cavern
  Hmax = 0; % shallowest point of cavern
  dHdx = (Hmax-Hmin)/0.4;
  icetopo = ones(Nx,Ny).*min(Hmax,Hmin + dHdx*(yc-Y(1))*1e-5);
  pload = -icetopo;
  
  
%   %%% tanh ice topography
%   slope = 0.035;
%   Hshelf = 0;
%   Hmax = 500;
%   hdiff = (Hmax-Hshelf);
%   ywidth = hdiff/(2*slope);
%   yshelf = Lmin/2;
%   offset = 5500;
%   icetopo = ones(Nx,Ny).*hdiff/2.*(tanh((yy-yshelf-offset)/ywidth + 1)) - Hmax/2;
%   pload = -icetopo;
%   %plot(yy,icetopo)


  %writeDataset(icetopo,fullfile(inputpath,'SHELFICEtopoFile'),ieee,prec);
  writeDataset(pload,fullfile(inputpath,'pload'),ieee,prec);

   
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% REFERENCE TEMPERATURE/SALINITY PROFILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
  %%% Quasi-tanh-shaped T/S profiles
%  Zpyc = -150;
%  Wpyc = 50;
%  Smin = 34.3;
%  Smax = 34.7;
%  Tmin = 0.0901-0.0575*Smin; %%% MITgcm surface freezing temperature
%  Tmax= 1;  
%  gam_h = 0.01;
%  tRef = Tmin + (Tmax-Tmin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2));
%  sRef = Smin + (Smax-Smin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2)); 
  
  %%% Quasi-linear near-surface stratification
%   Zpyc = 0;
%   Hpyc = 300;  
%   Smin = 34.3;
%   Smax = 34.7;
%   Tmin = 0.0901-0.0575*Smin; %%% MITgcm surface freezing temperature
%   Tmax= 1;  
%   gam_h = 0.01;
%   tRef = Tmin + (Tmax-Tmin)*0.5*(sqrt((1-zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2)-sqrt((1+zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2))/(1+gam_h)^(1/2);
%   sRef = Smin + (Smax-Smin)*0.5*(sqrt((1-zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2)-sqrt((1+zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2))/(1+gam_h)^(1/2);

%%% Linear T/S profiles           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDED
%  Zpyc = -150;
%  Wpyc = 50;
  Smin = 34.4;
  Smax = 34.6;
  dS = (Smax-Smin)/(Nr-1);
  Tmin = 0.0901-0.0575*Smin; %%% MITgcm surface freezing temperature
  Tmax= 0.5;
  dT = (Tmax-Tmin)/(Nr-1);
  sRef = [Smin : dS : Smax];
  tRef = repmat(Tmax,size(sRef)); 
  %tRef = [Tmin : dT : Tmax];
  
  %%% Plot the reference temperature
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(tRef,zz);
    xlabel('\theta_r_e_f');
    ylabel('z','Rotation',0);
    title('Reference temperature');
  end
  
  %%% Plot the reference salinity
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(sRef,zz);
    xlabel('S_r_e_f');
    ylabel('z','Rotation',0);
    title('Reference salinity');
  end
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFORMATION RADIUS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  %%% Check Brunt-Vaisala frequency using full EOS
  pp = - zz;
%   SA = gsw_SA_from_SP(sRef,pp,-50,-64);  
%   CT = gsw_CT_from_pt(sRef,tRef);
%   [N2 pp_mid] = gsw_Nsquared(SA,CT,pp);
  [N2 pp_mid] = gsw_Nsquared(sRef,tRef,pp);
  dzData = zz(1:end-1)-zz(2:end);

  %%% Calculate internal wave speed and first Rossby radius of deformation
  N = sqrt(max(N2,0)); 
  Cig = sum(N.*dzData);
  Rd = Cig./(pi*abs(f0));

  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    semilogx(N2,-pp_mid/1000);
    xlabel('N^2 (km)');
    ylabel('z (km)','Rotation',0);
    title('Buoyancy frequency');
  end
  
  if (showplots)
    ['First baroclinic Rossby deformation radius: ' num2str(Rd)]
  end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CALCULATE TIME STEP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  
  %%% These estimates are in no way complete, but they give at least some
  %%% idea of the time step needed to keep things stable. In complicated 
  %%% simulations, preliminary tests may be required to estimate the
  %%% parameters used to calculate these time steps.        
  
  %%% Gravity wave CFL

  %%% Upper bound for absolute horizontal fluid velocity (m/s)
  %%% At the moment this is just an estimate
  Umax = 0.5;  
  Wmax = 0.1; %%% NEEDS TO BE INCREASED FOR DEEP PYCNOCLINES
  %%% Max gravity wave speed 
  cmax = max(Cig);
  %%% Max gravity wave speed using total ocean depth
  cgmax = Umax + cmax;
  %%% Gravity wave CFL
  deltaT_gw = min([0.5*dx/cmax,0.5*dy/cmax]);
  %%% Advective CFL
  deltaT_adv= min([0.5*dx/Umax,0.5*dy/Umax]);
  %%% Vertical advective CFL
  deltaT_vadv = min(0.5*dz/Wmax);
  %%% CFL time step based on full gravity wave speed
  deltaT_fgw = min([0.5*dx/cgmax,0.5*dy/cgmax]);
    
  %%% Other stability conditions
  
  %%% Inertial CFL time step (Sf0<=0.5)
  deltaT_itl = 0.5/abs(f0);
  %%% Time step constraint based on horizontal diffusion 
  deltaT_Ah = 0.5*min([dx dy])^2/(4*viscAh);    
  %%% Time step constraint based on vertical diffusion
  deltaT_Ar = 0.5*min(dz)^2 / (4*viscAr);  
  %%% Time step constraint based on biharmonic viscosity 
  deltaT_A4 = 0.5*min([dx dy])^4/(32*viscA4);
  %%% Time step constraint based on horizontal diffusion of temp 
  deltaT_KhT = 0.4*min([dx dy])^2/(4*diffKhT);    
  %%% Time step constraint based on vertical diffusion of temp 
  deltaT_KrT = 0.4*min(dz)^2 / (4*diffKrT);
  
  %%% Time step size  
  deltaT = min([deltaT_fgw deltaT_gw deltaT_adv deltaT_itl deltaT_Ah deltaT_Ar deltaT_KhT deltaT_KrT deltaT_A4]);
  if (nonHydrostatic)
    deltaT = min([deltaT deltaT_vadv]);
  end
  deltaT = round(deltaT)
  nTimeSteps = ceil(simTime/deltaT);
  simTimeAct = nTimeSteps*deltaT;
  
  %%% Write end time time step size  
  parm03.addParm('endTime',nIter0*deltaT+simTimeAct,PARM_INT);
  parm03.addParm('deltaT',deltaT,PARM_REAL); 

  
  
  
  
  
  
  
  
 
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL DATA %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Random noise amplitude
  tNoise = 0.05;  
  sNoise = 0.005;
      
  %%% Align initial temp with background
  hydroTh = ones(Nx,Ny,Nr);
  hydroSa = ones(Nx,Ny,Nr);
  for k=1:1:Nr
    hydroTh(:,:,k) = squeeze(hydroTh(:,:,k))*tRef(k);
    hydroSa(:,:,k) = squeeze(hydroSa(:,:,k))*sRef(k);
  end
  
  %%% Add some random noise
  hydroTh = hydroTh + tNoise*(2*rand(Nx,Ny,Nr)-1);
  hydroSa = hydroSa + sNoise*(2*rand(Nx,Ny,Nr)-1);
  
  %%% Write to data files
  writeDataset(hydroTh,fullfile(inputpath,'hydrogThetaFile.bin'),ieee,prec); 
  parm05.addParm('hydrogThetaFile','hydrogThetaFile.bin',PARM_STR);
  writeDataset(hydroSa,fullfile(inputpath,'hydrogSaltFile.bin'),ieee,prec); 
  parm05.addParm('hydrogSaltFile','hydrogSaltFile.bin',PARM_STR);
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% TRACER DIFFUSION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  diffK4T = diffK4Tgrid * min([dx dy])^4 / (32*deltaT);
  diffK4S = diffK4T;
  parm01.addParm('diffK4T',diffK4T,PARM_REAL); 
  parm01.addParm('diffK4S',diffK4S,PARM_REAL); 
 
  
  
  
    
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SURFACE HEAT/SALT FLUXES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  %%% Fixed salt flux in the BW formation region  
  %salt_flux = zeros(Nx,Ny,nForcingPeriods);    
  %tt = zeros(1,nForcingPeriods);
  %for n=1:nForcingPeriods   
  %  t = (n-1)*externForcingPeriod;
  %  tt(n) = t;
  %  for j=1:Ny      
  %    for i=1:Nx             
  %      if (abs(yy(j)-Ypoly)<Wpoly/2)
          
          %%% Linear decrease of salt flux from salt_amp to zero after Tpoly days
%           salt_amp = -1e-1;
%           Tpoly = 3*t1day;
%           salt_flux(i,j,n) = salt_amp - salt_amp*min(1,t/Tpoly);             
          
          %%% Empirical ice growth equation from Anderson (1969)
%          C1 = 1e4; %%% Constants from Anderson (1969) 
%          C2 = 510;
%          C3 = 6.7/86400;
%          C4 = C2/(2*C1);
%          C5 = C3/C1;
%          Theta = 40; %%% Difference between atmospheric temperature and ocean freezing temperature 
%          rho_i = 920; %%% Density of ice
%          S_i = 5; %%% Salinity of ice
%          e = -C4 + sqrt(C4^2+C5*t*Theta); %%% Ice thickness (m) as a function of time
%          dedt = C3*Theta/(2*C1*e+C2); %%% Ice growth rate (m/s) as a function of time
%          salt_flux(i,j,n) = -dedt*rho_i*(Smin-S_i); %%% Equivalent salt flux (g/m^2/s) as a function of time
          
%        end
%      end
%    end         
%  end  
  
 
  
%  %%% Plot the surface salt flux
%  if (showplots)
%    figure(fignum);
%    fignum = fignum + 1;
%    clf;
%%     contourf(X,Y,squeeze(salt_flux(:,:,1)));
%    colorbar;
%%     plot(yy,squeeze(salt_flux(1,:,1)));    
%%     xlabel('y');
%    plot(tt,squeeze(salt_flux(1,Ny/2,:)));    
%    xlabel('t');
%v    ylabel('salt flux');
%    title('Surface salt flux');
%  end  
  
  %%% Save as parameters  
%  writeDataset(salt_flux,fullfile(inputpath,'saltFluxFile.bin'),ieee,prec);
%  parm05.addParm('saltFluxFile','saltFluxFile.bin',PARM_STR);
  












  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates the 'data' file
  write_data(inputpath,PARM,listterm,realfmt);
  
 
  %%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%
  %%%%%% KPP %%%%%%
  %%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%
  
  if (use_KPP)
      
      %%%%%%%%%%%%%%%%%%%%%%%
      %%%%% KPP SET-UP %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%
      
      
      %%% To store parameter names and values
      kpp_parm01 = parmlist;
      KPP_PARM = {kpp_parm01};
  
      kpp_parm01.addParm('KPPmixingMaps',false,PARM_BOOL);
      kpp_parm01.addParm('KPPwriteState',true,PARM_BOOL);
      % kpp_parm01.addParm('minKPPhbl',,PARM_REAL);
      kpp_parm01.addParm('Ricr',1,PARM_REAL);
      kpp_parm01.addParm('kpp_dumpFreq ',0,PARM_REAL);
      kpp_parm01.addParm('kpp_tavefreq ',0,PARM_REAL);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%% WRITE THE 'data.kpp' FILE %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      %%% Creates the 'data.ggl90' file
      write_data_kpp(inputpath,KPP_PARM,listterm,realfmt);
      
  end
  
  
  %%%%%%%%%%%%%%
  %%%%%%%%%%%%%%
  %%%%% GGL90 %%%%%
  %%%%%%%%%%%%%%
  %%%%%%%%%%%%%%
  
  if (use_GGL90)
      
      %%%%%%%%%%%%%%%%%%%%%%%
      %%%%% GGL SET-UP %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%
      
      
      %%% To store parameter names and values
      ggl90_parm01 = parmlist;
      GGL90_PARM = {ggl90_parm01};
      
      GGL90writeState=true; 
      mxlMaxFlag = 2; % impose mixing length lim
      mxlSurfFlag = true;
      GGL90TKEsurfMin = 1.5e-5; % minimum surface boundary condition
      
      
      ggl90_parm01.addParm('useGGL90',use_GGL90,PARM_BOOL);
      ggl90_parm01.addParm('mxlMaxFlag',mxlMaxFlag,PARM_REAL);
      ggl90_parm01.addParm('mxlSurfFlag',mxlSurfFlag,PARM_REAL);
      ggl90_parm01.addParm('GGL90TKEsurfMin',GGL90TKEsurfMin,PARM_REAL);
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%% WRITE THE 'data.ggl90' FILE %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      %%% Creates the 'data.ggl90' file
      write_data_ggl90(inputpath,GGL90_PARM,listterm,realfmt);
      
  end
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% RBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  
  
  
    
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  rbcs_parm01 = parmlist;
  rbcs_parm02 = parmlist;
  RBCS_PARM = {rbcs_parm01,rbcs_parm02};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  useRBCtemp = false;
  useRBCsalt = false;
  useRBCuVel = false;
  useRBCvVel = false;
  tauRelaxT = t1day;
  tauRelaxS = t1day;
  tauRelaxU = 0.1*t1day;
  tauRelaxV = 0.1*t1day;
  rbcs_parm01.addParm('useRBCtemp',useRBCtemp,PARM_BOOL);
  rbcs_parm01.addParm('useRBCsalt',useRBCsalt,PARM_BOOL);
  rbcs_parm01.addParm('useRBCuVel',useRBCuVel,PARM_BOOL);
  rbcs_parm01.addParm('useRBCvVel',useRBCvVel,PARM_BOOL);
  rbcs_parm01.addParm('tauRelaxT',tauRelaxT,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxS',tauRelaxS,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxU',tauRelaxU,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxV',tauRelaxV,PARM_REAL);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION TEMPERATURE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  %%% Set relaxation temp equal to surface freezing temperature -----> why
  %%% this value again??
  temp_relax = Tmin*ones(Nx,Ny,Nr);  
  salt_relax = Smin*ones(Nx,Ny,Nr); %%% should this be zero? currently useRBCsalt=false

  %%% Save as parameters
  writeDataset(temp_relax,fullfile(inputpath,'sponge_temp.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxTFile','sponge_temp.bin',PARM_STR);  
  writeDataset(salt_relax,fullfile(inputpath,'sponge_salt.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxSFile','sponge_salt.bin',PARM_STR);
  
  %%%%%%%%%%%%%%%%%%%%%  
  %%%%% RBCS MASK %%%%%
  %%%%%%%%%%%%%%%%%%%%%  

  
  %%% Mask is zero everywhere by default, i.e. no relaxation
  slope_msk=zeros(Nx,Ny,Nr); % 2D, for 3D slope_mask=zeros(Nx,Ny,Nr);
  
  for k=1:Nr %%% ----> there must be a better way to do this loop...
    for j=1:Ny %%% till ice shelf
      if zz(k)<icetopo(j)
        if zz(k)<icetopo(j)-dz(1) %%%% ------> should probably have this hFacC*deltaZ
          slope_msk(j,k) = 1;
        end
      end
    end
  end
  
  rbcs_msk = slope_msk;
  
  %%% mask nonzero in sponge layer
  rbcs_msk(Ny-nSponge:Ny,:,:) = 1;
        
  %%% Save as an input parameter
  writeDataset(rbcs_msk,fullfile(inputpath,'rbcs_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskFile(1)','rbcs_mask.bin',PARM_STR);  
  writeDataset(rbcs_msk,fullfile(inputpath,'rbcs_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskFile(2)','rbcs_mask.bin',PARM_STR);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.rbcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
  %%% Creates the 'data.rbcs' file
  write_data_rbcs(inputpath,RBCS_PARM,listterm,realfmt);
   
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% OBCS %%%%%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDED
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%    
 
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% To store parameter names and values
  obcs_parm01 = parmlist;
  obcs_parm02 = parmlist;
  obcs_parm03 = parmlist;
  OBCS_PARM = {obcs_parm01,obcs_parm02,obcs_parm03};
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OBCS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % north
  OB_Jnorth(1:Nx) = -1;
%   OB_Ieast(Ny) = -1;
%   OB_Iwest(Ny) = 1;
  %OB_Jsouth = 1;
  OBCS_monitorFreq = 1800;
  %OBCS_monSelect = 1;
  %obcs_parm01.addParm('OBCSfixTopo',false,PARM_BOOL); % This flag turns off checking and fixing problematic topography
  obcs_parm01.addParm('useOBCSprescribe',true,PARM_BOOL); % not prescribing any velocity defaults to zero flow
  obcs_parm01.addParm('useOBCSbalance',true,PARM_BOOL);
  obcs_parm01.addParm('useOBCSsponge',true,PARM_BOOL); % apply sponge layer
  obcs_parm01.addParm('OBCSbalanceSurf',false,PARM_BOOL);
  obcs_parm01.addParm('OB_Jnorth',OB_Jnorth,PARM_INTS);
 % obcs_parm01.addParm('OB_Jsouth',OB_Jsouth,PARM_INT);
%   obcs_parm01.addParm('OB_Iwest',OB_Iwest,PARM_INT);
%   obcs_parm01.addParm('OB_Ieast',OB_Ieast,PARM_INT);
  obcs_parm01.addParm('OBCS_monitorFreq',OBCS_monitorFreq,PARM_REAL); 
  %obcs_parm01.addParm('OBCS_monSelect',OBCS_monSelect,PARM_INT); % this flag tests balancing the surface flux %I think integer
  obcs_parm01.addParm('OBCS_balanceFacN',1,PARM_INT); % All adjustments are made through the northern boundary

  OBNt = repmat(tRef, [Nx 1]);
  OBNs = repmat(sRef,[Nx 1]);
  writeDataset(OBNt,fullfile(inputpath,'OBNtFile.bin'),ieee,prec);
  obcs_parm01.addParm('OBNtFile','OBNtFile.bin',PARM_STR); 
  writeDataset(OBNs,fullfile(inputpath,'OBNsFile.bin'),ieee,prec);
  obcs_parm01.addParm('OBNsFile','OBNsFile.bin',PARM_STR); 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SPONGE RELAXATION PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~(use_RBCS)
      Vrelaxobcsinner = 86400;
      Vrelaxobcsouter = 8640;
      Urelaxobcsinner = 86400;
      Urelaxobcsouter = 8640;
      spongeThickness = nSponge; % nsponge specified earlier in the meriodional grid resolution.
      obcs_parm03.addParm('Vrelaxobcsinner',Vrelaxobcsinner,PARM_REAL);
      obcs_parm03.addParm('Vrelaxobcsbound',Vrelaxobcsouter,PARM_REAL);
      obcs_parm03.addParm('Urelaxobcsinner',Urelaxobcsinner,PARM_REAL);
      obcs_parm03.addParm('Urelaxobcsbound',Urelaxobcsouter,PARM_REAL);
      obcs_parm03.addParm('spongeThickness',spongeThickness,PARM_REAL); % change to integer if grid point.
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.obcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Creates the 'data.obcs' file
  write_data_obcs(inputpath,OBCS_PARM,listterm,realfmt);
  
  
  

  
  
  
  
    
  
  
  %%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%
  %%%%% SHELF ICE %%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDED
  %%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%
    
  if (use_shelfIce)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SHELF ICE SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% To store parameter names and values
  shelfice_parm01 = parmlist;
  SHELFICE_PARM = {shelfice_parm01};
 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% ICE SHELF PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  shelfice_parm01.addParm('SHELFICEboundaryLayer',true,PARM_BOOL);
  shelfice_parm01.addParm('SHELFICEMassStepping',true,PARM_BOOL);
  shelfice_parm01.addParm('SHELFICEconserve',true,PARM_BOOL);
  %shelfice_parm01.addParm('SHELFICEadvDiffHeatFlux',true,PARM_BOOL);
  shelfice_parm01.addParm('SHELFICEuseGammaFrict',true,PARM_BOOL);
  shelfice_parm01.addParm('SHELFICEselectDragQuadr',1,PARM_INT);
  shelfice_parm01.addParm('rhoShelfIce',917,PARM_REAL);
  shelfice_parm01.addParm('SHELFICEwriteState',true,PARM_BOOL);
  shelfice_parm01.addParm('SHELFICEheatTransCoeff',0,PARM_REAL);
  shelfice_parm01.addParm('SHELFICEHeatCapacity_Cp',2000,PARM_REAL);
  % shelfice_parm01.addParm('SHI_withBL_uStarTopDz',true,PARM_BOOL); % ----> Add? might help with spurious variability. Used for shelfice remeshing but is it applied here too?

  
  
  %%% Writing icetopo.obcs file
   if (use_OBCS)
     yShelf = 400; % grid point where sloping shelf becomes horizontal, should think of a different way to define this with the variable gird calculation
     icetopo_obcs = icetopo;
     icetopo_obcs(yShelf:end,:) = 0;
     writeDataset(icetopo_obcs,fullfile(inputpath,'icetopo'),ieee,prec)
     shelfice_parm01.addParm('SHELFICEtopoFile','icetopo',PARM_STR);
  else
    writeDataset(icetopo,fullfile(inputpath,'icetopo'),ieee,prec); %%%%%%%%%%% How do I save this file correctly?
    shelfice_parm01.addParm('SHELFICEtopoFile','icetopo',PARM_STR);
  end

  
  %%% Creating phi0surf 
  
%   offsetZ = 200; % thickness of ice shelf extension (not sloping shelf). This thickness is not included in the domain but needs to be adjusted to the pressure in the calculations below
%   zgp1 = [0,cumsum(dz)];
%   zc = 0.5*(zgp1(1:end-1)+zgp1(2:end));
%   zg = zgp1(1:end-1);
%   dzm = abs([zg(1)-zc(1) 0.5*diff(zc)]);
%   dzp = abs([0.5*diff(zc) zc(end)-zg(end)]);
%   p = abs(zc)*g*rho0*1e-4;
%   dp = p;
%   kp = 0;
%   talpha = 2e-4;
%   sbeta = 7.4e-4; 
%   
%   while rms(dp)>1e-13
%     phiHydF(k) = 0;
%     p0 = p;
%     kp = kp+1;
%     for k=1:Nr
%       switch eos
%         case 'linear'
%           drho = rho0*(1-talpha*(t(k)-tref(k))+sbeta*(s(k)-sref(k))) - rho0;
%         case 'JMD95Z'
%           drho = densjmd95(sRef(k),tRef(k),p(k)+offsetZ) - rho0; %%% offset still needed?
%         case 'MDJWF'
%           drho = densmdjwf(s(k),t(k),p(k)+offsetZ) - rho0; %%% offset still needed?
%         otherwise
%         error(sprintf('unknown EOS: %s',eos))
%       end
%       phiHydC(k) = phiHydF(k) + dzm(k)*g*drho/rho0;
%       phiHydF(k+1) = phiHydC(k) + dzp(k)*g*drho/rho0;
%     end
%     switch eos
%       case 'MDJWF'
%         p = (g*rhoConst*abs(zc) + phiHydC*rhoConst) / g / rhoConst;
%     end
%     dp = p - p0; 
%   end
  
  
  
  
  
  %%% Compute hFac following the MITgcm source code
  hFacC = zeros(Nx,Ny,Nr);
  hFacC(:,1,101:Nr) = 1;
  for k=1:Nr
    hFacMnSz = max( hFacMin, min(hFacMinDr/dz(k),1) );
    for j=2:Ny
      for i =1:Nx
        % Non-dimensional distance between grid bound and domain
        % lower_R bound
        hFacCtmp = (min(zzf(k),icetopo(i,j))-max(zzf(k+1),h(i,j))) / dz(k);
        % Select between closed, open or partial ( 0, 1, 0-1)
        hFacCtmp = min( max(hFacCtmp,0), 1);
        % Impose minimal fraction and/or size
        if (hFacCtmp < hFacMnSz)
          if ( hFacCtmp < hFacMnSz*0.5)
            hFacC(i,j,k) = 0;
          else
            hFacC(i,j,k) = hFacMnSz;
          end
        else
          hFacC(i,j,k) = hFacCtmp;
        end
      end
    end
  end
  
  %Ice shelf load anomaly 
  salt_ref = 34.4; %% should I be using these values or Smin and Tmin?
  temp_ref = -1.9; 
  
  MIce_Anom = zeros(Nx,Ny);
  for i=1:Nx
      for j=1:yShelf
          MIce_Anom(i,j) = 0;
          for k = 1:Nr
              if ( zz(k) <icetopo(i,j) )
                  break;
              end
              idx(k) = j;
              Pressure = -rho0*g*zz(k);
              rhoShelfIce= densjmd95(salt_ref,temp_ref,Pressure/Pa1dbar);
              MIce_Anom(i,j) = MIce_Anom(i,j) + (g*(rhoShelfIce-rho0)*dz(k)*(1-hFacC(i,j,k)));
          end
      end
  end
 
  %plot(MIce_Anom)
  
  
  writeDataset(MIce_Anom,fullfile(inputpath,'SHELFICEloadAnomalyFile'),ieee,prec);
  shelfice_parm01.addParm('SHELFICEloadAnomalyFile','SHELFICEloadAnomalyFile',PARM_STR);
%   Calculating phi0surf  

%   msk = sum(hFacC,3); msk(msk>0)=1;
%   phi0surf = zeros(Nx,Ny);
%   for ix = 1:Nx
%     for iy = 1:Ny
%       k = max(find(abs(zg)<abs(icetopo(ix,iy))));
%       if isempty(k)
%         k=0;
%       end
%       if k>0
%         kp1 = min(k+1,Nr);
%         drloc = 1 - hFacC(ix,iy,k);
%         dphi = phiHydF(kp1)-phiHydF(k);
%         phi0surf(ix,iy) = (phiHydF(k)+drloc*dphi)*rho0*msk(ix,iy);
%       end
%     end
%   end
  
%  writeDataset(phi0surf,fullfile(inputpath,'phi0surf'),ieee,prec); %----> don't think this is saved right 
  %rbcs_parm01.addParm('relaxMaskFile(1)','rbcs_mask.bin',PARM_STR);  
  
%   %%% Creating iceShelf_Mass file
%   rSurf0 = rdmds('rSurfC');
%   fid = fopen(['phi0surf'],'r','b');
%   pLoadAn = fread(fid,[Nx Ny],'real*8'); fclose(fid); 
%   mIce = pLoadAn'/g-rSurf0*rho0;
%   if (use_OBCS)
%     mIce(yShelf:end,:) = 0;
%     writeDataset(mIce,fullfile(inputpath,'iceShelf_Mass'),ieee,prec);
%     shelfice_parm01.addParm('SHELFICEMassFile','iceShelf_Mass',PARM_STR);
%   else
%     writeDataset(mIce,fullfile(inputpath,'iceShelf_Mass'),ieee,prec);
%     shelfice_parm01.addParm('SHELFICEMassFile','iceShelf_Mass',PARM_STR);
%   end
%   
%   %%% Creating iceShelf_MassTend file  
%   mIceTend = zeros(size(mIce));
%   mIceTend(1:yShelf,:) = mIce*1e-5/deltaT;
%   writeDataset(mIceTend,fullfile(inputpath,'iceShelf_MassTend'),ieee,prec); 
%   shelfice_parm01.addParm('SHELFICEMassDynTendFile','iceShelf_MassTend',PARM_STR);
% 
%   %%% Creating iceMask file
%   mskIce = ones(size(mIce));
%   mskIce(yShelf:end,:) = 2;
%   writeDataset(mskIce,fullfile(inputpath,'iceMask'),ieee,prec)
%   %shelfice_parm01.addParm('SHELFICEtopoFile','iceMask',PARM_STR); 
  
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.shelfice' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Creates the 'data.shelfice' file
  write_data_shelfice(inputpath,SHELFICE_PARM,listterm,realfmt);
  
  end
  
  
  
  
  
  
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   %%%%%%%%% SEA ICE $ %%%%%%%%%%%
% %   
% %   
% %   % to store parameter names and values
%   seaice_parm01 = parmlist;
%   SEAICE_PARM = {seaice_parm01};
% % %   
% % %   
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%% SEA ICE  %%%%%%%%%%%%
%     %%%%%%%% PARAMETERS %%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Oringinal albedos from llc_1/48th  other values from llc_2160
% % 
%   SEAICEwriteState   = true;
%   SEAICEuseDYNAMICS  = true;
%   SEAICE_multDim     = 7;
% %   SEAICE_dryIceAlb   = 0.8783;
%   SEAICE_dryIceAlb   = 0.8509;
% %   SEAICE_wetIceAlb   = 0.7869;
%   SEAICE_wetIceAlb   = 0.7284;
% %   SEAICE_drySnowAlb  = 0.9482;
%   SEAICE_drySnowAlb  = 0.7754;
% %   SEAICE_wetSnowAlb  = 0.8216;
%   SEAICE_wetSnowAlb  = 0.7753;
%   SEAICE_waterDrag   = 5.5399;
%   SEAICE_drag        = 0.002;
% %   HO                 = 0.1;
%   HO                 = .05;
% 
% 
%   SEAICE_no_slip          = false;
% %   SEAICE_no_slip          = true;
% 
%   SEAICEadvScheme         = 7;
% %   SEAICEadvScheme         = 33;
% 
% 
%   %%%SOSEdoesn't have a seaice dataset for salinity, they used this value
%   %%%in their estimate
%   
%   LSR_ERROR               = 1.0e-4;
% %   LSR_ERROR               = 2.0e-4; 
%   MIN_ATEMP               = -40;
%   MIN_TICE                = -40;
%   SEAICE_area_reg         = 0.15;
%   SEAICE_hice_reg         = 0.1;
%   IMAX_TICE               = 6;
%   SEAICE_EPS		      = 1.0e-8;
% %   SEAICE_EPS              = 2.0e-9;
%   SEAICE_doOpenWaterMelt  = true;
%   SEAICE_areaLossFormula  = 1;
%   SEAICE_wetAlbTemp       = 0.0;
%   SEAICE_saltFrac         = 0.3;
% %   SEAICE_frazilFrac       = 0.003;
%  SEAICE_frazilFrac       = 0.01;
% %   SEAICE_frazilFrac       = 1.0;
% 
%   
%   seaice_parm01.addParm('LSR_ERROR',LSR_ERROR,PARM_REAL);
%   seaice_parm01.addParm('SEAICEwriteState',SEAICEwriteState,PARM_BOOL);
%   seaice_parm01.addParm('SEAICEuseDYNAMICS',SEAICEuseDYNAMICS,PARM_BOOL);
%   seaice_parm01.addParm('SEAICE_multDim',SEAICE_multDim,PARM_INT);
%   seaice_parm01.addParm('SEAICE_dryIceAlb',SEAICE_dryIceAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_wetIceAlb',SEAICE_wetIceAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_drySnowAlb',SEAICE_drySnowAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_wetSnowAlb',SEAICE_wetSnowAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_waterDrag',SEAICE_waterDrag,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_drag',SEAICE_drag,PARM_REAL);
%   seaice_parm01.addParm('HO',HO,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_dryIceAlb_south',SEAICE_dryIceAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_wetIceAlb_south',SEAICE_wetIceAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_drySnowAlb_south',SEAICE_drySnowAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_wetSnowAlb_south',SEAICE_wetSnowAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_waterDrag_south',SEAICE_waterDrag_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_drag_south',SEAICE_drag_south,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_no_slip',SEAICE_no_slip,PARM_BOOL);
% %   seaice_parm01.addParm('SEAICE_salinity',SEAICE_salinity,PARM_REAL);
%   seaice_parm01.addParm('SEAICEadvScheme',SEAICEadvScheme,PARM_INT);
%   seaice_parm01.addParm('MIN_ATEMP',MIN_ATEMP,PARM_REAL);
%   seaice_parm01.addParm('MIN_TICE',MIN_TICE,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_area_reg',SEAICE_area_reg,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_hice_reg',SEAICE_hice_reg,PARM_REAL);
%   seaice_parm01.addParm('IMAX_TICE',IMAX_TICE,PARM_INT);
%   seaice_parm01.addParm('SEAICE_EPS',SEAICE_EPS,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_doOpenWaterMelt',SEAICE_doOpenWaterMelt,PARM_BOOL);
%   seaice_parm01.addParm('SEAICE_areaLossFormula',SEAICE_areaLossFormula,PARM_INT);
%   seaice_parm01.addParm('SEAICE_wetAlbTemp',SEAICE_wetAlbTemp,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_saltFrac',SEAICE_saltFrac,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_frazilFrac',SEAICE_frazilFrac,PARM_REAL);
% 
%   
%   seaice_parm01.addParm('HeffFile',HeffFile,PARM_STR);
%   seaice_parm01.addParm('AreaFile',AreaFile,PARM_STR);
%   seaice_parm01.addParm('HsnowFile',HsnowFile,PARM_STR);
%   seaice_parm01.addParm('HsaltFile',HsaltFile,PARM_STR);
%   seaice_parm01.addParm('uIceFile',uIceFile,PARM_STR);
%   seaice_parm01.addParm('vIceFile',vIceFile,PARM_STR);
%   
%   
%      
% %   
% %   
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   %%%%% WRITE THE 'data.seaice' FILE %%%
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   
% %   
%   write_data_seaice(inputpath,SEAICE_PARM,listterm,realfmt);  
%   
%   
%   
%   
%   
%   
%   
%   
%   
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%
%     %%%%%%EXF PKG%%%%%
%     %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%
%     
%     
%   %%% To store parameter names and values
%     
%      EXF_NML_01 = parmlist;
%      EXF_NML_02 = parmlist;
%      EXF_NML_03 = parmlist;
%      EXF_NML_04 = parmlist;
%      EXF_NML_OBCS = parmlist;
% 
%      EXF_PARM = {EXF_NML_01,EXF_NML_02,EXF_NML_03,EXF_NML_04,EXF_NML_OBCS};  
%     
%     
% 
% 
% % %     EXF_NML_01
% 
%   	exf_albedo        = 0.15;
%  	exf_scal_BulkCdn  = 1.015;
%  	exf_iprec         = 64;  
%  	useExfYearlyFields= false;
%  	useExfCheckRange  = false;
%  	useRelativeWind   = true;
% %  	useRelativeWind   = false;
%     repeatPeriod      = 283996800.0;
% %     exf_offset_atemp =  273.16;
%     
%     
% %%%runoff from ERA is in hours, need to convert to seconds
% %     exf_inscal_runoff = 1.14e-04;
%     
%     
%     apressurefile     = 'apressurefile.bin';
%     atempfile         = 'atempfile.bin';
%     aqhfile           = 'aqhfile.bin';
%     uwindfile         = 'uwindfile.bin';
%     vwindfile         = 'vwindfile.bin';
%     precipfile        = 'precipfile.bin';
%     swdownfile        = 'swdownfile.bin';
%     lwdownfile        = 'lwdownfile.bin';
% %     runofffile        = 'runofffile.bin';
%     
%     
% % "*period=-12" specifies monthly-mean forcing
%    apressurestartdate1 = 20070101;
%    apressurestartdate2 = 000000;
% %    apressureperiod     = 86400.0;
%    apressureperiod     = 10800.0;
%     
%     aqhstartdate1 = 20070101;
%     aqhstartdate2 = 000000;
% %     aqhperiod = 86400.0;
%     aqhperiod           = 10800.0;
% 
%     atempstartdate1 = 20070101;
%     atempstartdate2 = 000000;
% %     atempperiod = 86400.0;
%     atempperiod         = 10800.0;
%  
% 
%     uwindstartdate1 = 20070101;
%     uwindstartdate2 = 000000;
% %    uwindperiod = 86400.0;
%     uwindperiod         = 10800.0;
%  
%     vwindstartdate1 = 20070101;
%     vwindstartdate2 = 000000;
% %     vwindperiod = 86400.0;
%     vwindperiod         = 10800.0; 
%  
%     precipstartdate1 = 20070101;
%     precipstartdate2 = 000000;
% %     precipperiod = 86400.0;
%     precipperiod        = 10800.0; 
%  
% 
%     swdownstartdate1 = 20070101;
%     swdownstartdate2 = 000000;
% %     swdownperiod = 86400.0;
%     swdownperiod        = 10800.0;
% % 
% 
%     lwdownstartdate1 = 20070101;
%     lwdownstartdate2 = 000000;
% %     lwdownperiod = 86400.0;
%     lwdownperiod        = 10800.0;
% 
% 
% %     runoffstartdate1 = 20070101;
% %     runoffstartdate2 = 000000;
% %     runoffperiod = 2592000.0;
% %    runoffperiod        = 10800.0;
%  
%    precip_lon0 = xmin;
%    precip_lon_inc = dmxg(1);
%    precip_lat0 = ymin;
%    precip_lat_inc = dmyg;
%    precip_nlon = Nx;
%    precip_nlat = Ny;
%    
%    atemp_lon0 = xmin;
%    atemp_lon_inc = dmxg(1);
%    atemp_lat0 = ymin;
%    atemp_lat_inc = dmyg;
%    atemp_nlon = Nx;
%    atemp_nlat = Ny;
%    
%    apressure_lon0 = xmin;
%    apressure_lon_inc = dmxg(1);
%    apressure_lat0 = ymin;
%    apressure_lat_inc = dmyg;
%    apressure_nlon = Nx;
%    apressure_nlat = Ny;
%     
%    aqh_lon0 = xmin;
%    aqh_lon_inc = dmxg(1);
%    aqh_lat0 = ymin;
%    aqh_lat_inc = dmyg;
%    aqh_nlon = Nx;
%    aqh_nlat = Ny;
%    
%    uwind_lon0 = xmin;
%    uwind_lon_inc = dmxg(1);
%    uwind_lat0 = ymin;
%    uwind_lat_inc = dmyg;
%    uwind_nlon = Nx;
%    uwind_nlat = Ny;
%    
%    vwind_lon0 = xmin;
%    vwind_lon_inc = dmxg(1);
%    vwind_lat0 = ymin;
%    vwind_lat_inc = dmyg;
%    vwind_nlon = Nx;
%    vwind_nlat = Ny;
%    
%    swdown_lon0 = xmin;
%    swdown_lon_inc = dmxg(1);
%    swdown_lat0 = ymin;
%    swdown_lat_inc = dmyg;
%    swdown_nlon = Nx;
%    swdown_nlat = Ny;
%    
%    lwdown_lon0 = xmin;
%    lwdown_lon_inc = dmxg(1);
%    lwdown_lat0 = ymin;
%    lwdown_lat_inc = dmyg;
%    lwdown_nlon = Nx;
%    lwdown_nlat = Ny;
% 
% %    runoff_lon0 = xmin;
% %    runoff_lon_inc = dmxg(1);
% %    runoff_lat0 = ymin;
% %    runoff_lat_inc = dmyg;
% %    runoff_nlon = Nx;
% %    runoff_nlat = Ny;
%   
%  

%      
% %%%%%%%%%%%%%%%%%%%% ADD PARAMETERS
% 
%   EXF_NML_01.addParm('exf_albedo',exf_albedo,PARM_INT);
%   EXF_NML_01.addParm('exf_scal_BulkCdn',exf_scal_BulkCdn,PARM_REAL);
%   EXF_NML_01.addParm('exf_iprec',exf_iprec,PARM_INT);
%   EXF_NML_01.addParm('useExfYearlyFields',useExfYearlyFields,PARM_BOOL);
%   EXF_NML_01.addParm('useExfCheckRange',useExfCheckRange,PARM_BOOL);
%   EXF_NML_01.addParm('useRelativeWind',useRelativeWind,PARM_BOOL);
%   EXF_NML_01.addParm('repeatPeriod',repeatPeriod,PARM_REAL);
% %   EXF_NML_03.addParm('exf_offset_atemp',exf_offset_atemp,PARM_REAL);
% %   EXF_NML_03.addParm('exf_inscal_runoff',exf_inscal_runoff,PARM_REAL);
%   
%   
%   EXF_NML_02.addParm('apressurefile',apressurefile,PARM_STR);
%   EXF_NML_02.addParm('atempfile',atempfile,PARM_STR);
%   EXF_NML_02.addParm('aqhfile',aqhfile,PARM_STR);
%   EXF_NML_02.addParm('uwindfile',uwindfile,PARM_STR);
%   EXF_NML_02.addParm('vwindfile',vwindfile,PARM_STR);
%   EXF_NML_02.addParm('precipfile',precipfile,PARM_STR);
%   EXF_NML_02.addParm('swdownfile',swdownfile,PARM_STR);
%   EXF_NML_02.addParm('lwdownfile',lwdownfile,PARM_STR);
% %   EXF_NML_02.addParm('runofffile',runofffile,PARM_STR);
% 
% 
%   EXF_NML_02.addParm('apressurestartdate1',apressurestartdate1,PARM_INT);
%   EXF_NML_02.addParm('apressurestartdate2',apressurestartdate2,PARM_INT);
%   EXF_NML_02.addParm('apressureperiod',apressureperiod,PARM_REAL);
% 
%   EXF_NML_02.addParm('atempstartdate1',atempstartdate1,PARM_INT);
%   EXF_NML_02.addParm('atempstartdate2',atempstartdate2,PARM_INT);
%   EXF_NML_02.addParm('atempperiod',atempperiod,PARM_REAL);
%   
%   EXF_NML_02.addParm('aqhstartdate1',aqhstartdate1,PARM_INT);
%   EXF_NML_02.addParm('aqhstartdate2',aqhstartdate2,PARM_INT);
%   EXF_NML_02.addParm('aqhperiod',aqhperiod,PARM_REAL);
%   
%   EXF_NML_02.addParm('uwindstartdate1',uwindstartdate1,PARM_INT);
%   EXF_NML_02.addParm('uwindstartdate2',uwindstartdate2,PARM_INT);
%   EXF_NML_02.addParm('uwindperiod',uwindperiod,PARM_REAL);
%   
%   EXF_NML_02.addParm('vwindperiod',vwindperiod,PARM_REAL);
%   EXF_NML_02.addParm('vwindstartdate1',vwindstartdate1,PARM_INT);
%   EXF_NML_02.addParm('vwindstartdate2',vwindstartdate2,PARM_INT);
%   
%   EXF_NML_02.addParm('precipperiod',precipperiod,PARM_REAL);
%   EXF_NML_02.addParm('precipstartdate1',precipstartdate1,PARM_INT);
%   EXF_NML_02.addParm('precipstartdate2',precipstartdate2,PARM_INT);
%   
%   EXF_NML_02.addParm('swdownperiod',swdownperiod,PARM_REAL);
%   EXF_NML_02.addParm('swdownstartdate1',swdownstartdate1,PARM_INT);
%   EXF_NML_02.addParm('swdownstartdate2',swdownstartdate2,PARM_INT);
%   
%   EXF_NML_02.addParm('lwdownperiod',lwdownperiod,PARM_REAL);
%   EXF_NML_02.addParm('lwdownstartdate1',lwdownstartdate1,PARM_INT);
%   EXF_NML_02.addParm('lwdownstartdate2',lwdownstartdate2,PARM_INT);
%   
% %   EXF_NML_02.addParm('runoffperiod',runoffperiod,PARM_REAL);
% 
%   EXF_NML_04.addParm('precip_lon0',precip_lon0,PARM_REAL);
%   EXF_NML_04.addParm('precip_lon_inc',precip_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('precip_lat0',precip_lat0,PARM_REAL);
%   EXF_NML_04.addParm('precip_lat_inc',precip_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('precip_nlon',precip_nlon,PARM_INT);
%   EXF_NML_04.addParm('precip_nlat',precip_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('atemp_lon0',atemp_lon0,PARM_REAL);
%   EXF_NML_04.addParm('atemp_lon_inc',atemp_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('atemp_lat0',atemp_lat0,PARM_REAL);
%   EXF_NML_04.addParm('atemp_lat_inc',atemp_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('atemp_nlon',atemp_nlon,PARM_INT);
%   EXF_NML_04.addParm('atemp_nlat',atemp_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('apressure_lon0',apressure_lon0,PARM_REAL);
%   EXF_NML_04.addParm('apressure_lon_inc',apressure_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('apressure_lat0',apressure_lat0,PARM_REAL);
%   EXF_NML_04.addParm('apressure_lat_inc',apressure_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('apressure_nlon',apressure_nlon,PARM_INT);
%   EXF_NML_04.addParm('apressure_nlat',apressure_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('aqh_lon0',aqh_lon0,PARM_REAL);
%   EXF_NML_04.addParm('aqh_lon_inc',aqh_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('aqh_lat0',aqh_lat0,PARM_REAL);
%   EXF_NML_04.addParm('aqh_lat_inc',aqh_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('aqh_nlon',aqh_nlon,PARM_INT);
%   EXF_NML_04.addParm('aqh_nlat',aqh_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('uwind_lon0',uwind_lon0,PARM_REAL);
%   EXF_NML_04.addParm('uwind_lon_inc',uwind_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('uwind_lat0',uwind_lat0,PARM_REAL);
%   EXF_NML_04.addParm('uwind_lat_inc',uwind_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('uwind_nlon',uwind_nlon,PARM_INT);
%   EXF_NML_04.addParm('uwind_nlat',uwind_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('vwind_lon0',vwind_lon0,PARM_REAL);
%   EXF_NML_04.addParm('vwind_lon_inc',vwind_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('vwind_lat0',vwind_lat0,PARM_REAL);
%   EXF_NML_04.addParm('vwind_lat_inc',vwind_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('vwind_nlon',vwind_nlon,PARM_INT);
%   EXF_NML_04.addParm('vwind_nlat',vwind_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('swdown_lon0',swdown_lon0,PARM_REAL);
%   EXF_NML_04.addParm('swdown_lon_inc',swdown_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('swdown_lat0',swdown_lat0,PARM_REAL);
%   EXF_NML_04.addParm('swdown_lat_inc',swdown_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('swdown_nlon',swdown_nlon,PARM_INT);
%   EXF_NML_04.addParm('swdown_nlat',swdown_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('lwdown_lon0',lwdown_lon0,PARM_REAL);
%   EXF_NML_04.addParm('lwdown_lon_inc',lwdown_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('lwdown_lat0',lwdown_lat0,PARM_REAL);
%   EXF_NML_04.addParm('lwdown_lat_inc',lwdown_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('lwdown_nlon',lwdown_nlon,PARM_INT);
%   EXF_NML_04.addParm('lwdown_nlat',lwdown_nlat,PARM_INT);
% 
% %   EXF_NML_04.addParm('runoff_lon0',runoff_lon0,PARM_REAL);
% %   EXF_NML_04.addParm('runoff_lon_inc',runoff_lon_inc,PARM_REAL);
% %   EXF_NML_04.addParm('runoff_lat0',runoff_lat0,PARM_REAL);
% %   EXF_NML_04.addParm('runoff_lat_inc',runoff_lat_inc,PARM_REALS);
% %   EXF_NML_04.addParm('runoff_nlon',runoff_nlon,PARM_INT);
% %   EXF_NML_04.addParm('runoff_nlat',runoff_nlat,PARM_INT);
% 
% 
%   %%z% Create the data.exf file
%   write_data_exf(inputpath,EXF_PARM,listterm,realfmt);






  
  
  
  
  
  
 
 
  
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
    
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  diag_parm01 = parmlist;
  diag_parm02 = parmlist;
  DIAG_PARM = {diag_parm01,diag_parm02};
  diag_matlab_parm01 = parmlist;
  DIAG_MATLAB_PARM = {diag_matlab_parm01}; %%% Matlab parameters need to be different
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%% Stores total number of diagnostic quantities
  ndiags = 0;
  
  diag_fields_avg = ...
  {...

%     'UV_VEL_Z','WU_VEL','WV_VEL'... %%% Momentum fluxes
%     'UVELSLT','VVELSLT','WVELSLT', ... %%% Salt fluxes

%     'PHIHYD', ... %%% Pressure
%     'LaUH1TH','LaVH1TH','LaHw1TH','LaHs1TH' ... %%% LAYERS fluxes
%     'ADVr_TH','ADVx_TH','ADVy_TH', ... %%%%% Advective fluxes of H.B.
%     'KPPg_TH', ... %%%% KPP non-local flux of P.T.
%     'DFrE_TH','DFrI_TH','DFxE_TH','DFyE_TH', ... %%%% Diffusive Fluxes
%     'oceQsw','TOTTTEND','WTHMASS','TFLUX'...
%      'UVELSQ','VVELSQ','WVELSQ',... %%% For Kinetic Energy
    %%'ETAN',...
%    'SHIgammT','SHIgammS','SHIuStar','SHI_mass', ... %%%%% Ice shelf melt
%      'SIarea','SIheff','SIhsnow','SIhsalt','SIuice','SIvice'...%%% Sea ice state
%      'UVEL','VVEL','WVEL'... %%% Velocities
%      'THETA', ... %%% Temperature
%      'SALT', ... %%% Salinity
%      'SHIfwFlx'...%%'SHIhtFlx','SHIForcT','SHIForcS', ...
%      'UVELTH','VVELTH','WVELTH', ... %%% Temperature fluxes
  };
  
  numdiags_avg = length(diag_fields_avg);  
  diag_freq_avg = 1*t1day;
  diag_phase_avg = 0;    
     
  diag_parm01.addParm('diag_mnc',true,PARM_BOOL);  
  for n=1:numdiags_avg    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields{1,',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames{',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL);  
    
  end
  
  diag_fields_inst = ...
  {...
    'UVEL','VVEL','WVEL', ...%%% velocities
    'THETA', ... %%% Temperature
    'SALT', ... %%% Salinity
    }; 
  
%   diag_fields_inst1 = ...
%   {...  
%     'ETAN', ...  %%% SSH
%     'ETANSQ', ... %%% square of surface height anomaly
%     'DETADT2', ... %%% square of surface height anomaly tendency
%     'SHIuStar', ... %%% frictional velocity?
%     'SHIgammT','SHIgammS',... %%%  shelfice heat and salt transfer coefficient
%     'TFLUX', 'SFLUX', ... %%% total heat and salt fluxes (match heat/salt-content variations), >0 incearses theta/salt 
%     'SHI_TauX','SHI_TauY' ... %%% ice shelf stresses
%     };
% 
%   diag_fields_inst2 = ...
%   {...
%     'UVEL','VVEL','WVEL', ...%%% velocities
%     'THETA', ... %%% Temperature
%     'SALT', ... %%% Salinity
%     'PHIHYD', ... %%% Hydrostatic pressure
%     'VVELMASS','UVELMASS', ... %%% Zonal and meridional mass-weighted comp of velocity
%     'WVELSQ' ... %%% square of vertical comp of velocity
%      %'KPPhbl' ... %%% BL depth
%     }; 
%   
%   diag_fields_inst3 = ...
%   {...
%     'GGL90TKE','GGL90Lmx','GGL90Prl','GGL90Kr' ...%%% 
%     }; 
% 
%   diag_fields_inst4 = ...
%   {...
%     'ADVx_TH','ADVy_TH','ADVr_TH', ... %%%
%     'DFxE_TH','DFyE_TH','DFrE_TH', ...%%%
%     'DFrI_TH' ...
%     }; 

%   numdiags_inst = length(diag_fields_inst1) + length(diag_fields_inst2) + length(diag_fields_inst3) + length(diag_fields_inst4);
  numdiags_inst = length(diag_fields_inst);
  diag_freq_inst = 1*t1day;
  
%   diag_freq_inst1= 1.*t1day;
%   n1 = length(diag_fields_inst1);
%   diag_freq_inst2 = 60.*t1day; %%% ??????
%   n2 = length(diag_fields_inst2) + n1;
%   diag_freq_inst3 = 20.*t1day;
%   n3 = length(diag_fields_inst3) + n2;
%   diag_freq_inst4 = 20.*t1day;
%   n4 = length(diag_fields_inst4) + n3;
% 
%   for n=1:numdiags_inst
%       if n<=n1
%           diag_fields_inst(n) = diag_fields_inst1(n);
%       elseif n>n1 && n<=n2
%           diag_fields_inst(n) = diag_fields_inst2(n-n1);
%       elseif n>n2 && n<=n3
%           diag_fields_inst(n) = diag_fields_inst3(n-n2);
%       elseif n>n3 && n<=n4
%           diag_fields_inst(n) = diag_fields_inst4(n-n3);
%       end
%   end
%   diag_freq_inst(1:n1) = diag_freq_inst1;
%   diag_freq_inst((1+n1) : n2) = diag_freq_inst2;
%   diag_freq_inst((1+n2) : n3) = diag_freq_inst3;
%   diag_freq_inst((1+n3) : n4) = diag_freq_inst4;

  diag_phase_inst = 0;
  
  for n=1:numdiags_inst    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);
    diag_parm01.addParm(['timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL);  
    
  end
  
%   nSetRegMskFile = 1;
%   set_regMask(1:2) = 1;
%   val_regMask(1) = 1; val_regMask(2) = 2;
%   ndiag2 = length(set_regMask);
%   
%   diag_parm02.addParm('diagSt_regMaskFile','iceMask',PARM_STR); 
%   
%   for n = 1:ndiag2
%     diag_parm02.addParm(['set_regMask(',num2str(ndiags),')'],set_regMask(n),PARM_REAL); 
%     diag_parm02.addParm(['val_regMask(',num2str(ndiags),')'],val_regMask(n),PARM_REAL); 
%   end
  
  
  %%% Create the data.diagnostics file
  write_data_diagnostics(inputpath,DIAG_PARM,listterm,realfmt);
  
  %%% Create the DIAGNOSTICS_SIZE.h file
  createDIAGSIZEh(codepath,ndiags,Nr);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%% PACKAGES %%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  
  packages = parmlist;
  PACKAGE_PARM = {packages};  
  
  packages.addParm('useDiagnostics',true,PARM_BOOL);    
  packages.addParm('useKPP',use_KPP,PARM_BOOL);
  packages.addParm('useRBCS',use_RBCS,PARM_BOOL);      
  packages.addParm('useSEAICE',use_seaIce,PARM_BOOL);  
  packages.addParm('useEXF',use_seaIce,PARM_BOOL);
  packages.addParm('useShelfIce',use_shelfIce,PARM_BOOL); 
  packages.addParm('useOBCS',use_OBCS,PARM_BOOL);
  packages.addParm('useGGL90',use_GGL90,PARM_BOOL);
  
  
  %%% Create the data.pkg file
  write_data_pkg(inputpath,PACKAGE_PARM,listterm,realfmt);
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE PARAMETERS TO A MATLAB FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates a matlab file defining all input parameters
  ALL_PARMS = [PARM PACKAGE_PARM DIAG_MATLAB_PARM];
  if (~nonHydrostatic)
    ALL_PARMS = [ALL_PARMS KPP_PARM];
  end
  if (use_shelfIce)
    ALL_PARMS = [ALL_PARMS SHELFICE_PARM];
  else 
    ALL_PARMS = [ALL_PARMS];
  end    
  write_matlab_params(inputpath,ALL_PARMS,realfmt);
  
end
