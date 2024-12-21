%%%
%%% diffusive_instability.m
%%%
%%% Solves the linear eigenvalue problem for 1D diffusive pycnocline
%%% instability. The structure of the solution vector is
%%% [u',v',w',b',phi'].
%%% 
%%% INPUT: Parameter structure 'params'
%%% Required fields:
%%%   M       Number of z-gridpoints 
%%%   Lz      Domain width in z (m)
%%%   Bz      Buoyancy frequency as a function of z (M-length vector, s^-2)
%%%   l       Along-slope wavenumber (1/m)
%%%
%%% Optional fields:
%%%   f               Coriolis parameter (1/s), default=-10^-4
%%%   s               Ice slope, default=10^-2
%%%
%%% OUTPUT: Result structure 'result'
%%% Output fields:    
%%%   M               As above
%%%   Lz              As above
%%%   dz              Grid spacing (m)
%%%   Bz              As above
%%%   Ri0             Reference buoyancy Richardson number (M length vector)
%%%   K0              Reference diffusivity (M length vector, m^2/s)
%%%   KRi             Reference diffusivity gradient (M length vector, m^2/s)
%%%   f               As above
%%%   s               As above
%%%   A               Discrete eigenvalue problem coefficient matrix
%%%   B               Discrete eigenvalue problem weight matrix
%%%   growth          Vector of sorted growth rates 
%%%   V               Matrix of eigenvectors 
%%%   D               Diagonal matrix of eigenvalues 
%%%   
%%%
function result = diffusive_instability (params)



  %%% Default physical parameters    
  f = -1e-4;
  s = 1e-2;  
  
  use_diff = false;
  
  periodic_bc = true;
  
  %%% Check required parameters are present  
  if (~isfield(params,'M'))
    error('Input struct params must specify number of x-gridpoints M');
  end
  M = params.M;
  if (~isfield(params,'Lz'))
    error('Input struct params must specify domain width Lx');
  end
  Lz = params.Lz;
  if (~isfield(params,'Bz'))
    error('Input struct params must specify buoyancy frequency Bz');
  end
  Bz = params.Bz;
  if (~isfield(params,'l'))
    error('Input struct params must specify along-slope wavenumber l');
  end
  l = params.l;  
  
  %%% Check for optional physical parameters
  if (isfield(params,'f'))
    f = params.f;
  end
  if (isfield(params,'s'))
    s = params.s;
  end                 
  
  %%% Grids 
  dz = Lz/M;
  zz = -Lz/2+0.5*dz:dz:Lz/2-0.5*dz;  

  %%% Reference Richardson number
  Ri0 = f^2./(s^2*Bz);
  
  %%% Vertical diffusivity and its gradient w.r.t. Ri
  [K0,KRi] = calc_kappa(Ri0);   
   
  %%% Initialize result structure and add computed parameters
  result = struct;  
  result.f = f;
  result.Lz = Lz;
  result.dz = dz;
  result.M = M;  
  result.s = s;
  result.zz = zz;  
  result.Bz = Bz;
  result.Ri0 = Ri0;
  result.K0 = K0;
  result.KRi = KRi;
    
  %%% Coefficient matrix
  A = zeros(5*M);

  %%% Weight matrix
  B = zeros(5*M);  



  %%% Coefficients in u equation
  for i=1:M

    iref = i;
    
    %%% Horizontal advection term
    A(iref,M+i) = -s^2*Bz(i)/f;
    
    %%% Coriolis term
    A(iref,M+i) = A(iref,M+i) - f;
    
    %%% Vertical advection term
    if (periodic_bc)
      A(iref,2*M+i) = A(iref,2*M+i)+ 0.5 * s/f * Bz(i); 
      if (i == M)
        A(iref,2*M+1) = A(iref,2*M+1) + 0.5 * s/f * Bz(i);    
      else
        A(iref,2*M+i+1) = A(iref,2*M+i+1) + 0.5 * s/f * Bz(i);    
      end
    end
    
    %%% Diffusive terms 
    if (use_diff)

      if (i > 1)
  
        %%% Terms involving u'
        Kap_tmp = 0.5*(K0(i)+K0(i-1)) - 2 * 0.5*(KRi(i)*Ri0(i) + KRi(i-1)*Ri0(i-1));
        A(iref,iref) = A(iref,iref) + Kap_tmp / dz^2;
        A(iref,iref-1) = A(iref,iref-1) - Kap_tmp / dz^2;
  
        %%% Terms involving b'
        Kap_tmp = 0.5*(KRi(i)*Ri0(i) + KRi(i-1)*Ri0(i-1)) * (s/f);
        A(iref,3*M+i) = A(iref,3*M+i) + Kap_tmp / dz^2;
        A(iref,3*M+i-1) = A(iref,3*M+i-1) - Kap_tmp / dz^2;
  
      end
      if (i < M)
  
        %%% Terms involving u'
        Kap_tmp = 0.5*(K0(i)+K0(i+1)) - 2 * 0.5*(KRi(i)*Ri0(i) + KRi(i+1)*Ri0(i+1));
        A(iref,iref) = A(iref,iref) + Kap_tmp / dz^2;
        A(iref,iref+1) = A(iref,iref+1) - Kap_tmp / dz^2;
  
        %%% Terms involving b'
        Kap_tmp = 0.5*(KRi(i)*Ri0(i) + KRi(i+1)*Ri0(i+1)) * (s/f);
        A(iref,3*M+i) = A(iref,3*M+i) + Kap_tmp / dz^2;
        A(iref,3*M+i+1) = A(iref,3*M+i+1) - Kap_tmp / dz^2;
  
      end 

    end

    B(iref,iref) = 1;
  
  end  


  %%% Coefficients in v equation
  for i=1:M
    
    iref = M+i;

    %%% Coriolis
    A(iref,i) = f; 

    %%% Along-slope pressure gradient
    A(iref,4*M+i) = 1i*l; 

    %%% Vertical pressure gradient
    if (i == M)
      if (periodic_bc)
        A(iref,4*M+1) = A(iref,4*M+1) + -s * 1/(2*dz);
      else
        A(iref,4*M+i) = A(iref,4*M+i) + -s * 1/(2*dz);
      end
    else
      A(iref,4*M+i+1) = A(iref,4*M+i+1) + -s * 1/(2*dz);
    end
    if (i == 1)
      if (periodic_bc)
        A(iref,4*M+M) = A(iref,4*M+M) + -s * -1/(2*dz);
      else
        A(iref,4*M+i) = A(iref,4*M+i) + -s * -1/(2*dz);
      end
    else
      A(iref,4*M+i-1) = A(iref,4*M+i-1) + -s * -1/(2*dz);
    end    
    
    %%% Diffusive terms            
    if (use_diff)
      if (i > 1)
        Kap_tmp = 0.5*(K0(i)+K0(i-1));
        A(iref,iref) = A(iref,iref) + Kap_tmp / dz^2;
        A(iref,iref-1) = A(iref,iref-1) - Kap_tmp / dz^2;
      end
      if (i < M)
        Kap_tmp = 0.5*(K0(i)+K0(i+1));
        A(iref,iref) = A(iref,iref) + Kap_tmp / dz^2;
        A(iref,iref+1) = A(iref,iref+1) - Kap_tmp / dz^2;
      end    
    end

    B(iref,iref) = 1;
  
  end 
  

  
  
  %%% Coefficients in continuity equation
  for i=1:M

    iref = 2*M+i;
    
    %%% Along-slope gradient of v
    A(iref,M+i) = 1i*l;

    %%% Vertical gradient of v
    if (i == M)
      if (periodic_bc)
        A(iref,M+1) = A(iref,M+1) + -s * 1/(2*dz);
      else
        A(iref,M+i) = A(iref,M+i) + -s * 1/(2*dz);
      end
    else
      A(iref,M+i+1) = A(iref,M+i+1) + -s * 1/(2*dz);
    end
    if (i == 1)
      if (periodic_bc)
        A(iref,M+M) = A(iref,M+M) + -s * -1/(2*dz);
      else
        A(iref,M+i) = A(iref,M+i) + -s * -1/(2*dz);
      end
    else
      A(iref,M+i-1) = A(iref,M+i-1) + -s * -1/(2*dz);
    end

    %%% Vertical gradient of w
%     if (i == M)
%       if (periodic_bc)
%         A(iref,2*M+1) = A(iref,2*M+1) + 1/(2*dz);
%       else
%         A(iref,2*M+i) = A(iref,2*M+i) + 1/(2*dz);
%       end
%     else
%       A(iref,2*M+i+1) = A(iref,2*M+i+1) + 1/(2*dz);
%     end
%     if (i == 1)
%       if (periodic_bc)
%         A(iref,2*M+M) = A(iref,2*M+M) + -1/(2*dz);
%       else
%         A(iref,2*M+i) = A(iref,2*M+i) + -1/(2*dz);
%       end
%     else
%       A(iref,2*M+i-1) = A(iref,2*M+i-1) + -1/(2*dz);
%     end

    %%% Vertical gradient of w
    if (i == M)
      if (periodic_bc)
        A(iref,2*M+1) = A(iref,2*M+1) + 1/(dz);
      else
        A(iref,2*M+i) = A(iref,2*M+i) + 1/(dz); %%% TODO incorrect
      end
    else
      A(iref,2*M+i+1) = A(iref,2*M+i+1) + 1/(dz);
    end    

    A(iref,2*M+i) = A(iref,2*M+i) + -1/(dz);    
    
  end

  

  %%% Coefficients in b equation
  for i=1:M

    iref = 3*M+i;
    
    %%% Horizontal advection term
    A(iref,M+i) = -s*Bz(i);
    
    %%% Vertical advection term
    if (periodic_bc)
      A(iref,2*M+i) = A(iref,2*M+i) + 0.5 * Bz(i); 
      if (i == M)
        A(iref,2*M+1) = A(iref,2*M+1) + 0.5 * Bz(i);    
      else
        A(iref,2*M+i+1) = A(iref,2*M+i+1) + 0.5 * Bz(i);    
      end
    end
    
    %%% Diffusive terms    
    if (use_diff)

      if (i > 1)
  
        %%% Terms involving b'
        Kap_tmp = 0.5*(K0(i)+K0(i-1)) + 0.5*(KRi(i)*Ri0(i) + KRi(i-1)*Ri0(i-1));
        A(iref,iref) = A(iref,iref) + Kap_tmp / dz^2;
        A(iref,iref-1) = A(iref,iref-1) - Kap_tmp / dz^2;
  
        %%% Terms involving u'
        Kap_tmp = 0.5*(KRi(i)*Ri0(i) + KRi(i-1)*Ri0(i-1)) * (-2*f/s);
        A(iref,i) = A(iref,i) + Kap_tmp / dz^2;
        A(iref,i-1) = A(iref,i-1) - Kap_tmp / dz^2;
  
      end
      if (i < M)
  
        %%% Terms involving b'
        Kap_tmp = 0.5*(K0(i)+K0(i+1)) + 0.5*(KRi(i)*Ri0(i) + KRi(i+1)*Ri0(i+1));
        A(iref,iref) = A(iref,iref) + Kap_tmp / dz^2;
        A(iref,iref+1) = A(iref,iref+1) - Kap_tmp / dz^2;
  
        %%% Terms involving u'
        Kap_tmp = 0.5*(KRi(i)*Ri0(i) + KRi(i+1)*Ri0(i+1)) * (-2*f/s);
        A(iref,i) = A(iref,i) + Kap_tmp / dz^2;
        A(iref,i+1) = A(iref,i+1) - Kap_tmp / dz^2;
  
      end 

    end

    B(iref,iref) = 1;
  
  end  

 

  
  
  %%% Coefficients in hydrostatic equation    
  for i=1:M

    iref = 4*M+i;
    
%     %%% Vertical pressure gradient
%     if (i == M)
%       if (periodic_bc)
%         A(iref,4*M+1) = A(iref,4*M+1) + 1/(2*dz);
%       else
%         A(iref,4*M+i) = A(iref,4*M+i) + 1/(2*dz);
%       end
%     else
%       A(iref,4*M+i+1) = A(iref,4*M+i+1) + 1/(2*dz);
%     end
%     if (i == 1)
%       if (periodic_bc)
%         A(iref,4*M+M) = A(iref,4*M+M) + -1/(2*dz);
%       else
%         A(iref,4*M+i) = A(iref,4*M+i) + -1/(2*dz);
%       end
%     else
%       A(iref,4*M+i-1) = A(iref,4*M+i-1) + -1/(2*dz);
%     end

    %%% Vertical pressure gradient
    if (i == M)
      if (periodic_bc)
        A(iref,4*M+1) = A(iref,4*M+1) + 1/(dz);
      else
        A(iref,4*M+i) = A(iref,4*M+i) + 1/(dz);
      end
    else
      A(iref,4*M+i+1) = A(iref,4*M+i+1) + 1/(dz);
    end

    A(iref,4*M+i) = A(iref,4*M+i) + -1/(dz);    

    %%% Buoyancy anomaly
    A(iref,3*M+i) = -1;

    %%% Buoyancy anomaly
    if (i == M)
      A(iref,3*M+1) = A(iref,3*M+1) + -0.5;
    else
      A(iref,3*M+i+1) = A(iref,3*M+i+1) + -0.5;
    end
    A(iref,3*M+i) = A(iref,3*M+i) + -0.5;

  end





  %%% Add the coefficient matrix to the result structure
  result.A = A;

  %%% Add the weight matrix to the result structure
  result.B = B;
  
  
  
  %%% Determine eigenvalues and sort growth rates
  omega = eig(A,B)/1i;
  growth = (omega - real(omega)).*-1i;
  growth = sort(growth,'descend');
  result.growth = growth;  
  
  %%% Extract eigenvectors   
  [V,D]=eig(A,B);
  result.V = V;
  result.D = D;  
  
end

