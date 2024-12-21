%%%
%%% Computes linear growth rates over a range of y and z wave numbers.
%%% 
function [gr,lamY,lamZ,Ri0] = analytical_soln (N2,s)

  %%% Parameters
  lamY = 10.^[1:.01:5];
  lamZ = 10.^[-1:.01:3];
  ll = 2*pi./lamY; %%% Wavenumber vector
  mm = 2*pi./lamZ;
  f = -1e-4;
  if s==0
    s = 1e-2;
  end
  Ri0 = f^2./(s^2*N2);
  [K0,KRi] = calc_kappa (Ri0);
  
  gr = zeros(length(ll),length(mm));
  for i=1:length(ll)
    for j=1:length(mm)
  
      l = ll(i);
      m = mm(j);
  
      c3 = 1;
      c2 = -KRi*Ri0*m^2;
      c1 = f^2 + N2*(l/m)^2;
      c0 = -KRi*Ri0*m^2*f^2*(1 -2*l/m/s);
      
      rr = roots([c3 c2 c1 c0]);
      gr(i,j) = max(real(rr) - m^2*K0);
  
  
      %%% Alternative version via full matrix solve
      %%% Solution vector [u,v,w,b,phi]
%       A = [(K0-2*KRi*Ri0)*m^2,    -f-s^2*N2/f,    s*N2/f,     (s/f)*KRi*Ri0*m^2,      0;
%            f,                     m^2*K0,         0,          0,                      1i*(l-s*m); ...
%            0,                     0,              0,          -1,                     1i*m; ...
%            -2*(f/s)*KRi*Ri0*m^2,  -s*N2,          N2,         (K0+KRi*Ri0)*m^2,       0;
%            0,                     l-s*m,          m,          0,                      0];
%       B = [1,     0,    0,    0,    0; ...
%            0,     1,    0,    0,    0; ...
%            0,     0,    0,    0,    0; ...
%            0      0,    0,    1,    0;
%            0      0,    0,    0,    0];
%   
%       the_eigs = eig(A,B);
%       the_eigs(isinf(the_eigs)) = [];
%       gr(i,j) = max(real(-the_eigs));
  
    end
  end

end






