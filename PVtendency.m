%%% Calculating PV tendency using MITgcm model diagnostics 

%%% need u,b and udot, bdot (dot = tendency).
%%% Um_Diss - U momentum tendency from Dissipation (Explicit part)
%%% Um_Impl - U momentum tendency from Dissipation (Implicit part)
%%% Um_Advec - U momentum tendency from Advection terms
%%% Um_Cori - U momentum tendency from Coriolis term
%%% Um_dPhiX - U momentum tendency from Pressure/Potential grad
%%% USidDrag - U momentum tendency from Side Drag
%%% ADVx_TH - Zonal Advective Flux of Pot.Temperature
%%% DFxE_TH - Zonal Diffusive Flux of Pot.Temperature
%%% ADVx_SLT - Zonal Advective Flux of Salinity
%%% DFxE_SLT - Zonal Diffusive Flux of Salinity
%%% TOTUTEND - Tendency of Zonal Component of Velocity
%%% TOTTTEND - Tendency of Potential Temperature
%%% TOTSTEND - Tendency of Salinity

%%% define b

%%% define bdot

%%% define u

%%% define udot

%%% term1 
A = mean(f0*bdot,2,'omitnan');
term1 = (diff(A,1,3))./dz;
% recenter
% interpolate first/last idx

%%% term2 
A = mean(diff(udot,1,2),3);
B = mean(diff(b,1,3),2);
term2 = -(A.*B)./dy./dz;

%%% term3
A = mean(diff(u,1,2),3);
B = mean(diff(bdot,1,3),2);
term3 = -(A.*B)./dy./dz;

%%% term4 
A = mean(diff(udot,1,3),2);
B = mean(diff(b,1,2),3);
term4 = (A.*B)./dy.dz;

%%% term5 
A = mean(diff(u,1,3),2);
B = mean(diff(bdot,1,2),3);
term5 = (A.*B)./dy./dz;

%%% qdot
qdot = term1 + term2 + term3 + term4 + term5; 








