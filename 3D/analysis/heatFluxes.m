%%% Calculate mean and eddy heat fluxes
% mean velocities and potential temperature
mean_uvel = mean(uvel,4);
mean_vvel = mean(vvel,4);
mean_wvel = mean(wvel,4);
mean_theta = mean(theta,4);
% mean advective heat flux
mean_uadv_hf = mean_uvel.*mean_theta;
mean_vadv_hf = mean_vvel.*mean_theta;
mean_wadv_hf = mean_wvel.*mean_theta;
% mean eddy heat flux / turbulent flux
eddy_uvel_hf = mean(uvel.*theta,4) - mean_uvel;
eddy_vvel_hf = mean(vvel.*theta,4) - mean_vvel;
eddy_wvel_hf = mean(wvel.*theta,4) - mean_wvel;
% mean diffusive flux -- make sure derivative is aligned with diffusivities
dmthetadx(2:Nx,:,:,:) = diff(mean_theta,1,1)./dx; dmthetadx(1,:,:,:) = 2*dmthetadx(2,:,:,:)-dmthetadx(3,:,:,:);
dmthetady(:,2:Ny,:,:) = diff(mean_theta,1,2)./dy; dmthetady(:,1,:,:) = 2*dmthetady(:,2,:,:)-dmthetady(:,3,:,:);
dmthetadz(:,:,:,1:Nr-1) = (mean_theta(:,:,:,1:Nr-1)-mean_theta(:,:,:,2:Nr))./dz; dmthetadz(:,:,:,Nr) = 2*dmthetadz(:,:,:,Nr-1)-dthetadz(:,:,:,Nr-2);
mean_udiff_hf = diffT.*dmthetadx;
mean_vdiff_hf = diffT.*dmthetady;
mean_wdiff_hf = diffT.*dmthetadz;
% mean heat flux
mean_u_hf = mean_uadv_hf + eddy_uvel_hf - mean_udiff_hf;
mean_v_hf = mean_vadv_hf + eddy_vvel_hf - mean_vdiff_hf;
mean_w_hf = mean_wadv_hf + eddy_wvel_hf - mean_wdiff_hf;
% total mean heat flux
tot_mean_hf = mean_u_hf + mean_v_hf + mean_w_hf;
