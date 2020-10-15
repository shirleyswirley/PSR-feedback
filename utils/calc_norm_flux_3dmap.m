function [norm_flux_map, z] = calc_norm_flux_3dmap(lon,lat,beta_map,beta_name,var_des)

% 1.) Constants from DeVries et al., 2014 mostly table 1
z_surf = 0:10:100; % [m]
z_deep = 250:250:1000; % [m]
z = [z_surf z_deep]; % water column depth range from bottom of euphotic zone to bottom of ocean
c_w = 2.2E5; % coeff in particle sinking rate equation as a fxn of D
c_r = 0.03; % degradation rate of sinking particles
c_1 = c_r/c_w;
DS = 20E-6; % smallest diameter of particles
DL_0 = 2000E-6; % largest diameter of particles in the euphotic zone
eta = 1.17; % exponent in particle sinking rate equation as a fxn of D
zeta = 2.28; % exponent in mass equation as a fxn of D

% 2.) Use DeVries et al., 2014 equations to create POC flux profiles from beta
norm_flux_map=nan(length(lat),length(lon),length(z));
    
for ilon=1:length(lon)
    for ilat=1:length(lat)
        if ~isnan(beta_map(ilat,ilon))
            for iz=1:length(z)
                zj=z(iz);
                % equation 11
                DL = (max(DL_0^eta-c_1*eta/zeta*zj,0))^(1/eta);
                % equation 13
                F = @(D)(D.^(zeta+eta-beta_map(ilat,ilon))).*(1+(c_1*eta/zeta).*(D.^-eta)*zj).^((1-beta_map(ilat,ilon))/eta);
                if (iz==1); F_0 = quad(F,DS,DL); end; % save the particle flux at surface to use for normalization 
                norm_flux_map(ilat,ilon,iz) = quad(F,DS,DL)/F_0;
            end
        end
    end
end

save(['norm_poc_flux_3dmap_' beta_name 'beta.mat'],'var_des','norm_flux_map','lon','lat','z','c_w','c_r','DS','DL_0','eta','zeta');
