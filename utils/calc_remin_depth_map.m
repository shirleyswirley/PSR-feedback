function remin_depth_map = calc_remin_depth_map(lon,lat,beta_map,beta_name,var_des)

% 1.) Constants from DeVries et al., 2014 mostly table 1
z_0=0; % [m]
z_f=5000; % [m]
z_res=50; % [m]
z=z_0:z_res:z_f; % water column depth range from bottom of euphotic zone to bottom of ocean
c_w = 2.2E5; % coeff in particle sinking rate equation as a fxn of D
c_r = 0.03; % degradation rate of sinking particles
c_1 = c_r/c_w;
DS = 20E-6; % smallest diameter of particles
DL_0 = 2000E-6; % largest diameter of particles in the euphotic zone
eta = 1.17; % exponent in particle sinking rate equation as a fxn of D
zeta = 2.28; % exponent in mass equation as a fxn of D

% 2.) Use DeVries et al., 2014 equations to create POC flux profiles from beta 
% Then using those profiles, find the depth at which the flux fraction of the
% surface is equal to what you specify
Fcurves=nan(length(lat),length(lon),length(z));
    
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
                Fcurves(ilat,ilon,iz) = quad(F,DS,DL)/F_0;
            end
        end
    end
end

% 3.) Then using those particle flux profiles,
% find the depth at which the flux fraction of the
% surface is equal to what you specify
flux_frac = 1/exp(1); %0.5;
remin_depth_map = nan(length(lat),length(lon));

for ilon=1:length(lon)
    for ilat=1:length(lat)
        if ~isnan(beta_map(ilat,ilon))
            [jnk vertidx]=min(abs(Fcurves(ilat,ilon,:)-flux_frac));
            if Fcurves(ilat,ilon,vertidx)>flux_frac
                remin_depth_map(ilat,ilon) = z(vertidx) - (Fcurves(ilat,ilon,vertidx)-flux_frac)*(z(vertidx)-z(vertidx+1))/(Fcurves(ilat,ilon,vertidx)-Fcurves(ilat,ilon,vertidx+1));
            elseif Fcurves(ilat,ilon,vertidx)<flux_frac
                remin_depth_map(ilat,ilon) = z(vertidx) + (flux_frac-Fcurves(ilat,ilon,vertidx))*(z(vertidx-1)-z(vertidx))/(Fcurves(ilat,ilon,vertidx-1)-Fcurves(ilat,ilon,vertidx));
            else % (Fcurves(ilat,ilon,vertidx)-0.5)==0
                remin_depth_map(ilat,ilon) = z(vertidx);
            end        
        end
    end
end

save(['remin_depth_map_' beta_name 'beta.mat'],'var_des','remin_depth_map','lon','lat','z','c_w','c_r','DS','DL_0','eta','zeta','flux_frac');
