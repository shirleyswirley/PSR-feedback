close all;
clear all;
setup_figs;

%-------------------------------
% Load data + grid variables
%-------------------------------
% - Load grid variables
reg_struct = load([data_path 'PRiSM_regions_2deg.mat'],'M3d','gridd');
M3d = reg_struct.M3d;
gridd = reg_struct.gridd;
lon2 = gridd.xt;
lat2 = gridd.yt;

% - Load Kost beta map
load([data_path 'Kost_beta_2deg.mat'],'beta_clim');
kost_beta_map = mean(beta_clim,3)';
kost_beta_map(M3d(:,:,1)==0)=NaN;

% - Load/create remin depth map based on Kost beta
beta_name = 'kost';
load_rd_map=1;
if load_rd_map==1
    load([data_path 'remin_depth_map_' beta_name 'beta.mat'],...
        'remin_depth_map');
    kost_rd_map = remin_depth_map;
else
    var_des = ...
        'remin depths in [m] calc using PRiSM from annual mean Kost beta map';
    kost_rd_map = calc_remin_depth_map(...
        lon2,lat2,kost_beta_map,beta_name,var_des);
end

%-------------------------------
% Plot figure
%-------------------------------
% - Define plot params
seqcmap=cbrewer('seq','YlGnBu',20,'linear');
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
labelfontsize = 10;
titlefontsize = 12;
titlefontwt = 'bold';
cbticklen = 0.03;

f=figure;
set(f,'color','white','units','inches',...
    'position',[1 1 8.5 2.5],'resize','off','paperpositionmode','auto');

% - Kost beta map
ax = subplot(121);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,kost_beta_map,'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
colormap(ax,seqcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('Annual mean {\itβ_ }','fontsize',titlefontsize,'fontweight',titlefontwt);

% - Remin depth map
ax = subplot(122);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,log10(kost_rd_map),'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
colormap(ax,flipud(seqcmap)); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
cb.Ticks = [2, 2.5, 3]; cb.TickLabels = {'10^2', '10^{2.5}', '10^3'};
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('Remineralization depth [m]_ ','fontsize',titlefontsize,'fontweight',titlefontwt);

print(f, [fig_save_path 'fig1_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'fig1_PSRfbpaper_final.png'], '-dpng', '-r300');
