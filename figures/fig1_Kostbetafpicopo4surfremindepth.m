close all;
clear all;

%-------------------------------
% Set up paths
%-------------------------------
utils_path = '/graid1/shirlleu/PSRfeedback/utils/';
data_path = '/graid1/shirlleu/PSRfeedback/data/';

%-------------------------------
% Set up warnings and utils
%-------------------------------
warning('off','all');
addpath(genpath(utils_path));

%-------------------------------
% Load data + grid variables
%-------------------------------

% - Load regions map + grid variables
reg_struct = load([data_path 'PRiSM_regions_2deg.mat'],'M3d','gridd','R2d','regnamesabbrev');
M3d = reg_struct.M3d;
gridd = reg_struct.gridd;
lon2 = gridd.xt;
lat2 = gridd.yt;
reg_map = reg_struct.R2d;
reg_names = reg_struct.regnamesabbrev;

% - Load Kost beta map
load([data_path 'Kost_beta_2deg.mat'],'beta_clim');
kost_beta_map = mean(beta_clim,3)';
kost_beta_map(M3d(:,:,1)==0)=NaN;

% - Load WOA po4 map
load([data_path 'WOA_2009_2deg.mat'],'po4an');
po4_surf_map = po4an(:,:,1);

% - Load fpico map
load([data_path 'fpico_monthly_1deg.mat'],'freqpftsmo1','lat1','lon1');
fpico_map = nanmean(freqpftsmo1(:,:,:,1),3);
fpico_map = [fpico_map(:,181:end) fpico_map(:,1:180)];
lon1 = [lon1(181:end)' lon1(1:180)'+360];

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
plotregcontours=0;
seqcmap=cbrewer('seq','YlGnBu',20,'linear');
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
labelfontsize = 10;
titlefontsize = 12;
titlefontwt = 'bold';
cbticklen = 0.03;

f=figure;
set(f,'color','white','units','inches','position',[1 1 8.5 5],'resize','off','paperpositionmode','auto');

% - Kost beta map
ax = subplot(221);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,kost_beta_map,'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
colormap(ax,seqcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('Annual mean β_ ','fontsize',titlefontsize,'fontweight',titlefontwt);

% - fpico map
ax = subplot(222);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat1) max(lat1)])
m_contourf(lon1,lat1,fpico_map,'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
colormap(ax,seqcmap); shading flat; cb = colorbar;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('Annual mean f_{pico}','fontsize',titlefontsize,'fontweight',titlefontwt);

% - PO4 map
ax=subplot(223);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,po4_surf_map,'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
%caxis([0 1.5]);
colormap(ax,flipud(seqcmap)); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('Annual mean surface PO_4 [mmol/m^3]','fontsize',titlefontsize,'fontweight',titlefontwt);

% - Remin depth map
ax = subplot(224);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,log10(kost_rd_map),'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
colormap(ax,flipud(seqcmap)); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('log(Remineralization depth [m])_ ','fontsize',titlefontsize,'fontweight',titlefontwt);

print(f, 'fig1_PSRfbpaper_final.pdf', '-dpdf', '-r300');
print(f, 'fig1_PSRfbpaper_final.png', '-dpng', '-r300');
