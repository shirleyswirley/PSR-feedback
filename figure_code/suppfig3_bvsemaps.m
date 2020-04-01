close all;
clear all;

%-------------------------------
% Set up paths
%-------------------------------
utils_path = '/graid1/shirlleu/PSRfeedback/utils/';
data_path = '/graid1/shirlleu/PSRfeedback/data/';
data_path_bvse = '/graid1/shirlleu/PSRfeedback/data/globalbetavsexportmaps_tempregress/';
fig_save_path = '/graid1/shirlleu/PSRfeedback/pdfs_pngs/';

%-------------------------------
% Set up warnings and utils
%-------------------------------
warning('off','all');
addpath(genpath(utils_path));

%-------------------------------
% Load all export maps + grid variables
% (everything on satellite grid since
% this is how everything was computed)
%-------------------------------

% - Define the eratio and npp algs you want to use
npp_algs = {'VGPM','VGPME','CbPM'};
eratio_algs = {'Laws2000','D2005PP','Laws2011D'};
eratio_algs_sh = {'L2000','D2005','L2011'};

% - Initialize matrix to hold all
% beta vs. export maps (bvse_maps)
eration = length(eratio_algs);
nppn = length(npp_algs);
timen = 160; % num of months in time series
load([data_path_bvse 'betavsexportn1tempregmo_VGPMVGPMECbPMcompiled.mat'],'lat1','lon1')
lon1 = [lon1(181:end)' lon1(1:180)'+360];
bvse_maps = nan(length(lat1),length(lon1),eration,nppn);

% - Load b vs e maps + populate bvse_maps
for ieratio = 1:eration
    for inpp = 1:nppn
        load([data_path_bvse npp_algs{inpp} 'based/Kostbetavs' ...
            eratio_algs{ieratio} 'exportn1tempregmo.mat'],...
            'betavsetempregresscoeffmapzeros');
        bvse_maps(:,:,ieratio,inpp) = betavsetempregresscoeffmapzeros;
    end
end

% - Put b vs. e maps on 2-deg PRiSM grid
load([data_path 'PRiSM_regions_2deg.mat'],'gridd','M3d');
lon2 = gridd.xt; lat2 = gridd.yt;
lon2grid = gridd.XT; lat2grid = gridd.YT;
[lon1grid,lat1grid] = meshgrid(lon1,lat1);
bvse_maps1 = [bvse_maps(:,181:end,:,:) bvse_maps(:,1:180,:,:)];
bvse_maps2 = nan(length(lat2),length(lon2),eration,nppn);
for ieratio = 1:eration
    for inpp = 1:nppn
        bvsetemp = bvse_maps1(:,:,ieratio,inpp);
        bvse_maps2(:,:,ieratio,inpp) = interp2(lon1grid,lat1grid,bvsetemp,lon2grid,lat2grid);
    end
end

%-------------------------------
% Plot figure
%-------------------------------

% - Define plot params
cmap = flipud(cbrewer('div','RdBu',21,'linear'));
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
labelfontsize = 10;
titlefontsize = 12;
titlefontwt = 'bold';
cbticklen = 0.03;
cmin = -1.4; cmax = 1.4;
qual9cmap = cbrewer('qual','Set1',9);
qual9bgcmap = ones(size(qual9cmap));
qual9bgcmap(6,:) = [0.8 0.8 0.8];
% gives yellow text a gray bg

% - Plot each b vs e map in turn
f=figure;
set(f,'color','white','units','inches','position',[0.5 0.5 11 6],'resize','off');
isubplot =  1;
for ieratio = 1:eration
    for inpp = 1:nppn
        subplot(3,3,isubplot);
        m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
        m_pcolor(lon2,lat2,bvse_maps2(:,:,ieratio,inpp));
        m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
        caxis([cmin cmax]); colormap(cmap); shading flat;
        cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
        mapvarnow = bvse_maps2(:,:,ieratio,inpp); extend_cbar_ticklabels; % needs cmin,cmax,mapvarnow,cb
        m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
        title([sprintf('%s{%f %f %f}','\color[rgb]',qual9cmap(isubplot,:)) ...
            '\beta_ vs. ' npp_algs{inpp} '+' eratio_algs_sh{ieratio} ' export'],...
            'backgroundcolor',qual9bgcmap(isubplot,:),'margin',1,...
            'verticalalignment','bottom','fontsize',titlefontsize,...
            'fontweight',titlefontwt);
        isubplot = isubplot+1;
    end
end

print(f, [fig_save_path 'suppfig3_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'suppfig3_PSRfbpaper_final.png'], '-dpng', '-r300');
