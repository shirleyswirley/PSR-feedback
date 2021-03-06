clear all;
close all;
setup_figs;

%-------------------------------
% Load all export maps + grid variables
%-------------------------------
% - List all npp and e-ratio algorithms
npp_algs = {'VGPM','VGPME','CbPM'};
eratio_algs = {'Laws2000','D2005PP','Laws2011D'};
eratio_algs_abbrevs = {'L2000','D2005','L2011'};

% - Initialize matrix to hold all monthly mean
% export maps (cexp_mon_maps)
eration = length(eratio_algs);
nppn = length(npp_algs);
timen = 160; % length of monthly mean export time series
load([data_path 'globalexportmaps_spatvareratioandezdepth/expCmo_D2005PPfromCbPM.mat'],'lat1','lon1');
lon1 = [lon1(181:end)' lon1(1:180)'+360];
cexp_mon_maps = nan(length(lat1),length(lon1),timen,nppn,eration);

% - Load export maps + populate cexp_mon_maps
for inpp = 1:nppn
    for ieratio = 1:eration
        load([data_path 'globalexportmaps_spatvareratioandezdepth/expCmo_' ...
            eratio_algs{ieratio} 'from' npp_algs{inpp} '.mat'],'expCmo1');
        disp([data_path 'globalexportmaps_spatvareratioandezdepth/expCmo_' ...
            eratio_algs{ieratio} 'from' npp_algs{inpp} '.mat']);
        cexp_mon_maps(:,:,:,inpp,ieratio) = ...
            [expCmo1(:,181:end,:) expCmo1(:,1:180,:)];
    end
end

%-------------------------------
% Plot figure
%-------------------------------
% - Define plot params
cmap = flipud(cbrewer('seq','YlGnBu',20,'linear'));
mapproj = 'gall-peters'; landcolor = [0.6 0.6 0.6];
labelfontsize = 10;
titlefontsize = 12; titlefontwt = 'bold';
cbticklen = 0.03; cbnumticks = 5;
cmin = 0; cmax = 8;

% - Plot each annual mean global export map in turn
f=figure;
set(f,'color','white','units','inches','position',[0.5 0.5 11 6]);
isubplot=1;
for ieratio = 1:eration
    for inpp = 1:nppn
        cexp_mapnow = nanmean(cexp_mon_maps(:,:,:,inpp,ieratio),3);
        subplot(3,3,isubplot);
        m_proj(mapproj,'lon',[0 360],'lat',[min(lat1) max(lat1)])
        m_pcolor(lon1,lat1,cexp_mapnow);
        m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
        caxis([cmin cmax]); colormap(cmap); shading flat;
        cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
        mapvarnow = cexp_mapnow; extend_cbar_ticklabels; % needs cmin,cmax,mapvarnow,cb,cbnumticks
        m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
        title([npp_algs{inpp} ' NPP+' eratio_algs_abbrevs{ieratio} ' e-ratio_ '],...
            'fontsize',titlefontsize,'fontweight',titlefontwt);
        isubplot=isubplot+1;
    end
end

print(f, [fig_save_path 'suppfig2_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'suppfig2_PSRfbpaper_final.png'], '-dpng', '-r300');
