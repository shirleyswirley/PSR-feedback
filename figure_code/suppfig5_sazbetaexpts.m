clear all;
close all;
setup_figs;

timen = 160; % length of monthly mean export and beta time series

%-------------------------------
% Load all export maps + grid variables
%-------------------------------
% - List all npp and e-ratio algorithms
npp_algs = {'VGPM','VGPME','CbPM'};
eratio_algs = {'Laws2000','D2005PP','Laws2011D'};
eratio_algs_sh = {'L2000','D2005','L2011'};

% - Initialize matrix to hold all monthly mean
% export maps (cexp_mon_maps)
eration = length(eratio_algs);
nppn = length(npp_algs);
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
% Load + regrid beta map
%-------------------------------
kost = load([data_path 'Kost_beta_monthly_1deg.mat']);
beta_mon_map = nan(length(lat1), length(lon1), timen); 
[lon1grid,lat1grid] = meshgrid(lon1,lat1);
lonk = [kost.lon1(181:end)' kost.lon1(1:180)'+360];
[lonkgrid,latkgrid] = meshgrid(lonk,kost.lat1);
for imo = 1:timen
    beta_mon_map(:,:,imo) = interp2(lonkgrid, latkgrid, kost.Ximo1(:,:,imo), lon1grid, lat1grid);
end 

%-------------------------------
% Regrid region map onto 1x1 satellite grid
%-------------------------------
reg_struct = load([data_path ...
    'PRiSM_regions_2deg.mat']);
lon2 = reg_struct.gridd.xt; lat2 = reg_struct.gridd.yt;
regmap = reg_struct.R2d;

%-------------------------------
% Plot figures
%-------------------------------
lonidx = 13; latidx = 134; regname = 'SAZ';

months = datetime(kost.Xitimemonames,'ConvertFrom','datenum');
qual9cmap = cbrewer('qual','Set1',9);
qual9bgcmap = ones(size(qual9cmap));
qual9bgcmap(6,:) = [0.8 0.8 0.8];
titlefontsize = 8; titlefontwt = 'bold';
sgtitlefontsize = 12; sgtitlefontwt = 'bold';

f=figure;
set(f,'color','white','units','inches','position',[0.5 0.5 12 7],'resize','off');
isubplot =  1;
for ieratio = 1:eration
    for inpp = 1:nppn
        ax=subplot(3,3,isubplot);
        yyaxis left;
        plot(months,...
            squeeze(cexp_mon_maps(latidx,lonidx,:,inpp,ieratio)...
            ./nanmean(cexp_mon_maps(latidx,lonidx,:,inpp,ieratio))));
        yyaxis right;
        plot(months, squeeze(beta_mon_map(latidx,lonidx,:)));
        xlim([datetime(1999,1,1), datetime(2005,1,1)]);
        xtickformat('MMM-yy'); xtickangle(45); grid on;
        ax.XAxis.FontSize = 8;
        title([sprintf('%s{%f %f %f}','\color[rgb]',qual9cmap(isubplot,:)) ...
            '\beta_ vs. ' npp_algs{inpp} '+' eratio_algs_sh{ieratio} ' export'],...
            'backgroundcolor',qual9bgcmap(isubplot,:),'margin',1,...
            'verticalalignment','bottom','fontsize',titlefontsize,...
            'fontweight',titlefontwt);
        isubplot = isubplot + 1;
    end
end
sgtitle(['SAZ example (lon=' num2str(round(lon1(lonidx))) '°, lat=' num2str(round(lat1(latidx))) '°)'], 'fontsize', sgtitlefontsize, 'fontweight', sgtitlefontwt);

print(f, [fig_save_path 'suppfig5_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'suppfig5_PSRfbpaper_final.png'], '-dpng', '-r300');
