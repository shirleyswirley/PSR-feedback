plotonly=0;

if plotonly==0

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
% Compute map of # of non-nan months
%-------------------------------
% NOTE: something is wrong with this...
be_nn_mon_maps = nan(length(lat1),length(lon1),nppn,eration);
for ieratio = 1:eration
    for inpp = 1:nppn
        for ilat = 1:length(lat1)
            for ilon = 1:length(lon1)
                bnow = squeeze(beta_mon_map(ilat,ilon,:));
                enow = squeeze(cexp_mon_maps(ilat,ilon,:,inpp,ieratio));
                nancount = sum(~isnan([bnow, enow]),2);
                be_nn_mon_maps(ilat,ilon,inpp,ieratio) = sum(nancount==2);
            end
        end
    end
end
bennmapnow = squeeze(be_nn_mon_maps(:,:,1,2));
figure;
pcolor(lon1,lat1,squeeze(be_nn_mon_maps(:,:,1,1)));
shading flat; colorbar;

%-------------------------------
% Regrid region map onto 1x1 satellite grid
%-------------------------------
reg_struct = load([data_path ...
    'PRiSM_regions_2deg.mat']);
lon2 = reg_struct.gridd.xt; lat2 = reg_struct.gridd.yt;
regmap = reg_struct.R2d;
%lont = [lon2(180)-360 lon2(181:end)-360 lon2(1:180)];
%[lontgrid,lattgrid]=meshgrid(lont,lat2);
%regmap = [regmap(:,180) regmap(:,181:end) regmap(:,1:180)];
%regmap = interp2(lontgrid,lattgrid,regmap,lon1grid,lat1grid,'nearest');

lonidxs = [263, 339, 52, 160, 13, 301];
latidxs = [95, 40, 91, 119, 134, 150];
figure;
%pcolor(lon1,lat1,nanmean(beta_mon_map,3));
pcolor(lon2,lat2,regmap);
hold on;
for iidx = 1:length(lonidxs)
    scatter(lon1(lonidxs(iidx)),lat1(latidxs(iidx)),'k.');
end
shading flat; colorbar;

%-------------------------------
% Count # of non-nan months at all lon at a given lat
% (or vice versa)
%-------------------------------
validlat = nan(1,length(lat1));
for ilat = 1:length(lat1)
    lonidx = 300; latidx = ilat;
    enow = squeeze(cexp_mon_maps(latidx,lonidx,:,1,1));
    bnow = squeeze(beta_mon_map(latidx,lonidx,:));
    nancount = sum(~isnan([bnow, enow]),2);
    validlat(ilat) = sum(nancount==2);
    %disp(ilat)
    %disp(sum(nancount==2))
end

validlon = nan(1,length(lon1));
for ilon = 1:length(lon1)
    lonidx = ilon; latidx = 150;
    enow = squeeze(cexp_mon_maps(latidx,lonidx,:,1,1));
    bnow = squeeze(beta_mon_map(latidx,lonidx,:));
    nancount = sum(~isnan([bnow, enow]),2);
    validlon(ilon) = sum(nancount==2);
    %disp(ilon)
    %disp(sum(nancount==2))
end
% --> latidx = 40 (51N): many
% --> latidx = 70 (21N): many
% --> latidx = 150: 116, 301 + others
% --> latidx = 153: 124, 304 + others
% --> latidx = 156: 237, 96
% --> latidx = 165: none

end

%-------------------------------
% Plot figures
%-------------------------------
lonidx = 304; latidx = 153; regname = 'AAZ';

months = datetime(kost.Xitimemonames,'ConvertFrom','datenum');
qual9cmap = cbrewer('qual','Set1',9);
qual9bgcmap = ones(size(qual9cmap));
qual9bgcmap(6,:) = [0.8 0.8 0.8];
titlefontsize = 8;
titlefontwt = 'bold';
sgtitlefontsize = 10;

figure;
%pcolor(lon1,lat1,nanmean(beta_mon_map,3));
pcolor(lon1,lat1,regmap);
hold on;
scatter(lon1(lonidx),lat1(latidx),'k.');
shading flat; colorbar;

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
sgtitle(['lon=' num2str(lon1(lonidx)) ', lat=' num2str(lat1(latidx))], 'fontsize', sgtitlefontsize);
print(f, [fig_save_path 'betaexporttimeseries_' regname '_lon1idx' num2str(lonidx) ...
    '_lat1idx' num2str(latidx) '.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'betaexporttimeseries_' regname '_lon1idx' num2str(lonidx) ...
    '_lat1idx' num2str(latidx) '.png'], '-dpng', '-r300');

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
        %xlim([datetime(1998,1,1), datetime(2005,1,1)]);
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
sgtitle(['lon=' num2str(lon1(lonidx)) ', lat=' num2str(lat1(latidx))], 'fontsize', sgtitlefontsize);
print(f, [fig_save_path 'betaexporttimeseriesfull_' regname '_lon1idx' num2str(lonidx) ...
    '_lat1idx' num2str(latidx) '.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'betaexporttimeseriesfull_' regname '_lon1idx' num2str(lonidx) ...
    '_lat1idx' num2str(latidx) '.png'], '-dpng', '-r300');

f=figure;
set(f,'color','white','units','inches','position',[0.5 0.5 10 7],'resize','off');
isubplot =  1;
for ieratio = 1:eration
    for inpp = 1:nppn
        subplot(3,3,isubplot);
        scatter(cexp_mon_maps(latidx,lonidx,:,inpp,ieratio)...
                ./nanmean(cexp_mon_maps(latidx,lonidx,:,inpp,ieratio)), ...
                beta_mon_map(latidx,lonidx,:));
        title([sprintf('%s{%f %f %f}','\color[rgb]',qual9cmap(isubplot,:)) ...
            '\beta_ vs. ' npp_algs{inpp} '+' eratio_algs_sh{ieratio} ' export'],...
            'backgroundcolor',qual9bgcmap(isubplot,:),'margin',1,...
            'verticalalignment','bottom','fontsize',titlefontsize,...
            'fontweight',titlefontwt);
        isubplot = isubplot + 1;
    end
end
sgtitle(['lon=' num2str(lon1(lonidx)) ', lat=' num2str(lat1(latidx))], 'fontsize', sgtitlefontsize);

