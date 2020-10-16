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
M3d = reg_struct.M3d;
reg_map = reg_struct.R2d;

%-------------------------------
% Define and plot map of example locations
%-------------------------------
regnames = {'IND','NA','ETP','STA','SAZ','AAZ'};
lonidxs = [263, 339, 52, 160, 13, 301];
latidxs = [95, 40, 91, 119, 134, 150];

% - Define plot params
seqcmap = cbrewer('seq','YlGnBu',10,'linear');
mapproj = 'gall-peters'; landcolor = [0.9 0.9 0.9];
reglinewidth = 2; reglinecolor = [0.5 0.5 0.5];
labelfontsize = 10;
titlefontsize = 10; titlefontwt = 'bold';
cbticklen = 0.03;

f=figure;
set(f,'color','white','units','inches',...
    'position',[1 1 4.25 2.5],'resize','off','paperpositionmode','auto');
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
hold on;
for ireg = 1:length(regnames)
    m_text(lon1(lonidxs(ireg)),lat1(latidxs(ireg)),regnames(ireg),'fontsize',8,'fontweight','bold');
end
m_pcolor(lon2,lat2,nan(size(reg_map)));
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
colormap(seqcmap); shading flat;
nreg = max(reg_map(:));
for i = 1:nreg
    RR = double(reg_map==i);
    RR(M3d(:,:,1)==0) = NaN;
    m_contour(lon2,lat2,RR,[.5 .5],'k',...
        'linewidth',reglinewidth,'linecolor',reglinecolor);
end
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('Example time series locations in various ocean regions',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

print(f, [fig_save_path 'suppfig4_1of2_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'suppfig4_1of2_PSRfbpaper_final.png'], '-dpng', '-r300');

%-------------------------------
% Plot figures
%-------------------------------
months = datetime(kost.Xitimemonames,'ConvertFrom','datenum');
qual9cmap = cbrewer('qual','Set1',9);
qual9bgcmap = ones(size(qual9cmap));
qual9bgcmap(6,:) = [0.8 0.8 0.8];
titlefontsize = 12; titlefontwt = 'bold';
sgtitlefontsize = 14; sgtitlefontwt = 'bold';

ieratio = 2;
inpp = 1;

% for ieratio = 1:eration
%     for inpp = 1:nppn
f=figure;
set(f,'color','white','units','inches','position',[0.5 0.5 12 4.67],'resize','off');
isubplot =  1;
for ireg = 1:length(regnames)
    ax=subplot(2,3,ireg);
    yyaxis left;
    plot(months,...
        squeeze(cexp_mon_maps(latidxs(ireg),lonidxs(ireg),:,inpp,ieratio)...
        ./nanmean(cexp_mon_maps(latidxs(ireg),lonidxs(ireg),:,inpp,ieratio))));
    yyaxis right;
    plot(months, squeeze(beta_mon_map(latidxs(ireg),lonidxs(ireg),:)));
    xlim([datetime(1999,1,1), datetime(2005,1,1)]);
    xtickformat('MMM-yy'); xtickangle(45); grid on;
    ax.XAxis.FontSize = 8;
    title(regnames{ireg},'fontsize',titlefontsize,'fontweight',titlefontwt);
    isubplot = isubplot + 1;
end
sgtitle(['\beta_ and ' npp_algs{inpp} '+' eratio_algs_sh{ieratio} ' export time series'],...
    'fontsize',sgtitlefontsize,'fontweight',sgtitlefontwt);

print(f, [fig_save_path 'suppfig4_2of2_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'suppfig4_2of2_PSRfbpaper_final.png'], '-dpng', '-r300');
%     end
% end
