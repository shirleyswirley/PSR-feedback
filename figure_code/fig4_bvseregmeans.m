clear all;
close all;
setup_figs;

%-------------------------------
% Load grid variables + regions map
%-------------------------------
reg_struct = load([data_path ...
    'PRiSM_regions_2deg.mat']);
M3d = reg_struct.M3d; gridd = reg_struct.gridd;
lon2 = gridd.xt; lat2 = gridd.yt;
isurf = find(M3d(:,:,1)==1);
area_ocn_only = nan(size(gridd.Areat)); % m^2
area_ocn_only(isurf) = gridd.Areat(isurf);
orig_reg_map = reg_struct.R2d;
% Original region numbers:
% 1 = AAZ, 2 = SAZ, 3 = STA, 4 = STP,
% 5 = ETA, 6 = ETP, 7 = NA, 8 = NP,
% 9 = IND, 10 = ARC, 11 = MED
reg_map = nan(size(orig_reg_map));
reg_map(orig_reg_map==1) = 1;
reg_map(orig_reg_map==2) = 2;
reg_map(orig_reg_map==3) = 3;
reg_map(orig_reg_map==4) = 4;
reg_map(orig_reg_map==9) = 5;
reg_map(orig_reg_map==5) = 6;
reg_map(orig_reg_map==6) = 7;
reg_map(orig_reg_map==7) = 8;
reg_map(orig_reg_map==8) = 9;
reg_map(orig_reg_map==10) = 10;
reg_map(orig_reg_map==11) = 11;
% Reordered/abridged region numbers:
% 1 = AAZ, 2 = SAZ, 3 = STA, 4 = STP,
% 5 = IND, 6 = ETA, 7 = ETP, 8 = NA,
% 9 = NP, 10 = ARC, 11 = MED
reg_names = [reg_struct.regnamesabbrev(1:4) ...
    reg_struct.regnamesabbrev(9) ...
    reg_struct.regnamesabbrev(5:8) ...
    reg_struct.regnamesabbrev(10:11)];
nregs = length(reg_names);
%figure; pcolor(reg_map); shading flat; colorbar;

% - Define the eratio and npp algs you want to use
eratio_algs = {'Laws2000','D2005PP','Laws2011D'};
eratio_algs_sh = {'L2000','D2005','L2011'};
npp_algs = {'VGPM','VGPME','CbPM'};
eration = length(eratio_algs);
nppn = length(npp_algs);
timen = 160; % num of months in time series

%-------------------------------
% Load and regrid beta vs. export maps
%-------------------------------
load([data_path 'globalbetavsexportmaps_tempregress/betavsexportn1tempregmo_VGPMVGPMECbPMcompiled.mat'],'lat1','lon1')
lon1 = [lon1(181:end)' lon1(1:180)'+360];

% - Load b vs e maps + populate bvse_maps
bvse_maps = nan(length(lat1),length(lon1),eration,nppn);
for ieratio = 1:eration
    for inpp = 1:nppn
        load([data_path 'globalbetavsexportmaps_tempregress/' ...
            npp_algs{inpp} 'based/Kostbetavs' ...
            eratio_algs{ieratio} 'exportn1tempregmo.mat'],...
            'betavsetempregresscoeffmapzeros');
        bvse_maps(:,:,ieratio,inpp) = betavsetempregresscoeffmapzeros;
    end
end

% - Put b vs. e maps on 2-deg PRiSM grid
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
% Weight separate export maps using regional weights
% from Weber et al. (2016) supplementary table 2
%-------------------------------
% - Recreate weights table from
% Weber et al. (2016) supplementary table 2
nonwtedwt=1/(eration*nppn); %1/9 
regwtstable = nan(nregs,eration,nppn);
regwtstable(:,1,1) = [0.1139 0.3207 0.2308 0.0504 nonwtedwt 0.0656 0.0656 0.0478 0      nonwtedwt nonwtedwt];
regwtstable(:,2,1) = [0.1508 0.2328 0.1677 0.0300 nonwtedwt 0.0729 0.0729 0.0697 0.0026 nonwtedwt nonwtedwt];
regwtstable(:,3,1) = [0.0927 0.0454 0.0975 0.0208 nonwtedwt 0.0445 0.0445 0.1169 0.1855 nonwtedwt nonwtedwt];
regwtstable(:,1,2) = [0.1507 0.0420 0.1419 0.0663 nonwtedwt 0.1213 0.1213 0.1184 0.1197 nonwtedwt nonwtedwt];
regwtstable(:,2,2) = [0.1349 0.0212 0.0993 0.0435 nonwtedwt 0.1516 0.1516 0.1294 0.2379 nonwtedwt nonwtedwt];
regwtstable(:,3,2) = [0.0622 0.0036 0.0636 0.0292 nonwtedwt 0.1080 0.1080 0.1308 0.1211 nonwtedwt nonwtedwt];
regwtstable(:,1,3) = [0.0478 0.2014 0.0900 0.2688 nonwtedwt 0.1667 0.1667 0.1263 0.0107 nonwtedwt nonwtedwt];
regwtstable(:,2,3) = [0.1215 0.1141 0.0640 0.2695 nonwtedwt 0.1047 0.1047 0.1322 0.0978 nonwtedwt nonwtedwt];
regwtstable(:,3,3) = [0.1255 0.0188 0.0451 0.2216 nonwtedwt 0.1648 0.1648 0.1286 0.2247 nonwtedwt nonwtedwt];

% - Construct weight maps
regwtsmap = nan(length(lat2),length(lon2),eration,nppn);
for inpp = 1:nppn
    for ieratio = 1:eration
        regwtsmaptemp = nan(size(reg_map));
        for ireg = 1:nregs
            regwtsmaptemp(reg_map==ireg) = regwtstable(ireg,ieratio,inpp); 
        end
        regwtsmap(:,:,ieratio,inpp) = regwtsmaptemp;
    end
end

%-------------------------------
% Calculate wted reg mean b vs. e
%-------------------------------
regwtedmeanbvsemap = nan(length(lat2),length(lon2));
regwtedmean1stdbvsemap = nan(length(lat2),length(lon2));
regwtedmeanbvsermat = nan(nregs,eration,nppn);
regwtedmeanbvser = nan(nregs,1);
regwted1stdbvser = nan(nregs,1);
bvseregblocksmap = nan(length(lat2),length(lon2));

% - Compute weighted mean
% global resolved beta vs. export slope map
regwtedmeanbvsemap = sum(sum(regwtsmap.*bvse_maps2,4),3);

% - Compute weighted stdev above/below mean
% global resolved beta vs. export slope map
regwted1stdbvsemap = zeros(length(lat2),length(lon2));
for inpp = 1:nppn
    for ieratio = 1:eration
        regwted1stdbvsemap = regwted1stdbvsemap + regwtsmap(:,:,ieratio,inpp).*(bvse_maps2(:,:,ieratio,inpp)-regwtedmeanbvsemap).^2;
    end
end
regwted1stdbvsemap = regwted1stdbvsemap.^(1/2);

% - Compute regional mean temporally regressed
% beta vs. export slope for each NPP/e-ratio combo 
for ieratio = 1:eration
    for inpp = 1:nppn
        %needs mapnow,nregs,reg_map,area_ocn_only
        mapnow = bvse_maps2(:,:,ieratio,inpp); calc_reg_means;
        regwtedmeanbvsermat(:,ieratio,inpp) = mapnowrm'; 
    end
end
regwtedmeanbvser = sum(sum(regwtedmeanbvsermat.*regwtstable,2),3);

% - Compute regional mean temporally regressed
% beta vs. export slope stdevs for each NPP/e-ratio combo 
for iregion = 1:nregs
    regwted1stdbvser(iregion) = sum(sum( ...
        squeeze(regwtstable(iregion,:,:)).* ...
        (squeeze(regwtedmeanbvsermat(iregion,:,:))-regwtedmeanbvser(iregion)).^2));
end
regwted1stdbvser = regwted1stdbvser.^(1/2);

% - Create final regional block
% beta vs. export slope maps 
bvseregblocksmap = nan(size(reg_map));
for iregion = 1:nregs
    bvseregblocksmap(reg_map==iregion) = regwtedmeanbvser(iregion);  
end

%-------------------------------
% Plot figure
%-------------------------------
regnumsnow = 1:9;
divcmap = flipud(cbrewer('div','RdYlBu',21,'linear'));
qual9 = cbrewer('qual','Set1',9);
mapproj = 'gall-peters'; landcolor = [0.6 0.6 0.6];
reglinewidth = 2; reglinecolor = [0 0 0];
labelfontsize = 11; labelfontwt = 'normal';
titlefontsize = 13; titlefontwt = 'bold';
cbticklen = 0.03;

f=figure;
set(f,'color','white',...
    'units','inches','position',[0.5 0.5 6.5 13]);

% - Reg mean beta vs. norm export slopes
% for all 9 combos of NPP + e-ratio
ax = subplot(3,1,1);
regwtedmeanbvsermattemp1 = squeeze(regwtedmeanbvsermat);
regwtedmeanbvsermattemp = ...
    [regwtedmeanbvsermattemp1(regnumsnow,:,1) ...
     regwtedmeanbvsermattemp1(regnumsnow,:,2) ...
     regwtedmeanbvsermattemp1(regnumsnow,:,3)];
bh = bar(1:length(reg_names(regnumsnow)),...
     regwtedmeanbvsermattemp(regnumsnow,:),1,...
     'FaceColor','flat','EdgeColor',[0.35 0.35 0.35]);
for k = 1:9
    bh(k).CData = k-1;
end
exportlabels = {[npp_algs{1} '+' eratio_algs_sh{1}],...
    [npp_algs{2} '+' eratio_algs_sh{1}],[npp_algs{3} '+' eratio_algs_sh{1}],...
    [npp_algs{1} '+' eratio_algs_sh{2}],[npp_algs{2} '+' eratio_algs_sh{2}],...
    [npp_algs{3} '+' eratio_algs_sh{2}],[npp_algs{1} '+' eratio_algs_sh{3}],...
    [npp_algs{2} '+' eratio_algs_sh{3}],[npp_algs{3} '+' eratio_algs_sh{3}]};
cb = colorbar; cb.FontSize = 8;
cb.Ticks = 0.5:length(exportlabels)+0.5;
cb.TickLabels = exportlabels; cb.TickLength = cbticklen;
colormap(ax, qual9); caxis([0 9]);
xlim([0.5 length(reg_names(regnumsnow))+0.5]);
set(gca,'XTickLabel',reg_names(regnumsnow),...
    'fontsize',labelfontsize,'fontweight',labelfontwt);
ylim([-0.65 0.2]); ylabel('\beta vs. export slope');
%title('All regional mean $\frac{d\beta_{obs}}{dE_{n,obs}}$','Interpreter','latex');
title('All regional mean_   ',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

% - Reg mean beta vs. norm export slopes
% bar chart w/ stdev error bars
regwtedmeanbvsernow = regwtedmeanbvser(regnumsnow);
ax = subplot(3,1,2);
bh = bar(1:length(reg_names(regnumsnow)),...
     [regwtedmeanbvsernow],'FaceColor','flat');hold on;
cmapnow = interp1(linspace(-max(abs(regwtedmeanbvsernow)),...
    max(abs(regwtedmeanbvsernow)),size(divcmap,1)),...
    divcmap,regwtedmeanbvsernow);
for k = 1:length(regnumsnow)
    bh.CData(k,:) = cmapnow(k,:);
end
errorbar(1:length(reg_names(regnumsnow)),[regwtedmeanbvsernow],...
    [regwted1stdbvser(regnumsnow)],'k.');
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
colormap(ax,divcmap);
caxis([-max(abs(regwtedmeanbvser)) max(abs(regwtedmeanbvser))]);
xlim([0.5 length(reg_names(regnumsnow))+0.5]);
set(gca,'XTickLabel',reg_names(regnumsnow),...
    'fontsize',labelfontsize,'fontweight',labelfontwt);
ylim([-0.6 0.2]); ylabel('\beta vs. export slope');
%title('Regionally weighted, regional mean $\frac{d\beta_{obs}}{dE_{n,obs}}$','Interpreter','latex');
title('Regionally weighted, regional mean_   ',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

% - Final regional blocks global map
ax = subplot(3,1,3);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_pcolor(lon2,lat2,bvseregblocksmap);
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
colormap(ax,divcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
caxis([-max(max(abs(regwtedmeanbvser))) max(max(abs(regwtedmeanbvser)))]);
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
hold on;
nreg = max(reg_map(:));
for i = 1:nreg
    RR = double(reg_map==i);
    RR(M3d(:,:,1)==0) = NaN;
    m_contour(lon2,lat2,RR,[.5 .5],'k',...
        'linewidth',reglinewidth,'linecolor',reglinecolor);
end
%title('Final map of $\frac{d\beta_{obs}}{dE_{n,obs}}$','Interpreter','latex');
title('Final map of_   ',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

print(f, [fig_save_path 'fig4_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'fig4_PSRfbpaper_final.png'], '-dpng', '-r300');
