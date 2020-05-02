clear all;
close all;
setup_figs;

%-------------------------------
% Load grid variables + regions map
%-------------------------------
% - Load grid variables + regions map
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
erationamearray = {'Laws2000','D2005PP','Laws2011D'};
eratioshortnamearray = {'L2000','D2005','L2011'};
NPPnamearray = {'VGPM','VGPME','CbPM'};
eration = length(erationamearray);
nppn = length(NPPnamearray);
timen = 160; % num of months in time series

%--------------------------------------------------
% Load and regrid beta vs. export
%--------------------------------------------------
% Remember that this grid is ever so slightly diff from the Kost and Guidi beta grids,
% but it's totally ok b/c both are just lower res from orig sat data + close enough
load('/graid1/shirlleu/NPPslopeproj/matfiles/globalexportmaps_spatvareratioandezdepth/expCann_D2005PPfromphysatchl1.mat','lat1','lon1');

% - Temporally regressed beta vs. normalized export maps
% Initialize beta vs. export filenames
fnamebase = ['/graid1/shirlleu/NPPslopeproj/matfiles/globalbetavsexporttempregressmaps/'];
% Initialize matrix to hold all beta vs. export maps
bvse = nan(length(lat1),length(lon1),eration,nppn);
% Load beta vs. export maps
for ieratio = 1:eration
    for inpp = 1:nppn
        load([fnamebase NPPnamearray{inpp} 'based/Kostbetavs' erationamearray{ieratio} 'exportn1tempregmo.mat'],'betavsetempregresscoeffmapzeros');
        bvse(:,:,ieratio,inpp) = betavsetempregresscoeffmapzeros;
    end
end

% - Put b vs. e on PRiSM grid so we can use those approximate areas to calculate area-wted reg means
lon1p = [lon1(181:end)' lon1(1:180)'+360];
[lon1pgrid,lat1grid] = meshgrid(lon1p,lat1);
bvsepgrid1 = [bvse(:,181:end,:,:,:) bvse(:,1:180,:,:,:)];
bvsepgrid = nan(length(lat2),length(lon2),eration,nppn);
for ieratio = 1:eration
    for inpp = 1:nppn
        bvsetemp = bvsepgrid1(:,:,ieratio,inpp);
        bvsepgrid(:,:,ieratio,inpp) = interp2(lon1pgrid,lat1grid,bvsetemp,gridd.XT,gridd.YT);
    end
end

%--------------------------------------------------
% Weight separate CbPM, VGPM, and VGPM-Eppley maps using
% Tom's regional weights from the PNAS transfer efficiency paper 
%--------------------------------------------------
% - Recreate weights table
% from Weber et al. (2016) supplementary table 2
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

%--------------------------------------------------
% CALCULATE OBS WEIGHTED MEAN BETA VS. EXPORT GLOBAL MAP
% + OBS WEIGHTED REGIONAL MEANS USING 
% METHOD 1: TEMPORALLY REGRESS BY GRID PT AND
% THEN SPATIALLY AVERAGE REGIONAL BETA VS. EXPORT SLOPES 
%--------------------------------------------------
regwtedmeanbvsemap = nan(length(lat2),length(lon2));
regwtedmean1stdbvsemap = nan(length(lat2),length(lon2));
regwtedmeanbvsermat = nan(nregs,eration,nppn);
regwtedmeanbvser = nan(nregs,1);
regwted1stdbvser = nan(nregs,1);
bvseregblocksmap = nan(length(lat2),length(lon2));

% - Compute weighted mean beta vs. export slope map
regwtedmeanbvsemap = sum(sum(regwtsmap.*bvsepgrid,4),3);

% - Compute weighted stdev above/below mean beta vs. export slope map
regwted1stdbvsemap = zeros(length(lat2),length(lon2));
for inpp = 1:nppn
    for ieratio = 1:eration
        regwted1stdbvsemap = regwted1stdbvsemap + regwtsmap(:,:,ieratio,inpp).*(bvsepgrid(:,:,ieratio,inpp)-regwtedmeanbvsemap).^2;
    end
end
regwted1stdbvsemap = regwted1stdbvsemap.^(1/2);

% - Compute regional mean temporally regressed beta vs. export slope for each NPP/e-ratio/beta combo 
for ieratio = 1:eration
    for inpp = 1:nppn
        bvsetemp = bvsepgrid(:,:,ieratio,inpp);
        %needs mapnow,nregs,reg_map,area_ocn_only
        mapnow = bvsetemp; calc_reg_means;
        regwtedmeanbvsermat(:,ieratio,inpp) = mapnowrm'; 
    end
end
regwtedmeanbvser = sum(sum(regwtedmeanbvsermat.*regwtstable,2),3);

% - Compute regional mean temporally regressed beta vs. export slope stdevs for each NPP/e-ratio/beta combo 
for iregion = 1:nregs
    regwted1stdbvser(iregion) = sum(sum( ...
        squeeze(regwtstable(iregion,:,:)).* ...
        (squeeze(regwtedmeanbvsermat(iregion,:,:))-regwtedmeanbvser(iregion)).^2));
end
regwted1stdbvser = regwted1stdbvser.^(1/2);

% - Create final regional block beta vs. export slope maps 
bvseregblocksmaptemp = nan(size(reg_map));
for iregion = 1:nregs
    bvseregblocksmaptemp(reg_map==iregion) = regwtedmeanbvser(iregion);  
end
bvseregblocksmap = bvseregblocksmaptemp;

%--------------------------------------------------
% FIGURE 4 AND SUPPLEMENTARY FIGURE 5
% Plot beta vs. export regionally-weighted global map, all beta vs. export
% regional means, regionally-weighted beta vs. export regional means,
% final global map of beta vs. export regionally-weighted regional means
%--------------------------------------------------
redblue = flipud(cbrewer('div','RdYlBu',100,'linear'));
lon2shift = [lon2(91:180)-360 lon2(1:90)];
linewidth=2; % for regional contour lines
linecolor=[0 0 0]; % for regional contour lines
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
labelfontsize = 14;

% - Plot global maps of b vs. e with chosen significance hatching
% and w/ all 9 regional mean b vs. e values,
% w/ Kost and Guidi in diff figs with 2 by 2 subplots
qual9 = cbrewer('qual','Set1',9);

f=figure;
set(f,'color','white','units','inches','position',[0.5 0.5 13 7]);

% - Reg wted mean beta vs. norm export slope map
cmax = 1.06;
% Kost min, max = -1.0541, 0.2229
% Guidi min, max = -0.9772, 0.9740
ax=subplot(2,2,1);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,regwtedmeanbvsemap,30);
%m_pcolor(lon2,lat2,regwtedmeanbvsemap(:,:,ibeta));
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([-cmax cmax]);
colormap(ax,redblue);
shading flat;
h = colorbar;
set(h,'fontsize',labelfontsize);
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
hold on;

% - Reg mean beta vs. norm export slopes for all 9 combos of NPP + e-ratio
ax=subplot(2,2,2);
regwtedmeanbvsermattemp1 = squeeze(regwtedmeanbvsermat);
% end-2 is to exclude Med and Arc region from the final plot
regwtedmeanbvsermattemp = [regwtedmeanbvsermattemp1(1:end-2,:,1) regwtedmeanbvsermattemp1(1:end-2,:,2) regwtedmeanbvsermattemp1(1:end-2,:,3)];
% Reorder Indian Ocean
%regwtedmeanbvsermattemp = [regwtedmeanbvsermattemp(1:4,:); regwtedmeanbvsermattemp(9,:); regwtedmeanbvsermattemp(5:8,:)]; 
bh=bar(1:length(reg_names(1:end-2)),regwtedmeanbvsermattemp(1:9,:),1);
set(bh,'edgecolor',[0.35 0.35 0.35]);
colormap(ax,qual9);
exportlabels = {[NPPnamearray{1} '+' eratioshortnamearray{1}],[NPPnamearray{2} '+' eratioshortnamearray{1}],[NPPnamearray{3} '+' eratioshortnamearray{1}],...
    [NPPnamearray{1} '+' eratioshortnamearray{2}],[NPPnamearray{2} '+' eratioshortnamearray{2}],[NPPnamearray{3} '+' eratioshortnamearray{2}],...
    [NPPnamearray{1} '+' eratioshortnamearray{3}],[NPPnamearray{2} '+' eratioshortnamearray{3}],[NPPnamearray{3} '+' eratioshortnamearray{3}]};
caxis([0 9]);
%colorbar;
%h=colorbar('YTick',0.5:length(exportlabels)+0.5,'YTickLabels',exportlabels);
%cbfreeze(colorbar('YTickLabels',exportlabels));
%set(h,'fontsize',12,'fontweight','bold');
xlim([0.5 length(reg_names(1:end-2))+0.5]);
set(gca,'XTickLabel',reg_names(1:end-2),'fontsize',labelfontsize-2,'fontweight','bold');
ylim([-0.65 0.2]);
ylabel('\beta vs. export slope');
title('Reg mean \beta vs. norm export slope');

% - Reg mean beta vs. norm export slopes bar chart 
% end-2 is to exclude Med and Arc regions from the final plot
% Reorder Indian Ocean
ax=subplot(2,2,3);
bh=bar(1:length(reg_names(1:end-2)),...
   [regwtedmeanbvser(1:9)]);hold on;
   %[regwtedmeanbvser(1:4,ibeta); regwtedmeanbvser(9,ibeta); regwtedmeanbvser(5:8,ibeta)]);hold on;
bc=get(bh,'Children');
%set(bc,'CData',[regwtedmeanbvser(1:4,ibeta); regwtedmeanbvser(9,ibeta); regwtedmeanbvser(5:8,ibeta)]);
set(bc,'CData',[regwtedmeanbvser(1:9)]);
set(h,'fontsize',labelfontsize);
errorbar(1:length(reg_names(1:end-2)),[regwtedmeanbvser(1:9)],...
    [regwted1stdbvser(1:9)],'k.');
caxis([-max(max(abs(regwtedmeanbvser))) max(max(abs(regwtedmeanbvser)))]);
colormap(ax,redblue);
%errorbar(1:length(reg_names(1:end-2)),[regwtedmeanbvser(1:4,ibeta); regwtedmeanbvser(9,ibeta); regwtedmeanbvser(5:8,ibeta)],...
    %[regwted1stdbvser(1:4); regwted1stdbvser(9); regwted1stdbvser(5:8)],'k.');
xlim([0.5 length(reg_names(1:end-2))+0.5]);
set(gca,'XTickLabel', reg_names(1:end-2),'fontsize',labelfontsize-2,'fontweight','bold');
%set(gca, 'XTick', 1:length(reg_names(1:end-2)), 'XTickLabel', reg_names(1:end-2));
% % OR: 
%NumTicks = length(reg_names(1:end-2));
%L = get(gca,'XLim');
%set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%set(gca,'XTickLabel',[' ',reg_names(1:end-2),' ']);
ylim([-0.6 0.2]);
ylabel('\beta vs. export slope');
title('Reg mean \beta vs. norm export slope');
%grid on;

% - Final regional blocks global map
ax=subplot(2,2,4);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_pcolor(lon2,lat2,bvseregblocksmap);
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([-max(max(abs(regwtedmeanbvser))) max(max(abs(regwtedmeanbvser)))]);
colormap(ax,redblue);
shading flat;
set(h,'fontsize',labelfontsize);
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
hold on;
nreg = max(reg_map(:))
for i = 1:nreg
    RR = double(reg_map==i);
    RR(M3d(:,:,1)==0) = NaN;
    m_contour(lon2,lat2,RR,[.5 .5],'k','linewidth',linewidth,'linecolor',linecolor);hold on;
end
title('Final global \beta vs. norm export slope map');

