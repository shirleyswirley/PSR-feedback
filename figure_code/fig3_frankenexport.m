close all;
clear all;

%-------------------------------
% Set up paths
%-------------------------------
utils_path = '/graid1/shirlleu/PSRfeedback/utils/';
data_path = '/graid1/shirlleu/PSRfeedback/data/';
fig_save_path = '/graid1/shirlleu/PSRfeedback/pdfs_pngs/';

%-------------------------------
% Set up warnings and utils
%-------------------------------
warning('off','all');
addpath(genpath(utils_path));

%-------------------------------
% Load all export and beta maps + grid variables
%-------------------------------

% - Define the eratio and npp algs you want to use
npp_algs = {'VGPM','VGPME','CbPM'};
eratio_algs = {'Laws2000','D2005PP','Laws2011D'};
eratio_algs_sh = {'L2000','D2005','L2011'};

% - Initialize matrices to hold all monthly
% export maps (cexp_mon_maps)
eration = length(eratio_algs);
nppn = length(npp_algs);
timen = 160; % num of months in time series
load([data_path_bvse 'betavsexportn1tempregmo_VGPMVGPMECbPMcompiled.mat'],...
    'lat1','lon1')
lon1 = [lon1(181:end)' lon1(1:180)'+360];
cexp_mon_maps = nan(length(lat1),length(lon1),timen,nppn,eration);

% - Load b vs e maps + populate bvse_maps
for ieratio = 1:eration
    for inpp = 1:nppn
        load([data_path_bvse npp_algs{inpp} 'based/Kostbetavs' ...
            eratio_algs{ieratio} 'exportn1tempregmo.mat'],...
            'betavsetempregresscoeffmapzeros');
        bvse_maps(:,:,ieratio,inpp) = betavsetempregresscoeffmapzeros;
    end
end

% - Load export maps + populate cexp_mon_maps 
for inpp = 1:nppn
    for ieratio = 1:eration
        load([data_path_cexp 'expCmo_' eratio_algs{ieratio} 'from' ...
            npp_algs{inpp} '.mat'],'expCmo1');
        disp([data_path_cexp 'expCmo_' eratio_algs{ieratio} 'from' ...
            npp_algs{inpp} '.mat']);
        cexp_mon_maps(:,:,:,inpp,ieratio) = ...
            [expCmo1(:,181:end,:) expCmo1(:,1:180,:)];
    end
end

% - Initialize matrix to hold all monthly
% beta maps (beta_mon_maps)
beta_mon_maps = nan(length(lat1),length(lon1),timen);
 
% - Load beta maps + populate beta_mon_maps
kost = load([data_path 'Kost_beta_monthly_1deg.mat'],'Ximo1','lon1','lat1');
kost.lon1 = [kost.lon1(181:end)' kost.lon1(1:180)'+360];
beta_mon_maps1 = [kost.Ximo1(:,181:end,:) kost.Ximo1(:,1:180,:)];
beta_mon_maps = nan(size(beta_mon_maps1));
for itime = 1:timen
    beta_mon_maps(:,:,ieratio,inpp) = ...
        interp2(kost.lon1,kost.lat1,beta_mon_maps1(:,:,itime),lon1,lat1);
end

%-------------------------------
% Weight separate NPP+e-ratio=export combos
% using Weber et al. (2016)'s regional weights
%-------------------------------

% - Load frankenstein C export map + export regions map
cexp_franken = load([data_path 'franken_export_1deg.mat']);
fgrid = cexp_franken.grid; reg_map = cexp_franken.R2d;
reg_map = interp2([fgrid.xt(end)-360 fgrid.xt],...
    fgrid.yt',[reg_map(:,end) reg_map],lon1,lat1,'nearest');
reg_names = cexp_franken.reg_names;
reg_names = [reg_names,'Indian','Arctic','MedSeaPlus'];
reg_names_sh = {'AAZ','SAZ','STA','STP','ETA','ETP','NA','NP','IND','ARC','MED'};
[lon1grid,lat1grid] = meshgrid(lon1,lat1);
reg_map(reg_map==0 & lat1grid<30) = 9; % Indian Ocean region
reg_map(reg_map==0 & lat1grid>60) = 10; % Arctic Ocean region
reg_map(reg_map==0) = 11; % Mediterannean Sea + other seas

% - Create regional weights table
% (from Weber et al. 2016 supplementary table 2)
nonwtedwt=1/(eration*nppn); %1/9 
regwts_table = nan(length(reg_names),eration,nppn);
regwts_table(:,1,1) = [0.1139 0.3207 0.2308 0.0504 0.0656 0.0656 0.0478 0      nonwtedwt nonwtedwt nonwtedwt];
regwts_table(:,2,1) = [0.1508 0.2328 0.1677 0.0300 0.0729 0.0729 0.0697 0.0026 nonwtedwt nonwtedwt nonwtedwt];
regwts_table(:,3,1) = [0.0927 0.0454 0.0975 0.0208 0.0445 0.0445 0.1169 0.1855 nonwtedwt nonwtedwt nonwtedwt];
regwts_table(:,1,2) = [0.1507 0.0420 0.1419 0.0663 0.1213 0.1213 0.1184 0.1197 nonwtedwt nonwtedwt nonwtedwt];
regwts_table(:,2,2) = [0.1349 0.0212 0.0993 0.0435 0.1516 0.1516 0.1294 0.2379 nonwtedwt nonwtedwt nonwtedwt];
regwts_table(:,3,2) = [0.0622 0.0036 0.0636 0.0292 0.1080 0.1080 0.1308 0.1211 nonwtedwt nonwtedwt nonwtedwt];
regwts_table(:,1,3) = [0.0478 0.2014 0.0900 0.2688 0.1667 0.1667 0.1263 0.0107 nonwtedwt nonwtedwt nonwtedwt];
regwts_table(:,2,3) = [0.1215 0.1141 0.0640 0.2695 0.1047 0.1047 0.1322 0.0978 nonwtedwt nonwtedwt nonwtedwt];
regwts_table(:,3,3) = [0.1255 0.0188 0.0451 0.2216 0.1648 0.1648 0.1286 0.2247 nonwtedwt nonwtedwt nonwtedwt];

% - Construct regional weights maps
regwts_map = nan(length(lat1),length(lon1),eration,nppn);
for inpp = 1:nppn
    for ieratio = 1:eration
        regwts_maptemp = nan(size(reg_map));
        for ireg = 1:length(reg_names)
            regwts_maptemp(reg_map==ireg) = regwts_table(ireg,ieratio,inpp);
        end
        regwts_map(:,:,ieratio,inpp) = regwts_maptemp;
    end
end

% - Check to make sure weight maps were constructed properly
plotnow=0;
if plotnow==1
    f=figure;
    set(f,'color','white','units','inches',...
        'position',[0.5 0.5 11 6],'resize','off');
    isubplot =  1; % subplot idx
    for inpp = 1:nppn
        for ieratio = 1:eration
            subplot(3,3,isubplot);
            pcolor(lon1,lat1,regwts_map(:,:,ieratio,inpp));
            shading flat;colorbar;
            caxis([0 max(max(max(max(regwts_table))))]);
            title(['All beta vs.' 10 npp_algs{inpp} ...
                ' npp+' eratio_algs{ieratio} ' e-ratio reg wts']);
            isubplot = isubplot+1;
        end
    end
end

%--------------------------------------------------
% Create franken export monthly time series
%--------------------------------------------------
cexpf_mon_maps = zeros(length(lat1),length(lon1),timen); 
for itime = 1:timen 
    for inpp = 1:nppn
        for ieratio = 1:eration
            cexpf_mon_maps(:,:,itime) = ...
                cexpf_mon_maps(:,:,itime) + ...
                cexp_mon_maps(:,:,itime,inpp,ieratio).*regwts_map(:,:,ieratio,inpp);
        end
    end
end

cmap=cbrewer('seq','YlGnBu',100,'linear');
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
linewidth=2; % for regional contour lines
linecolor=[0 0 0]; % for regional contour lines
labelfontsize = 14;

f=figure;
set(f,'color','white','units','inches',...
    'position',[0.5 0.5 3.5 2],'resize','off');

plotregcontours=1;
expCfrank = interp2([fgrid.xt(end)-360 fgrid.xt],fgrid.yt,[franken.Cfrank(:,end) franken.Cfrank],lon1,lat1);
% Next few lines are to fill in the small nan islands east of Australia
latn = -10; lats = -30; lonl = 160; lonr = 184;
latidxnan=find(lat1<latn&lat1>lats);
lonidxnan=find(lon1<lonr&lon1>lonl);
expCfrank1 = inpaint_nans(expCfrank(latidxnan,lonidxnan)); % interpolate over small nan islands
expCfrank(latidxnan,lonidxnan)=expCfrank1; % fill in the small nan islands
%figure; pcolor(flipud(expCfrank)); shading flat; colorbar;
% Next few lines are to fill in the small nan islands south east of Australia
latn = -35; lats = -48; lonl = 160; lonr = 184;
latidxnan=find(lat1<latn&lat1>lats);
lonidxnan=find(lon1<lonr&lon1>lonl);
expCfrank1 = inpaint_nans(expCfrank(latidxnan,lonidxnan)); % interpolate over small nan islands
expCfrank(latidxnan,lonidxnan)=expCfrank1; % fill in the small nan islands
%figure; pcolor(flipud(expCfrank)); shading flat; colorbar;
% Next few lines are to fill in the small nan islands at the NA/STA boundary
latn = 40; lats = 35; lonl = 333; lonr = 337;
latidxnan=find(lat1<latn&lat1>lats);
lonidxnan=find(lon1<lonr&lon1>lonl);
expCfrank1 = inpaint_nans(expCfrank(latidxnan,lonidxnan)); % interpolate over small nan islands
expCfrank(latidxnan,lonidxnan)=expCfrank1; % fill in the small nan islands
%figure; pcolor(flipud(expCfrank)); shading flat; colorbar;
% Next few lines are to fill in the small nan islands in the Trop Pac
latn = 23; lats = 18; lonl = 200; lonr = 210;
latidxnan=find(lat1<latn&lat1>lats);
lonidxnan=find(lon1<lonr&lon1>lonl);
expCfrank1 = inpaint_nans(expCfrank(latidxnan,lonidxnan)); % interpolate over small nan islands
expCfrank(latidxnan,lonidxnan)=expCfrank1; % fill in the small nan islands
%figure; pcolor(flipud(expCfrank)); shading flat; colorbar;
subplot(122);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat1) max(lat1)])
m_contourf(lon1,lat1,log10(expCfrank),30);
%m_pcolor(lon1,lat1,log10(expCfrank));
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([-0.5 1]);
%caxis([0 8]);
colormap(flipud(cmap));
shading flat;
[arrb,arrt]=computepointyColorbarlims(log10(expCfrank),-0.5,1);
pointyColorbar(arrb,arrt);
h=pointyColorbar(arrb,arrt);
set(h,'fontsize',labelfontsize);
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
if plotregcontours==1
  hold on;
  %nreg = max(reg_map(:))
  nreg = 8;
  for i = 1:nreg
      RR = double(reg_map==i);
      RR(isnan(reg_map)) = NaN;
      m_contour(lon1,lat1,RR,[.5 .5],'k','linewidth',linewidth,'linecolor',linecolor);hold on;
  end
end
hold on;
tropex = zeros(size(reg_map));
tropex(latidxt,lonidxt) = 1;
tropex(isnan(reg_map)) = NaN;
m_contour(lon1,lat1,tropex,[.5 .5],'linewidth',linewidth,'linecolor',[1 1 1]);hold on;
hold on;
subtropex = zeros(size(reg_map));
subtropex(latidxst,lonidxst) = 1;
subtropex(isnan(reg_map)) = NaN;
m_contour(lon1,lat1,subtropex,[.5 .5],'linewidth',linewidth,'linecolor',[1 1 1]);hold on;
title('franken export');
