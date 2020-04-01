addpath(genpath('/graid1/shirlleu/matlabroutines'))

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define desired circ rate change
%%%%%%%%%%%%%%%%%%%%%%%%%%
%circfactor = 0.9;
circfactor = 1.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get relevant matfiles
%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------
% Tom's regions map (same grid as all other matfiles below)
%---------------------
regmapst = load('/graid1/shirlleu/PRiSM_GCM/data/PRiSMregions.mat');
regmapgrid = regmapst.gridd;
regmaporig = regmapst.R2d;
% 1 = AAZ, 2 = SAZ, 3 = STA, 4 = STP,
% 5 = ETA, 6 = ETP, 7 = NA, 8 = NP,
% 9 = IND, 10 = ARC, 11 = MED
regmap = nan(size(regmaporig));
regmap(regmaporig==1) = 1;
regmap(regmaporig==2) = 2;
regmap(regmaporig==3) = 3;
regmap(regmaporig==4) = 4;
regmap(regmaporig==9) = 5;
regmap(regmaporig==5) = 6;
regmap(regmaporig==6) = 7;
regmap(regmaporig==7) = 8;
regmap(regmaporig==8) = 9;
% 1 = AAZ, 2 = SAZ, 3 = STA, 4 = STP,
% 5 = IND, 6 = ETA, 7 = ETP, 8 = NA,
% 9 = NP
regnames = [regmapst.regnamesabbrev(1:4) regmapst.regnamesabbrev(9) regmapst.regnamesabbrev(5:8)];
%figure; pcolor(regmap); shading flat; colorbar;

%---------------------
% Baseline 1*circ
%---------------------
fnamebl = '/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/newtonmethod/spiralintoeq/prog/PSDfboff/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_Kostbetainit_PSDfboff_epsnolims_1.00circ_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_spiralintoeqiternum0.mat';

%---------------------
% PSD fb off
%---------------------
fnameoff = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/PSDfboff/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfboffKostbeta_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor);

%---------------------
% PSD fb on, all regions turned on
%---------------------
fnameon = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocksm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor);

%---------------------
% PSD fb on, all regions turned on,
% beta vs. export - 1stdev
%---------------------
fnameonms = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdbelowm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor);

%---------------------
% PSD fb on, all regions turned on,
% beta vs. export + 1stdev
%---------------------
fnameonps = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdabovem1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load relevant variables
%%%%%%%%%%%%%%%%%%%%%%%%%%

% - Common to all files
load(fnamebl,'output');
gridd = output.grid;
lats2 = gridd.yt;
lons2 = gridd.xt;
M3d = output.M3d;
isurf = find(M3d(:,:,1)==1);
area = gridd.Areat; % m^2
areaocnonly = nan(size(area));
areaocnonly(isurf) = area(isurf);

% - Load Kost beta map
kost2 = load('/graid1/shirlleu/PRiSM_GCM/sran/matfiles/PSD_beta_v2.mat');
beta_kost2 = mean(kost2.beta_clim,3)';
% -->the above goes w/ lats2, lons2
kost1 = load('/graid1/shirlleu/NPPslopeproj/matfiles/betamaps/Xi1mol.mat');
beta_kost1 = mean(kost1.Xiclim1,3)';
lons1 = kost1.lon1'; lats1 = flipud(kost1.lat1)';
lons1 = [lons1(181:end) lons1(1:180)+360];
beta_kost1 = flipud([beta_kost1(181:end,:); beta_kost1(1:180,:)]');
[beta_kost2_red,~,lons2_red,lats2_red] = RedGridRes(beta_kost1,2,2,lons1,lats1,0);
beta_kost2_red = inpaint_nans(beta_kost2_red);
beta_kost2_red = interp2(lons2_red,lats2_red',beta_kost2_red,lons2,lats2');

% - Load Kost beta remin depth map
load('remindepthmap_Kostbeta.mat','remindepthmap'); % fluxfrac = 1/e
rd_kost = remindepthmap;

% - Baseline 
load(fnamebl,'output');
ec_bl = output.expCmapSS; % don't need here tho

% - Feedback off, 1 degree
load(fnameoff,'output');
ec_off = output.expCmapnow(:,:,end); % don't need here tho

% - Feedback on
load(fnameon,'output');
ec_on = output.expCmapnow(:,:,end); % don't need here tho
beta_on = output.betamapnow(:,:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beta maps, 1 degree
% (regridded obs resolution,
% baseline + fboff const PSD map only)
% --> Conclusion is that 2-deg PSD_beta_v2.mat time mean used in PRiSM
% is very very close to 1-deg Xi time mean regridded to 2 deg,
% but not exactly so. It's good enough, but still dunno where
% PSD_beta_v2.mat came from.
%%%%%%%%%%%%%%%%%%%%%%%%%%
plotnow=0;
if plotnow==1
divnow = flipud(cbrewer('div','RdBu',21,'linear'));
seqnow = parula(20);
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
divplotseqnow = [0,1,0,1,0,1];
f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 8.5 6]);
ax=subplot(221);
m_proj(mapproj,'lon',[0 360],'lat',[min(lats1) max(lats1)]);
m_pcolor(lons1,lats1,beta_kost1);
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
c=colorbar; c.TickLength = 0.05; shading flat;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
caxis([3.5 6]);
colormap(ax,seqnow); hold on;
nreg = max(regmap(:));
for i = 1:nreg
    RR = double(regmap==i);
    RR(M3d(:,:,1)==0) = NaN;
    m_contour(lons2,lats2,RR,[.5 .5],'k','linewidth',1.5);
end
title('1 deg Kost beta, red from obs');
ax=subplot(222);
m_proj(mapproj,'lon',[0 360],'lat',[min(lats2) max(lats2)]);
m_pcolor(lons2,lats2,beta_kost2);
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
c=colorbar; c.TickLength = 0.05; shading flat;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
caxis([3.5 6]);
colormap(ax,seqnow); hold on;
nreg = max(regmap(:));
for i = 1:nreg
    RR = double(regmap==i);
    RR(M3d(:,:,1)==0) = NaN;
    m_contour(lons2,lats2,RR,[.5 .5],'k','linewidth',1.5);
end
title('2 deg Kost beta used in PRiSM');
ax=subplot(223);
m_proj(mapproj,'lon',[0 360],'lat',[min(lats2) max(lats2)]);
m_pcolor(lons2,lats2,beta_kost2_red);
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
c=colorbar; c.TickLength = 0.05; shading flat;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
caxis([3.5 6]);
colormap(ax,seqnow); hold on;
nreg = max(regmap(:));
for i = 1:nreg
    RR = double(regmap==i);
    RR(M3d(:,:,1)==0) = NaN;
    m_contour(lons2,lats2,RR,[.5 .5],'k','linewidth',1.5);
end
title('2 deg Kost beta, red from 1 deg, regrid to prism grid');
ax=subplot(224);
m_proj(mapproj,'lon',[0 360],'lat',[min(lats2) max(lats2)]);
m_pcolor(lons2,lats2,beta_kost2_red-beta_kost2);
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
c=colorbar; c.TickLength = 0.05; shading flat;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
caxis([-0.25 0.25]);
colormap(ax,divnow);
hold on;
nreg = max(regmap(:));
for i = 1:nreg
    RR = double(regmap==i);
    RR(M3d(:,:,1)==0) = NaN;
    m_contour(lons2,lats2,RR,[.5 .5],'k','linewidth',1.5);
end
title('2 deg Kost beta, red - prism');
saveas(f,'compare1and2degKostbeta.png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beta maps, 2 degree
% (model resolution)
%%%%%%%%%%%%%%%%%%%%%%%%%%
plotnow=1;
if plotnow==1
divnow = flipud(cbrewer('div','RdBu',21,'linear'));
seqnow = parula(20);
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
divplotseqnow = [0,1,0,1,0,1];
titlesnow = {'Baseline beta',[num2str(circfactor) ' circ on-off'],...
    [num2str(circfactor) ' circ off'],[num2str(circfactor) ' circ off-baseline'],...
    [num2str(circfactor) 'circ on'],[num2str(circfactor) 'circ on-baseline']};
plotsnow = nan(length(lats2),length(lons2),6);
plotsnow(:,:,1) = beta_kost2;
plotsnow(:,:,2) = beta_on-beta_kost2;
plotsnow(:,:,3) = beta_kost2;
plotsnow(:,:,4) = beta_kost2-beta_kost2;
plotsnow(:,:,5) = beta_on;
plotsnow(:,:,6) = beta_on-beta_kost2;
f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 8.5 7]);
for isp=1:6
    ax=subplot(3,2,isp);
    m_proj(mapproj,'lon',[0 360],'lat',[min(lats2) max(lats2)]);
    m_pcolor(lons2,lats2,plotsnow(:,:,isp));
    m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
    c=colorbar; c.TickLength = 0.05; shading flat;
    m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
    if divplotseqnow(isp)==1
        %caxis([-max(max(abs(plotsnow(:,:,isp)))) ...
        %        max(max(abs(plotsnow(:,:,isp))))]);
        caxis([-0.1 0.1]);
        colormap(ax,divnow);
    elseif divplotseqnow(isp)==0
        caxis([3.5 6]);
        colormap(ax,seqnow);
    end
    hold on;
    nreg = max(regmap(:));
    for i = 1:nreg
        RR = double(regmap==i);
        RR(M3d(:,:,1)==0) = NaN;
        m_contour(lons2,lats2,RR,[.5 .5],'k','linewidth',1.5);
    end
    title(titlesnow{isp});
end
saveas(f,sprintf('betaendmaps_%.2fcirc_narrowcaxis.png',circfactor));
%saveas(f,sprintf('betaendmaps_%.2fcirc_widecaxis.png',circfactor));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate global mean beta
%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('global mean beta, Kost:');
nansum(nansum(beta_kost2.*areaocnonly))./nansum(nansum(areaocnonly))
disp(['global mean beta, 100 yr fb on, ' num2str(circfactor) ' circ:']);
nansum(nansum(beta_on.*areaocnonly))./nansum(nansum(areaocnonly))

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate regional mean beta
%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_kost_regs = nan(1,length(regnames));
for ireg=1:length(regnames)
    kostbetanow = beta_kost2(regmap==ireg);
    areaocnonlynow = areaocnonly(regmap==ireg);
    beta_kost_regs(ireg) = nansum(nansum(kostbetanow.*areaocnonlynow))...
        ./nansum(nansum(areaocnonlynow));
end
beta_on_regs = nan(1,length(regnames));
for ireg=1:length(regnames)
    beta_onnow = beta_on(regmap==ireg);
    areaocnonlynow = areaocnonly(regmap==ireg);
    beta_on_regs(ireg) = nansum(nansum(beta_onnow.*areaocnonlynow))...
        ./nansum(nansum(areaocnonlynow));
end

plotnow=1;
if plotnow==1
f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 6 6]);
subplot(211);
bar(1:length(regnames),[beta_kost_regs' beta_on_regs']);
set(gca, 'XTickLabel', regnames);
ylabel('\beta');
legend({'Kost',['100-yr ' num2str(circfactor) ' circ']},'location','northoutside');
subplot(212);
bar(1:length(regnames),beta_on_regs-beta_kost_regs);
ylabel('\Delta\beta');
title(['100-yr ' num2str(circfactor) ' circ - baseline']);
set(gca, 'XTickLabel', regnames);
saveas(f,sprintf('betaendregbarchart_%.2fcirc.png',circfactor));
end

disp('reg mean beta, Kost:');
beta_kost_regs
disp(['reg mean beta, 100 yr fb on - Kost, ' num2str(circfactor) ' circ:']);
beta_on_regs-beta_kost_regs

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remin depth maps, 2 degree
% (model resolution)
%%%%%%%%%%%%%%%%%%%%%%%%%%
vardesnow = 'remin depths in [m] calc using PRiSM from 0.9 circ fb on, 100-yr beta map';
betamapnow = beta_on; 
lonsnow = lons2; 
latsnow = lats2; 
betanamenow = ['100yr' sprintf('%.2f',circfactor) 'circ'];
loadremindepthmap=1;
if loadremindepthmap==1
    load(['remindepthmap_' betanamenow 'beta.mat'],'remindepthmap');
else
    calc_remin_depths;
    % --> var you want is remindepthmap
end
rd_on = remindepthmap;

plotnow=1;
if plotnow==1
divnow = flipud(cbrewer('div','RdBu',21,'linear'));
seqnow = parula(20);
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
divplotseqnow = [0,1,0,1,0,1];
titlesnow = {'Baseline remin depth (=1/e of ez)',[num2str(circfactor) ' circ on-off'],...
    [num2str(circfactor) ' circ off'],[num2str(circfactor) ' circ off-baseline'],...
    [num2str(circfactor) ' circ on'],[num2str(circfactor) ' circ on-baseline']};
plotsnow = nan(length(lats2),length(lons2),6);
plotsnow(:,:,1) = rd_kost;
plotsnow(:,:,2) = rd_on-rd_kost;
plotsnow(:,:,3) = rd_kost;
plotsnow(:,:,4) = rd_kost-rd_kost;
plotsnow(:,:,5) = rd_on;
plotsnow(:,:,6) = rd_on-rd_kost;
f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 8.5 7]);
for isp=1:6
    ax=subplot(3,2,isp);
    m_proj(mapproj,'lon',[0 360],'lat',[min(lats2) max(lats2)]);
    m_pcolor(lons2,lats2,plotsnow(:,:,isp));
    m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
    c=colorbar; c.TickLength = 0.05; shading flat;
    m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
    if divplotseqnow(isp)==1
        caxis([-max(max(abs(plotsnow(:,:,isp)))) ...
                max(max(abs(plotsnow(:,:,isp))))]);
        caxis([-75 75]);
        colormap(ax,divnow);
    elseif divplotseqnow(isp)==0
        caxis([0 1650]);
        colormap(ax,seqnow);
    end
    hold on;
    nreg = max(regmap(:));
    for i = 1:nreg
        RR = double(regmap==i);
        RR(M3d(:,:,1)==0) = NaN;
        m_contour(lons2,lats2,RR,[.5 .5],'k','linewidth',1.5);
    end
    title(titlesnow{isp});
end
saveas(f,sprintf('remindepthendmaps_%.2fcirc_narrowcaxis.png',circfactor));
%saveas(f,sprintf('remindepthendmaps_%.2fcirc_widecaxis.png',circfactor));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate global mean remin depth
%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('global mean remin depth (=1/e of ez), Kost:')
nansum(nansum(rd_kost.*areaocnonly))./nansum(nansum(areaocnonly))
disp(['global mean remin depth (=1/e of ez) 100 yr fb on, ' num2str(circfactor) ' circ:']);
nansum(nansum(rd_on.*areaocnonly))./nansum(nansum(areaocnonly))

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate regional mean remin depth
%%%%%%%%%%%%%%%%%%%%%%%%%%
rd_kost_regs = nan(1,length(regnames));
for ireg=1:length(regnames)
    rd_kostnow = rd_kost(regmap==ireg);
    areaocnonlynow = areaocnonly(regmap==ireg);
    rd_kost_regs(ireg) = nansum(nansum(rd_kostnow.*areaocnonlynow))...
        ./nansum(nansum(areaocnonlynow));
end
rd_on_regs = nan(1,length(regnames));
for ireg=1:length(regnames)
    rd_onnow = rd_on(regmap==ireg);
    areaocnonlynow = areaocnonly(regmap==ireg);
    rd_on_regs(ireg) = nansum(nansum(rd_onnow.*areaocnonlynow))...
        ./nansum(nansum(areaocnonlynow));
end

plotnow=1;
if plotnow==1
f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 6 6]);
subplot(211);
bar(1:length(regnames),[rd_kost_regs' rd_on_regs']);
set(gca, 'XTickLabel', regnames);
ylabel('Remin depth [m]');
legend({'Kost',['100-yr ' num2str(circfactor) ' circ']},'location','northoutside');
subplot(212);
bar(1:length(regnames),rd_on_regs-rd_kost_regs);
ylabel('\Delta(remin depth) [m]');
title(['100-yr ' num2str(circfactor) ' circ - baseline']);
set(gca, 'XTickLabel', regnames);
saveas(f,sprintf('remindepthendregbarchart_%.2fcirc.png',circfactor));
end

disp('reg mean remin depth (depth = 1/e of ez), Kost:');
rd_kost_regs
disp(['reg mean remin depth (depth = 1/e of ez), 100 yr fb on - Kost, ' num2str(circfactor) ' circ:']);
rd_on_regs-rd_kost_regs
