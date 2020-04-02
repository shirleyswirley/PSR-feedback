addpath(genpath('/graid1/shirlleu/matlabroutines'))

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get relevant matfiles
%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------
% Tom's regions map (same grid as all other matfiles below)
%---------------------
regmapst = load('/graid1/shirlleu/PRiSM_GCM/data/PRiSMregions.mat');
regmapgrid = regmapst.gridd;
regmap = regmapst.R2d;
%figure; pcolor(regmap); shading flat; colorbar;

%---------------------
% Baseline 1*circ
%---------------------
fnamebl = '/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/newtonmethod/spiralintoeq/prog/PSDfboff/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_Kostbetainit_PSDfboff_epsnolims_1.00circ_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_spiralintoeqiternum0.mat';

%---------------------
% PSD fb off
%---------------------
fnameoff = '/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/PSDfboff/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfboffKostbeta_epsbtwn2and6p5_0.90circ_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_1ov1000dt_ts100yrs_prevNMinit.mat';

%---------------------
% PSD fb on, all regions turned on
%---------------------
fnameon = '/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocksm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_0.90circ_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat';

%---------------------
% PSD fb on, all regions turned on,
% beta vs. export - 1stdev
%---------------------
fnameonms = '/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdbelowm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_0.90circ_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat';

%---------------------
% PSD fb on, all regions turned on,
% beta vs. export + 1stdev
%---------------------
fnameonps = '/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdabovem1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_0.90circ_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load relevant variables
%%%%%%%%%%%%%%%%%%%%%%%%%%

% - Common to all files
load(fnamebl,'output');
gridd = output.grid;
lats = gridd.yt;
lons = gridd.xt;
M3d = output.M3d;
iocn = find(M3d==1);
isurf = find(M3d(:,:,1)==1);

area = gridd.Areat; % m^2
areaocnonly = nan(size(area));
areaocnonly(isurf) = area(isurf);
areazonalsum = nansum(areaocnonly,2);

depths = gridd.zt;
[~,tcdepthidx]=min(abs(depths-200));
[~,deepdepthidx]=min(abs(depths-2000));
% --> gets 187.5m depth

% - Baseline 
p_bl = M3d*NaN;
p_bl(iocn) = output.po4;
ps_bl = squeeze(p_bl(:,:,1)); % surface
ptc_bl = squeeze(p_bl(:,:,tcdepthidx)); % thermocline
pzm_bl = squeeze(nanmean(p_bl,2)); % zonal mean
ec_bl = output.expCmapSS;

% - Feedback off
load(fnameoff,'output');
p_off = M3d*NaN;
p_off(iocn) = output.po4end;
ps_off = squeeze(p_off(:,:,1)); % surface
ptc_off = squeeze(p_off(:,:,tcdepthidx)); % thermocline
pzm_off = squeeze(nanmean(p_off,2)); % zonal mean
ec_off = output.expCmapnow(:,:,end);

% - Feedback on
load(fnameon,'output');
p_on = M3d*NaN;
p_on(iocn) = output.po4end;
ps_on = squeeze(p_on(:,:,1)); % surface
ptc_on = squeeze(p_on(:,:,tcdepthidx)); % thermocline
pzm_on = squeeze(nanmean(p_on,2)); % zonal mean
ec_on = output.expCmapnow(:,:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% PO4 zonal mean transects
%%%%%%%%%%%%%%%%%%%%%%%%%%
divnow = flipud(cbrewer('div','RdBu',20,'linear'));
seqnow = parula(20);
divplotseqnow = [0,1,0,1,0,1];
titlesnow = {'Baseline zonal mean PO4','On-off',...
    'Off','Off-baseline','On','On-baseline'};
didxsnow = 1:deepdepthidx;
plotsnow = nan(length(depths(didxsnow)),length(lats),6); 
plotsnow(:,:,1) = pzm_bl(:,didxsnow)';
plotsnow(:,:,2) = (pzm_on(:,didxsnow)-pzm_off(:,didxsnow))';
plotsnow(:,:,3) = pzm_off(:,didxsnow)';
plotsnow(:,:,4) = (pzm_off(:,didxsnow)-pzm_bl(:,didxsnow))';
plotsnow(:,:,5) = pzm_on(:,didxsnow)';
plotsnow(:,:,6) = (pzm_on(:,didxsnow)-pzm_bl(:,didxsnow))';
f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 8.5 7]);
for isp=1:6
    ax=subplot(3,2,isp);
    pcolor(lats,depths(didxsnow),plotsnow(:,:,isp));
    c=colorbar; c.TickLength = 0.05; shading flat;
    xticks(-80:20:80); yticks(250:250:2000); grid on;
    set(gca,'YDir','reverse','TickLength',[0.03 0.035],'layer','top');
    xlabel('Lat'); ylabel('Depth [m]');
    if divplotseqnow(isp)==1
        caxis([-max(max(abs(plotsnow(:,:,isp)))) ...
                max(max(abs(plotsnow(:,:,isp))))]);
        colormap(ax,divnow);
    elseif divplotseqnow(isp)==0
        colormap(ax,seqnow);
    end
    title(titlesnow{isp});
end
saveas(f,'po4zonalmeans.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%
% PO4 surface maps
%%%%%%%%%%%%%%%%%%%%%%%%%%
divnow = flipud(cbrewer('div','RdBu',20,'linear'));
seqnow = parula(20);
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
divplotseqnow = [0,1,0,1,0,1];
titlesnow = {'Baseline surface PO4','On-off',...
    'Off','Off-baseline','On','On-baseline'};
plotsnow = nan(length(lats),length(lons),6);
plotsnow(:,:,1) = ps_bl;
plotsnow(:,:,2) = (ps_on-ps_off);
plotsnow(:,:,3) = ps_off;
plotsnow(:,:,4) = (ps_off-ps_bl);
plotsnow(:,:,5) = ps_on;
plotsnow(:,:,6) = (ps_on-ps_bl);
f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 8.5 7]);
for isp=1:6
    ax=subplot(3,2,isp);
    m_proj(mapproj,'lon',[0 360],'lat',[min(lats) max(lats)]);
    m_pcolor(lons,lats,plotsnow(:,:,isp));
    m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
    c=colorbar; c.TickLength = 0.05; shading flat;
    m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
    if divplotseqnow(isp)==1
        caxis([-max(max(abs(plotsnow(:,:,isp)))) ...
                max(max(abs(plotsnow(:,:,isp))))]);
        colormap(ax,divnow);
    elseif divplotseqnow(isp)==0
        colormap(ax,seqnow);
    end
    hold on;
    nreg = max(regmap(:));
    for i = 1:nreg
        RR = double(regmap==i);
        RR(M3d(:,:,1)==0) = NaN;
        m_contour(lons,lats,RR,[.5 .5],'k','linewidth',1.5);
    end
    title(titlesnow{isp});
end
saveas(f,'po4surfmaps.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%
% PO4 thermocline maps
%%%%%%%%%%%%%%%%%%%%%%%%%%
divnow = flipud(cbrewer('div','RdBu',20,'linear'));
seqnow = parula(20);
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
divplotseqnow = [0,1,0,1,0,1];
titlesnow = {'Baseline 187.5m PO4','On-off',...
    'Off','Off-baseline','On','On-baseline'};
plotsnow = nan(length(lats),length(lons),6);
plotsnow(:,:,1) = ptc_bl;
plotsnow(:,:,2) = (ptc_on-ptc_off);
plotsnow(:,:,3) = ptc_off;
plotsnow(:,:,4) = (ptc_off-ptc_bl);
plotsnow(:,:,5) = ptc_on;
plotsnow(:,:,6) = (ptc_on-ptc_bl);
f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 8.5 7]);
for isp=1:6
    ax=subplot(3,2,isp);
    m_proj(mapproj,'lon',[0 360],'lat',[min(lats) max(lats)]);
    m_pcolor(lons,lats,plotsnow(:,:,isp));
    m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
    c=colorbar; c.TickLength = 0.05; shading flat;
    m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
    if divplotseqnow(isp)==1
        caxis([-max(max(abs(plotsnow(:,:,isp)))) ...
                max(max(abs(plotsnow(:,:,isp))))]);
        colormap(ax,divnow);
    elseif divplotseqnow(isp)==0
        colormap(ax,seqnow);
    end
    hold on;
    nreg = max(regmap(:));
    for i = 1:nreg
        RR = double(regmap==i);
        RR(M3d(:,:,1)==0) = NaN;
        m_contour(lons,lats,RR,[.5 .5],'k','linewidth',1.5);
    end
    title(titlesnow{isp});
end
saveas(f,'po4tcmaps.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%
% deltaE terms, fb on
%%%%%%%%%%%%%%%%%%%%%%%%%%
divnow = flipud(cbrewer('div','RdBu',10,'linear'));
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
divplotseqnow = [1,1,1,1,1,1];
titlesnow = {'\Delta(w/w_{bl}) P_{tc,bl}','\DeltaP_{tc,on}',...
    '\Delta(w/w_{bl}) \DeltaP_{tc,on}',...
    '(1+2+3)-\DeltaE_{on}','1+2+3',...
    '\DeltaE_{on}'};
plotsnow = nan(length(lats),length(lons),6);
plotsnow(:,:,1) = -0.1*ptc_bl;
plotsnow(:,:,2) = ptc_on-ptc_bl; 
plotsnow(:,:,3) = -0.1*(ptc_on-ptc_bl);
plotsnow(:,:,4) = -0.1*ptc_bl + ptc_on-ptc_bl + -0.1*(ptc_on-ptc_bl) - (ec_on-ec_bl);
plotsnow(:,:,5) = -0.1*ptc_bl + ptc_on-ptc_bl + -0.1*(ptc_on-ptc_bl);
plotsnow(:,:,6) = ec_on-ec_bl;
f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 8.5 7]);
for isp=1:6
    ax=subplot(3,2,isp);
    m_proj(mapproj,'lon',[0 360],'lat',[min(lats) max(lats)]);
    m_pcolor(lons,lats,plotsnow(:,:,isp));
    m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
    c=colorbar; c.TickLength = 0.05; shading flat;
    m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
    if divplotseqnow(isp)==1
        caxis([-max(max(abs(plotsnow(:,:,isp)))) ...
                max(max(abs(plotsnow(:,:,isp))))]);
        colormap(ax,divnow);
    elseif divplotseqnow(isp)==0
        colormap(ax,seqnow);
    end
    hold on;
    nreg = max(regmap(:));
    for i = 1:nreg
        RR = double(regmap==i);
        RR(M3d(:,:,1)==0) = NaN;
        m_contour(lons,lats,RR,[.5 .5],'k','linewidth',1.5);
    end
    title(titlesnow{isp});
end
saveas(f,'deltaEtermfbonmaps.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%
% deltaE terms, fb off
%%%%%%%%%%%%%%%%%%%%%%%%%%
divnow = flipud(cbrewer('div','RdBu',10,'linear'));
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
divplotseqnow = [1,1,1,1,1,1];
titlesnow = {'\Delta(w/w_{bl}) P_{tc,bl}','\DeltaP_{tc,off}',...
    '\Delta(w/w_{bl}) \DeltaP_{tc,off}',...
    '(1+2+3)-\DeltaE_{off}','1+2+3',...
    '\DeltaE_{off}'};
plotsnow = nan(length(lats),length(lons),6);
plotsnow(:,:,1) = -0.1*ptc_bl;
plotsnow(:,:,2) = ptc_off-ptc_bl; 
plotsnow(:,:,3) = -0.1*(ptc_off-ptc_bl);
plotsnow(:,:,4) = -0.1*ptc_bl + ptc_off-ptc_bl + -0.1*(ptc_off-ptc_bl) - (ec_off-ec_bl);
plotsnow(:,:,5) = -0.1*ptc_bl + ptc_off-ptc_bl + -0.1*(ptc_off-ptc_bl);
plotsnow(:,:,6) = ec_off-ec_bl;
f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 8.5 7]);
for isp=1:6
    ax=subplot(3,2,isp);
    m_proj(mapproj,'lon',[0 360],'lat',[min(lats) max(lats)]);
    m_pcolor(lons,lats,plotsnow(:,:,isp));
    m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
    c=colorbar; c.TickLength = 0.05; shading flat;
    m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
    if divplotseqnow(isp)==1
        caxis([-max(max(abs(plotsnow(:,:,isp)))) ...
                max(max(abs(plotsnow(:,:,isp))))]);
        colormap(ax,divnow);
    elseif divplotseqnow(isp)==0
        colormap(ax,seqnow);
    end
    hold on;
    nreg = max(regmap(:));
    for i = 1:nreg
        RR = double(regmap==i);
        RR(M3d(:,:,1)==0) = NaN;
        m_contour(lons,lats,RR,[.5 .5],'k','linewidth',1.5);
    end
    title(titlesnow{isp});
end
saveas(f,'deltaEtermfboffmaps.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%
% (deltaEon-deltaEoff)/deltaEoff terms
%%%%%%%%%%%%%%%%%%%%%%%%%%
divnow = flipud(cbrewer('div','RdBu',20,'linear'));
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
divplotseqnow = [1,1,1,1,1,1,1,1,1];
titlesnow = {'(\DeltaE_{on}-\DeltaE_{off}) / \DeltaE_{off}',...
    '\DeltaE_{on}-\DeltaE_{off}', '\DeltaE_{off}',...
    '5/6',...
    '(1 + (\Deltaw/w_{bl}))(\DeltaP_{tc,on}-\DeltaP_{tc,off})', ...
    '\Deltaw/w_{bl} P_{tc,bl} + \DeltaP_{tc,off} + \Deltaw/w_{bl}\DeltaP_{tc,off}',...
    '1-4',...
    '\DeltaP_{tc,on}-\DeltaP_{tc,off}',...
    '(\Deltaw/w_{bl})(\DeltaP_{tc,on}-\DeltaP_{tc,off})'};
plotsnow = nan(length(lats),length(lons),6);
plotsnow(:,:,1) = (ec_on-ec_off)./(ec_off-ec_bl);
plotsnow(:,:,2) = ec_on-ec_off;
plotsnow(:,:,3) = ec_off-ec_bl;
plotsnow(:,:,4) = 0.9*(ptc_on-ptc_off)./...
    (-0.1*ptc_bl + ptc_off-ptc_bl + -0.1*(ptc_off-ptc_bl));
plotsnow(:,:,5) = 0.9*(ptc_on-ptc_off);
plotsnow(:,:,6) = -0.1*ptc_bl + ptc_off-ptc_bl + -0.1*(ptc_off-ptc_bl);
plotsnow(:,:,7) = (ec_on-ec_off)./(ec_off-ec_bl) - ...
    0.9*(ptc_on-ptc_off)./...
    (-0.1*ptc_bl + ptc_off-ptc_bl + -0.1*(ptc_off-ptc_bl));
plotsnow(:,:,8) = ptc_on-ptc_off;
plotsnow(:,:,9) = -0.1*(ptc_on-ptc_off);
f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 12.75 7]);
for isp=1:9
    ax=subplot(3,3,isp);
    m_proj(mapproj,'lon',[0 360],'lat',[min(lats) max(lats)]);
    m_pcolor(lons,lats,plotsnow(:,:,isp));
    m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
    c=colorbar; c.TickLength = 0.05; shading flat;
    m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
    if divplotseqnow(isp)==1
        if isp==1 | isp==4 | isp==7 
            cminmaxnow = 0.3;
        else 
            cminmaxnow = max(max(abs(plotsnow(:,:,isp))));
        end
        caxis([-cminmaxnow cminmaxnow]);
        colormap(ax,divnow);
    elseif divplotseqnow(isp)==0
        colormap(ax,seqnow);
    end
    hold on;
    nreg = max(regmap(:));
    for i = 1:nreg
        RR = double(regmap==i);
        RR(M3d(:,:,1)==0) = NaN;
        m_contour(lons,lats,RR,[.5 .5],'k','linewidth',1.5);
    end
    title(titlesnow{isp});
end
saveas(f,'fbtermsmaps.png')

