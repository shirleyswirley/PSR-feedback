plotonly=0;

if plotonly==0

close all;
clear all;
setup_figs;

%-------------------------------
% Define desired circ rate change
% and analyzed flux depth
%-------------------------------
circfactor = 0.9;
fluxdepth = 90;

%-------------------------------
% Load grid variables + regions map
% + model output
%-------------------------------
% - Load grid variables + regions map
reg_struct = load([data_path ...
    'PRiSM_regions_2deg.mat']);
M3d = reg_struct.M3d; gridd = reg_struct.gridd;
lon2 = gridd.xt; lat2 = gridd.yt;
isurf = find(M3d(:,:,1)==1);
area_ocn_only = nan(size(gridd.Areat)); % m^2
area_ocn_only(isurf) = gridd.Areat(isurf);
area_oo_zsum = nansum(area_ocn_only,2);
tot_ocn_area = sum(gridd.Areat(isurf)); % [m^2]
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
% Reordered/abridged region numbers:
% 1 = AAZ, 2 = SAZ, 3 = STA, 4 = STP,
% 5 = IND, 6 = ETA, 7 = ETP, 8 = NA,
% 9 = NP
reg_names = [reg_struct.regnamesabbrev(1:4) ...
    reg_struct.regnamesabbrev(9) ...
    reg_struct.regnamesabbrev(5:8)];
nregs = length(reg_names);
%figure; fcolor(reg_map); shading flat; colorbar;

%----------Load beta and export maps
% - Baseline 1*circ, steady-state
bl_fname = [model_output_path 'BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_Kostbetainit_PSDfboff_epsnolims_1.00circ_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_spiralintoeqiternum0.mat'];
load(bl_fname,'output');
bl_ce = output.expCmapSS;
% Zonal mean
bl_cezm = squeeze(nanmean(bl_ce,2));
% Regional + global means
%needs mapnow,nregs,reg_map,area_ocn_only,tot_ocn_area
mapnow = bl_ce; calc_reg_global_means; bl_cerm = mapnowrm;

% - PSR fb off, new circ, 100 yr model output
off_fname = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfboffKostbeta_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor)];
load(off_fname,'output');
off_ce = output.expCmapnow(:,:,end);
% Zonal mean
off_cezm = squeeze(nanmean(off_ce,2));
% Regional + global means
%needs mapnow,nregs,reg_map,area_ocn_only,tot_ocn_area
mapnow = off_ce; calc_reg_global_means; off_cerm = mapnowrm;

% - PSD fb on in all regions,
% new circ, 100 yr model output 
on_fname = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocksm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor)];
load(on_fname,'output');
on_ce = output.expCmapnow(:,:,end);
on_beta_map = output.betamapnow(:,:,end);
% Zonal mean
on_cezm = squeeze(nanmean(on_ce,2));
% Regional + global means
%needs mapnow,nregs,reg_map,area_ocn_only,tot_ocn_area
mapnow = on_ce; calc_reg_global_means; on_cerm = mapnowrm;

% - Kost beta
% (applies to baseline 1*circ and
% any circ PSD fb off runs)
load([data_path 'Kost_beta_2deg.mat'],'beta_clim');
kost_beta_map = mean(beta_clim,3)';
kost_beta_map(M3d(:,:,1)==0)=NaN;

%----------Load/create normalized POC flux 3d maps
% - PSD fb on in all regions,
% new circ, 100 yr model output 
beta_name = ['kostfbon_100yr' sprintf('%.2f',circfactor) 'circ'];
load_flux_map=1;
if load_flux_map==1
    load([data_path 'norm_poc_flux_3dmap_' beta_name 'beta.mat'],...
        'norm_flux_map','z');
    on_nflux_map = norm_flux_map;
    on_nflux_map_z = z;
else
    var_des = ...
        ['Normalized POC flux profiles calc using PRiSM from 0.9 circ fb on, 100-yr beta map; model output used was' on_fname];
    [on_nflux_map, on_nflux_map_z] = calc_norm_flux_3dmap(...
        lon2,lat2,on_beta_map,beta_name,var_des);
end

% - Kost beta
% (applies to baseline 1*circ and
% any circ PSD fb off runs)
beta_name = 'kost';
load_flux_map=1;
if load_flux_map==1
    load([data_path 'norm_poc_flux_3dmap_' beta_name 'beta.mat'],...
        'norm_flux_map','z');
    kost_nflux_map = norm_flux_map;
    kost_nflux_map_z = z;
else
    var_des = ...
        'Normalized POC flux profiles calc using PRiSM from annual mean Kost beta map';
    [kost_nflux_map, kost_nflux_map_z] = calc_norm_flux_3dmap(...
        lon2,lat2,kost_beta_map,beta_name,var_des);
end

%-------------------------------
% Compute flux maps and means + fb strength
%-------------------------------
kost_flux_zidx = find(kost_nflux_map_z==fluxdepth);
on_flux_zidx = find(on_nflux_map_z==fluxdepth);

bl_flux = bl_ce.*kost_nflux_map(:,:,kost_flux_zidx);
off_flux = off_ce.*kost_nflux_map(:,:,kost_flux_zidx);
on_flux = on_ce.*on_nflux_map(:,:,on_flux_zidx); 

bl_fluxzm = squeeze(nanmean(bl_flux,2));
off_fluxzm = squeeze(nanmean(off_flux,2));
on_fluxzm = squeeze(nanmean(on_flux,2));

mapnow = bl_flux; calc_reg_global_means; bl_fluxrm = mapnowrm;
mapnow = off_flux; calc_reg_global_means; off_fluxrm = mapnowrm;
mapnow = on_flux; calc_reg_global_means; on_fluxrm = mapnowrm;

% - Calculate export actual reg/global, zonal mean fb strength
roundple = -3;
fbst_cerm = 100*roundn(on_cerm-off_cerm,roundple)./(off_cerm-bl_cerm);
fbst_cezm = 100*roundn(on_cezm-off_cezm,roundple)./(off_cezm-bl_cezm);

% - Calculate flux actual reg/global, zonal mean fb strength
roundplf = -3;
fbst_fluxrm = 100*roundn(on_fluxrm-off_fluxrm,roundplf)./(off_fluxrm-bl_fluxrm);
fbst_fluxzm = 100*roundn(on_fluxzm-off_fluxzm,roundplf)./(off_fluxzm-bl_fluxzm);

end

%-------------------------------
% Plot map figure
%-------------------------------
% - Define plot params
seqcmap = cbrewer('seq','YlGnBu',20,'linear');
divcmap = cbrewer('div','RdBu',21,'linear'); 
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
labelfontsize = 10;
titlefontsize = 12;
titlefontwt = 'bold';
cbticklen = 0.03;
cmin = -4.5; cmax = 1;
cminmax = 2.25;

f=figure;
set(f,'color','white','units','inches',...
    'position',[0.5 0.5 13 7],'resize','off','paperpositionmode','auto');

% - Baseline flux map
ax = subplot(331);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,log10(bl_flux),'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([cmin cmax]);
colormap(ax,seqcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('log10(Baseline flux [molC m^{-2} yr^{-1}])','fontsize',titlefontsize,'fontweight',titlefontwt);

% - fb off
ax = subplot(332);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,log10(off_flux),'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([cmin cmax]);
colormap(ax,seqcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('fb off','fontsize',titlefontsize,'fontweight',titlefontwt);

% - fb on
ax = subplot(333);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,log10(on_flux),'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([cmin cmax]);
colormap(ax,seqcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('fb on','fontsize',titlefontsize,'fontweight',titlefontwt);

% - fb on minus off
ax = subplot(334);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,on_flux-off_flux,'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([-0.2 0.2]);
colormap(ax,divcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('on-off [molC m^{-2} yr^{-1}]','fontsize',titlefontsize,'fontweight',titlefontwt);

% - fb off minus bl
ax = subplot(335);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,off_flux-bl_flux,'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([-cminmax cminmax]);
colormap(ax,divcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('off-bl','fontsize',titlefontsize,'fontweight',titlefontwt);

% - fb on minus bl
ax = subplot(336);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,on_flux-bl_flux,'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([-cminmax cminmax]);
colormap(ax,divcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('on-bl','fontsize',titlefontsize,'fontweight',titlefontwt);

% - (fb on minus off)/(fb off minus bl) = fb strength
ax = subplot(337);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,(on_flux-off_flux)./(off_flux-bl_flux),'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([-0.3 0.3]);
colormap(ax,divcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('(on-off)/(off-bl)=fb strength','fontsize',titlefontsize,'fontweight',titlefontwt);

% - (fb off minus bl)/bl
ax = subplot(338);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,(off_flux-bl_flux)./bl_flux,'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([-0.3 0.3]);
colormap(ax,divcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('(off-bl)/bl','fontsize',titlefontsize,'fontweight',titlefontwt);

% - (fb on minus bl)/bl
ax = subplot(339);
m_proj(mapproj,'lon',[0 360],'lat',[min(lat2) max(lat2)]);
m_contourf(lon2,lat2,(on_flux-bl_flux)./bl_flux,'linestyle','none');
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([-0.3 0.3]);
colormap(ax,divcmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
title('(on-bl)/bl','fontsize',titlefontsize,'fontweight',titlefontwt);

sgtitle(['Depth = ' num2str(fluxdepth) ' m']);

%-------------------------------
% Plot zonal, regional mean figure
%-------------------------------
tickfontsize = 9;
labelfontsize = 10; labelfontwt = 'normal';
titlefontsize = 11; titlefontwt = 'bold';
fbofflinestyle='-'; fbonlinestyle='--';
fbofflinewidth=2; fbonlinewidth=2;
cecol = [0 0 0.8]; fcol = [1 0.4 0];

f=figure; set(f,'Color','White',...
    'Units','Inches','Position',[1 1 8.5 10]);

%---------ROW 1
cemin = 0; cemax = 10;
ceticks = linspace(cemin,cemax,6);
fmin = 0; fmax = 10;
fticks = linspace(fmin,fmax,6);

% - Zonal avgs for baseline 1*circ case
subplot(421);
[ax,hce,hp]=plotyy(lat2,bl_cezm,lat2,bl_fluxzm);
set(hce,'Color',cecol,'LineStyle','-','LineWidth',2);
set(hp,'Color',fcol,'LineStyle','-','Linewidth',2);
set(ax(1),'YColor',cecol,'XLim',[-80 80],...
    'XTick',linspace(-80,80,9),'box','off',...
    'YLim',[cemin cemax],'YTick',ceticks,...
    'fontsize',tickfontsize);
set(ax(2),'YColor',fcol,'XLim',[-80 80],...
    'YLim',[fmin fmax],'YTick',fticks,...
    'fontsize',tickfontsize);
rl=refline(ax(1),0,ax(1).YLim(2)); rl.Color='k';
ylabel(ax(1),'Export [molC m^{-2} yr^{-1}]',...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
ylabel(ax(2),['Flux_{' num2str(fluxdepth) 'm} [molC m^{-2} yr^{-1}]'],...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
set(ax,'TickDir','out');
xlabel('Latitude','fontsize',labelfontsize);
title(['Baseline zonal mean '...
    sprintf('%s{%f %f %f}','\color[rgb]',cecol) 'export' ...
    '\color{black} & ' ...
    sprintf('%s{%f %f %f}','\color[rgb]',fcol) 'flux_{' num2str(fluxdepth) 'm}'],...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

% - Regional avgs for baseline 1*circ case
barwidth = 0.3; offset = barwidth/2;
subplot(422);
[ax,hce,hp]=plotyy(...
    (1:nregs+1)-offset,bl_cerm,...
    (1:nregs+1)+offset,bl_fluxrm,...
    @(x,y) bar(x,y,barwidth), @(x,y) bar(x,y,barwidth));
set(hce,'FaceColor',cecol);
set(hp,'FaceColor',fcol);
set(ax(1),'YColor',cecol,'YLim',[cemin cemax],'box','off',...
    'YTick',ceticks,'XLim',[0.5 nregs+1.5],...
    'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'fontsize',tickfontsize);
set(ax(2),'YColor',fcol,'YLim',[fmin fmax],...
    'YTick',fticks,'XLim',[0.5 nregs+1.5],...
    'XTick',[],...
    'fontsize',tickfontsize);
rl1=refline(ax(1),0,ax(1).YLim(2)); rl1.Color='k'; % adds box top edge
rl2=refline(ax(2),0,ax(2).YLim(1)); rl2.Color='k'; % make x-axis black
ylabel(ax(1),'Export [molC m^{-2} yr^{-1}]','FontSize',labelfontsize,...
    'FontWeight',labelfontwt);
ylabel(ax(2),['Flux_{' num2str(fluxdepth) 'm} [molC m^{-2} yr^{-1}]'],...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
set(ax,'TickDir','out');
xlabel('Region','fontsize',labelfontsize);
title(['Baseline regional mean '...
    sprintf('%s{%f %f %f}','\color[rgb]',cecol) 'export' ...
    '\color{black} & ' ...
    sprintf('%s{%f %f %f}','\color[rgb]',fcol) 'flux_{' num2str(fluxdepth) 'm}'],...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

%---------ROW 2
ylimsnow = sort([-18 1]);
%ylimsnow = sort(circadj*[-5.5 2]);
% - Relative change in zonal mean flux from baseline case
ax = subplot(423);
plot(lat2,100*roundn(on_fluxzm-bl_fluxzm,roundplf)...
    ./bl_fluxzm,...
    'color',fcol,'linestyle',fbonlinestyle,...
    'linewidth',fbonlinewidth); hold on;
plot(lat2,100*roundn(off_fluxzm-bl_fluxzm,roundplf)...
    ./bl_fluxzm,...
    'color',fcol,'linestyle',fbofflinestyle,...
    'linewidth',fbofflinewidth);
plot(lat2,zeros(size(lat2)),'k');
xlim([-80 80]); ylim(ylimsnow);
set(ax,'XTick',linspace(-80,80,9),...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Latitude','fontsize',labelfontsize);
ylabel('Relative change [%]','fontsize',labelfontsize);
title(['100-yr change in regional mean flux_{' num2str(fluxdepth) 'm}'],...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

% - Rel change in regional mean flux from baseline case
barwidth = 0.8;
ax = subplot(424);
solidbars = bar(...
    (1:nregs+1)-0.15,100*roundn(off_fluxrm-bl_fluxrm,roundplf)...
    ./bl_fluxrm,barwidth/3); hold on;
hatchedbars = bar(...
    (1:nregs+1)+0.15,100*roundn(on_fluxrm-bl_fluxrm,roundplf)...
    ./bl_fluxrm,barwidth/3,'facecolor','white');
colornow = fcol; format_fbonoff_bars; % needs colornow,hatchedbars,solidbars:
xlim([0.5 nregs+1.5]); ylim(ylimsnow);
set(ax,'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Region','fontsize',labelfontsize);
ylabel('Relative change [%]','fontsize',labelfontsize);
title(['100-yr change in regional mean flux_{' num2str(fluxdepth) 'm}'],...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

%---------ROW 3
ylimsnow = sort([-17 1]);
%ylimsnow = sort(circadj*[-20 1]);
% - Relative change in zonal mean ce from baseline case
ax = subplot(425);
plot(lat2,100*roundn(off_cezm-bl_cezm,roundplf)./bl_cezm,...
    'color',cecol,'linestyle',fbofflinestyle,...
    'linewidth',fbofflinewidth); hold on;
plot(lat2,100*roundn(on_cezm-bl_cezm,roundplf)./bl_cezm,...
    'color',cecol,'linestyle',fbonlinestyle,...
    'linewidth',fbonlinewidth);
plot(lat2,zeros(size(lat2)),'k');
xlim([-80 80]); ylim(ylimsnow);
set(ax,'XTick',linspace(-80,80,9),...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Latitude','fontsize',labelfontsize);
ylabel('Relative change [%]','fontsize',labelfontsize);
title('100-yr change in zonal mean export_ ',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

% - Rel change in regional mean ce from baseline case
barwidth = 0.8;
ax = subplot(426);
solidbars = bar(...
    (1:nregs+1)-0.15,100*roundn(off_cerm-bl_cerm,roundple)...
    ./bl_cerm,barwidth/3); hold on;
hatchedbars = bar(...
    (1:nregs+1)+0.15,100*roundn(on_cerm-bl_cerm,roundple)...
    ./bl_cerm,barwidth/3,'facecolor','white');
colornow = cecol; format_fbonoff_bars; % needs colornow,hatchedbars,solidbars:
xlim([0.5 nregs+1.5]); ylim(ylimsnow);
set(ax,'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Region','fontsize',labelfontsize);
ylabel('Relative change [%]','fontsize',labelfontsize);
title('100-yr change in regional mean export_ ',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

%---------ROW 4
ylimsnow = [-30 30];
% - Zonal mean fb strength, export and flux
ax = subplot(427);
plot(lat2,fbst_cezm,...
    'Color',cecol,'LineStyle','-','LineWidth',2); hold on;
plot(lat2,fbst_fluxzm,...
    'Color',fcol,'LineStyle','-','LineWidth',2);
xlim([-80 80]); ylim(ylimsnow);
set(ax,'XTick',linspace(-80,80,9),...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Latitude','fontsize',labelfontsize);
ylabel('Feedback strength [%]','fontsize',labelfontsize);
title('Zonal mean global feedback strength_ ',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

% - Regional mean fb strength, export and flux
barwidth=1; offset = 0.125;
ax = subplot(428);
b = bar(1:nregs+1,[fbst_cerm' fbst_fluxrm'],barwidth);
set(b(1),'facecolor',cecol);
set(b(2),'facecolor',fcol);
xlim([0.5 nregs+1.5]); ylim(ylimsnow);
set(ax,'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Region','fontsize',labelfontsize);
ylabel('Feedback strength [%]','fontsize',labelfontsize);
title('Regional mean global feedback strength_ ',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

print(f, [fig_save_path 'flux_at_' num2str(fluxdepth) 'm_depth.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'flux_at_' num2str(fluxdepth) 'm_depth.png'], '-dpng', '-r300');
