close all;
clear all;
setup_figs;

%-------------------------------
% Define desired circ rate change
%-------------------------------
circfactor = 0.9;
circadj = 1;
numyrs = 100;
figname = 'fig7';
%numyrs = 500;
%figname = 'fig7_500yr';

%circfactor = 1.1;
%circadj = -1;
%numyrs = 100;
%figname = 'fig7_1.1circ';
%numyrs = 500;
%figname = 'fig7_500yr_1.1circ';

%-------------------------------
% Load grid variables + regions map
%-------------------------------
reg_struct = load([data_path ...
    'PRiSM_regions_2deg.mat']);
M3d = reg_struct.M3d; gridd = reg_struct.gridd;
lon2 = gridd.xt; lat2 = gridd.yt;
depth = gridd.zt; iocn = find(M3d==1);
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
%figure; pcolor(reg_map); shading flat; colorbar;

% - Get depth layers above and below 200 m 
[~,didx]=min(abs(depth-200));
if depth(didx)>200
    didxdp = didx; didxsh = didx-1;
elseif depth(didx)<200
    didxsh = didx; didxdp = didx+1;
end

%-------------------------------
% Load and calc po4, export model output
%-------------------------------

% ----- Baseline 1*circ, steady-state
bl_fname = [model_output_path 'BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_Kostbetainit_PSDfboff_epsnolims_1.00circ_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_spiralintoeqiternum0.mat'];
load(bl_fname,'output');
bl_p = M3d*NaN; bl_p(iocn) = output.po4;
bl_ce = output.expCmapSS;
% - 200 m po4
newdepth = 200; map3dnow = bl_p;
%needs didxsh,didxdp,map3dnow,newdepth
interp_to_depth; bl_p200 = map2dnow;
% - Zonal means
bl_pzm = squeeze(nanmean(bl_p,2));
bl_p200zm = squeeze(nanmean(bl_p200,2));
bl_cezm = squeeze(nanmean(bl_ce,2));
% - Regional + global means
%needs mapnow,nregs,reg_map,area_ocn_only,tot_ocn_area
mapnow = bl_p200; calc_reg_global_means; bl_p200rm = mapnowrm;
mapnow = bl_ce; calc_reg_global_means; bl_cerm = mapnowrm;

% ----- PSR fb off, new circ, 100 yr model output
off_fname = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfboffKostbeta_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_1ov1000dt_ts%dyrs_prevNMinit.mat',circfactor,numyrs)];
load(off_fname,'output');
off_p = M3d*NaN; off_p(iocn) = output.po4end;
off_ce = output.expCmapnow(:,:,end);
% - 200 m po4
newdepth = 200; map3dnow = off_p;
%needs didxsh,didxdp,map3dnow,newdepth
interp_to_depth; off_p200 = map2dnow;
% - Zonal means
off_pzm = squeeze(nanmean(off_p,2));
off_p200zm = squeeze(nanmean(off_p200,2));
off_cezm = squeeze(nanmean(off_ce,2));
% - Regional + global means
%needs mapnow,nregs,reg_map,area_ocn_only,tot_ocn_area
mapnow = off_p200; calc_reg_global_means; off_p200rm = mapnowrm;
mapnow = off_ce; calc_reg_global_means; off_cerm = mapnowrm;

% ----- Global PSR fb on, new circ, 100 yr model output 
on_fname = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocksm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts%dyrs_prevNMinit.mat',circfactor,numyrs)];
load(on_fname,'output');
on_p = M3d*NaN; on_p(iocn) = output.po4end;
on_ce = output.expCmapnow(:,:,end);
% - 200 m po4
newdepth = 200; map3dnow = on_p;
%needs didxsh,didxdp,map3dnow,newdepth
interp_to_depth; on_p200 = map2dnow;
% - Zonal means
on_pzm = squeeze(nanmean(on_p,2));
on_p200zm = squeeze(nanmean(on_p200,2));
on_cezm = squeeze(nanmean(on_ce,2));
% - Regional + global means 
%needs mapnow,nregs,reg_map,area_ocn_only,tot_ocn_area
mapnow = on_p200; calc_reg_global_means; on_p200rm = mapnowrm;
mapnow = on_ce; calc_reg_global_means; on_cerm = mapnowrm;

%-------------------------------
% Calculate actual export-based fb strength
% + predicted po4-based fb strength
%-------------------------------
deltacircovcirc = (circfactor-1)/1;

% - Calculate actual reg/global, zonal mean fb strength
roundple = -3;
fbst_cerm = 100*roundn(on_cerm-off_cerm,roundple)./(off_cerm-bl_cerm);
fbst_cezm = 100*roundn(on_cezm-off_cezm,roundple)./(off_cezm-bl_cezm);

% - Calculate predicted reg/global, zonal mean fb strength
roundplp = -3;
fbst_p200rm = 100*...
    roundn(on_p200rm-off_p200rm,roundplp)./bl_p200rm...
    ./(deltacircovcirc + roundn(off_p200rm-bl_p200rm,roundplp)./bl_p200rm);
fbst_p200zm = 100*...
    roundn(on_p200zm-off_p200zm,roundplp)./bl_p200zm...
    ./(deltacircovcirc + roundn(off_p200zm-bl_p200zm,roundplp)./bl_p200zm);

%-------------------------------
% Plot figure
%-------------------------------
tickfontsize = 9;
labelfontsize = 10; labelfontwt = 'normal';
titlefontsize = 11; titlefontwt = 'bold';
fbofflinestyle='-'; fbonlinestyle='--';
fbofflinewidth=2; fbonlinewidth=2;
cecol = [0 0 0.8]; pcol = [1 0.4 0];

f=figure; set(f,'Color','White',...
    'Units','Inches','Position',[1 1 8.5 10]);

%---------ROW 1
cemin = 0; cemax = 10;
ceticks = linspace(cemin,cemax,6);
pmin = 0; pmax = 3;
pticks = linspace(pmin,pmax,6);

% - Zonal avgs for baseline 1*circ case
subplot(421);
[ax,hce,hp]=plotyy(lat2,bl_cezm,lat2,bl_p200zm);
set(hce,'Color',cecol,'LineStyle','-','LineWidth',2);
set(hp,'Color',pcol,'LineStyle','-','Linewidth',2);
set(ax(1),'YColor',cecol,'XLim',[-80 80],...
    'XTick',linspace(-80,80,9),'box','off',...
    'YLim',[cemin cemax],'YTick',ceticks,...
    'fontsize',tickfontsize);
set(ax(2),'YColor',pcol,'XLim',[-80 80],...
    'YLim',[pmin pmax],'YTick',pticks,...
    'fontsize',tickfontsize);
rl=refline(ax(1),0,ax(1).YLim(2)); rl.Color='k';
ylabel(ax(1),'Export [molC m^{-2} yr^{-1}]',...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
ylabel(ax(2),'P_{200m} [mmolP m^{-3}]',...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
set(ax,'TickDir','out');
xlabel('Latitude','fontsize',labelfontsize);
title(['Baseline zonal mean '...
    sprintf('%s{%f %f %f}','\color[rgb]',cecol) 'export' ...
    '\color{black} & ' ...
    sprintf('%s{%f %f %f}','\color[rgb]',pcol) 'P_{200m}'],...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

% - Regional avgs for baseline 1*circ case
barwidth = 0.3; offset = barwidth/2;
subplot(422);
[ax,hce,hp]=plotyy(...
    (1:nregs+1)-offset,bl_cerm,...
    (1:nregs+1)+offset,bl_p200rm,...
    @(x,y) bar(x,y,barwidth), @(x,y) bar(x,y,barwidth));
set(hce,'FaceColor',cecol);
set(hp,'FaceColor',pcol);
set(ax(1),'YColor',cecol,'YLim',[cemin cemax],'box','off',...
    'YTick',ceticks,'XLim',[0.5 nregs+1.5],...
    'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'fontsize',tickfontsize);
set(ax(2),'YColor',pcol,'YLim',[pmin pmax],...
    'YTick',pticks,'XLim',[0.5 nregs+1.5],...
    'XTick',[],...
    'fontsize',tickfontsize);
rl1=refline(ax(1),0,ax(1).YLim(2)); rl1.Color='k'; % adds box top edge
rl2=refline(ax(2),0,ax(2).YLim(1)); rl2.Color='k'; % make x-axis black
ylabel(ax(1),'Export [molC m^{-2} yr^{-1}]','FontSize',labelfontsize,...
    'FontWeight',labelfontwt);
ylabel(ax(2),'P_{200m} [mmolP m^{-3}]',...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
set(ax,'TickDir','out');
xlabel('Region','fontsize',labelfontsize);
title(['Baseline regional mean '...
    sprintf('%s{%f %f %f}','\color[rgb]',cecol) 'export' ...
    '\color{black} & ' ...
    sprintf('%s{%f %f %f}','\color[rgb]',pcol) 'P_{200m}'],...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

%---------ROW 2
ylimsnow = sort(circadj*[-4.2 1]);
%ylimsnow = sort(circadj*[-5.5 2]);
% - Relative change in zonal mean p200 from baseline case
ax = subplot(423);
plot(lat2,100*roundn(on_p200zm-bl_p200zm,roundplp)...
    ./bl_p200zm,...
    'color',pcol,'linestyle',fbonlinestyle,...
    'linewidth',fbonlinewidth); hold on;
plot(lat2,100*roundn(off_p200zm-bl_p200zm,roundplp)...
    ./bl_p200zm,...
    'color',pcol,'linestyle',fbofflinestyle,...
    'linewidth',fbofflinewidth);
plot(lat2,zeros(size(lat2)),'k');
xlim([-80 80]); ylim(ylimsnow);
set(ax,'XTick',linspace(-80,80,9),...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Latitude','fontsize',labelfontsize);
ylabel('Relative change [%]','fontsize',labelfontsize);
title('100-yr change in zonal mean P_{200m}',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

% - Rel change in regional mean p200 from baseline case
barwidth = 0.8;
ax = subplot(424);
solidbars = bar(...
    1:nregs+1,100*roundn(off_p200rm-bl_p200rm,roundplp)...
    ./bl_p200rm,barwidth); hold on;
hatchedbars = bar(...
    1:nregs+1,100*roundn(on_p200rm-bl_p200rm,roundplp)...
    ./bl_p200rm,barwidth,'facecolor','white');
colornow = pcol; format_fbonoff_bars; % needs colornow,hatchedbars,solidbars:
xlim([0.5 nregs+1.5]); ylim(ylimsnow);
set(ax,'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Region','fontsize',labelfontsize);
ylabel('Relative change [%]','fontsize',labelfontsize);
title('100-yr change in regional mean P_{200m}',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

%---------ROW 3
ylimsnow = sort(circadj*[-17 1]);
%ylimsnow = sort(circadj*[-20 1]);
% - Relative change in zonal mean ce from baseline case
ax = subplot(425);
plot(lat2,100*roundn(off_cezm-bl_cezm,roundplp)./bl_cezm,...
    'color',cecol,'linestyle',fbofflinestyle,...
    'linewidth',fbofflinewidth); hold on;
plot(lat2,100*roundn(on_cezm-bl_cezm,roundplp)./bl_cezm,...
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
    1:nregs+1,100*roundn(off_cerm-bl_cerm,roundple)...
    ./bl_cerm,barwidth); hold on;
hatchedbars = bar(...
    1:nregs+1,100*roundn(on_cerm-bl_cerm,roundple)...
    ./bl_cerm,barwidth,'facecolor','white');
colornow = cecol; format_fbonoff_bars; % needs colornow,hatchedbars,solidbars:
xlim([0.5 nregs+1.5]); ylim(ylimsnow);
set(ax,'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Region','fontsize',labelfontsize);
ylabel('Relative change [%]','fontsize',labelfontsize);
title('100-yr change in regional mean export_ ',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

%---------ROW 4
ylimsnow = [0 24];
%ylimsnow = [0 30];
% - Zonal mean fb strength, actual and predicted
ax = subplot(427);
plot(lat2,abs(fbst_cezm),...
    'Color',cecol,'LineStyle','-','LineWidth',2); hold on;
plot(lat2,abs(fbst_p200zm),...
    'Color',pcol,'LineStyle','-','LineWidth',2);
xlim([-80 80]); ylim(ylimsnow);
set(ax,'XTick',linspace(-80,80,9),...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Latitude','fontsize',labelfontsize);
ylabel('Feedback strength [%]','fontsize',labelfontsize);
title('Zonal mean global feedback strength_ ',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

% - Regional mean fb strength, actual and predicted
barwidth=1; offset = 0.125;
ax = subplot(428);
b = bar(1:nregs+1,abs([fbst_cerm' fbst_p200rm']),barwidth);
set(b(1),'facecolor',cecol);
set(b(2),'facecolor',pcol);
xlim([0.5 nregs+1.5]); ylim(ylimsnow);
set(ax,'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Region','fontsize',labelfontsize);
ylabel('Feedback strength [%]','fontsize',labelfontsize);
title('Regional mean global feedback strength_ ',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

fbst_cerm
print(f, [fig_save_path figname '_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path figname '_PSRfbpaper_final.png'], '-dpng', '-r300');
