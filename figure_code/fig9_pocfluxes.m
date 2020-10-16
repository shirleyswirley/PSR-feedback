plotonly=0;

if plotonly==0

close all;
clear all;
setup_figs;

%-------------------------------
% Define desired circ rate change
% and analyzed flux depth
%-------------------------------
figname = 'fig9';
circfactor = 0.9;
fluxdepth = 725; % m below bottom of euphotic zone (75 m)

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
        ['Normalized POC flux profiles calc using PRiSM from ' sprintf('%.2f',circfactor) ' circ fb on, 100-yr beta map; model output used was' on_fname];
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

bl_flux = nan(size(kost_nflux_map));
off_flux = nan(size(kost_nflux_map));
for iz=1:length(kost_nflux_map_z)
    bl_flux(:,:,iz) = bl_ce.*kost_nflux_map(:,:,iz);
    off_flux(:,:,iz) = off_ce.*kost_nflux_map(:,:,iz);
end
on_flux = nan(size(on_nflux_map));
for iz=1:length(on_nflux_map_z)
    on_flux(:,:,iz) = on_ce.*on_nflux_map(:,:,iz);
end

bl_fluxzm = nan(length(kost_nflux_map_z),length(lat2));
off_fluxzm = nan(length(kost_nflux_map_z),length(lat2));
for iz=1:length(kost_nflux_map_z)
    bl_fluxzm(iz,:) = nanmean(bl_flux(:,:,iz),2);
    off_fluxzm(iz,:) = nanmean(off_flux(:,:,iz),2);
end
on_fluxzm = nan(length(on_nflux_map_z),length(lat2));
for iz=1:length(on_nflux_map_z)
    on_fluxzm(iz,:) = nanmean(on_flux(:,:,iz),2);
end

bl_fluxrm = nan(length(kost_nflux_map_z),nregs+1);
off_fluxrm = nan(length(kost_nflux_map_z),nregs+1);
for iz=1:length(kost_nflux_map_z)
    mapnow = bl_flux(:,:,iz); calc_reg_global_means; bl_fluxrm(iz,:) = mapnowrm;
    mapnow = off_flux(:,:,iz); calc_reg_global_means; off_fluxrm(iz,:) = mapnowrm;
end
on_fluxrm = nan(length(on_nflux_map_z),nregs+1);
for iz=1:length(on_nflux_map_z)
    mapnow = on_flux(:,:,iz); calc_reg_global_means; on_fluxrm(iz,:) = mapnowrm;
end

% - Calculate flux actual reg/global, zonal mean fb strength
% --> assume kost_nflux_map_z == on_nflux_map_z
fbst_fluxzm = nan(length(kost_nflux_map_z),length(lat2)); 
fbst_fluxrm = nan(length(kost_nflux_map_z),nregs+1); 
roundplf = -3;
for iz=1:length(kost_nflux_map_z)
    fbst_fluxzm(iz,:) = 100*roundn(on_fluxzm(iz,:)-off_fluxzm(iz,:),roundplf)...
        ./(off_fluxzm(iz,:)-bl_fluxzm(iz,:));
    fbst_fluxrm(iz,:) = 100*roundn(on_fluxrm(iz,:)-off_fluxrm(iz,:),roundplf)...
        ./(off_fluxrm(iz,:)-bl_fluxrm(iz,:));
end

% - Calculate reg/global, zonal mean depths at which fb strength
% switches from negative to positive
% (fbstsd = feedback strength switch depth)
fbstsd_fluxzm = nan(nregs+1,1);
for ilat = 1:length(lat2)
    fbstnow = fbst_fluxzm(:,ilat);
    idx_below = find(fbstnow>=0, 1);
    if isempty(idx_below) | idx_below==1
        fbstsd_fluxzm(ilat) = NaN;
    elseif fbstnow(idx_below)==0
        fbstsd_fluxzm(ilat) = kost_nflux_map_z(idx_below);
    else
        fbstsd_fluxzm(ilat) = interp1([fbstnow(idx_below-1) fbstnow(idx_below)],...
                        [kost_nflux_map_z(idx_below-1) kost_nflux_map_z(idx_below)],0);
    end
end
fbstsd_fluxrm = nan(nregs+1,1); 
for ireg = 1:nregs+1
    fbstnow = fbst_fluxrm(:,ireg);
    idx_below = find(fbstnow>=0, 1);
    if isempty(idx_below) | idx_below==1
        fbstsd_fluxrm(ireg) = NaN;
    elseif fbstnow(idx_below)==0
        fbstsd_fluxrm(ireg) = kost_nflux_map_z(idx_below);
    else
        fbstsd_fluxrm(ireg) = interp1([fbstnow(idx_below-1) fbstnow(idx_below)],...
                        [kost_nflux_map_z(idx_below-1) kost_nflux_map_z(idx_below)],0);
    end
end

end % end if plotonly==0

%-------------------------------
% Plot zonal, regional mean figure
%-------------------------------
tickfontsize = 9;
labelfontsize = 10; labelfontwt = 'normal';
titlefontsize = 11; titlefontwt = 'bold';
fbofflinestyle='-'; fbonlinestyle='--';
fbofflinewidth=2; fbonlinewidth=2;
fcol = [0 0 0.8];

f=figure; set(f,'Color','White',...
    'Units','Inches','Position',[1 1 8.5 5]);

%---------ROW 1
ylimsnow = sort([-18 1]);
%ylimsnow = sort(circadj*[-5.5 2]);

% - Relative change in zonal mean flux from baseline case
ax = subplot(221);
plot(lat2,100*...
        roundn(on_fluxzm(on_flux_zidx,:)-bl_fluxzm(kost_flux_zidx,:),roundplf)...
        ./bl_fluxzm(kost_flux_zidx,:),...
    'color',fcol,'linestyle',fbonlinestyle,...
    'linewidth',fbonlinewidth); hold on;
plot(lat2,100*...
        roundn(off_fluxzm(kost_flux_zidx,:)-bl_fluxzm(kost_flux_zidx,:),roundplf)...
        ./bl_fluxzm(kost_flux_zidx,:),...
    'color',fcol,'linestyle',fbofflinestyle,...
    'linewidth',fbofflinewidth);
plot(lat2,zeros(size(lat2)),'k');
xlim([-80 80]); ylim(ylimsnow);
set(ax,'XTick',linspace(-80,80,9),...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Latitude','fontsize',labelfontsize);
ylabel('Relative change [%]','fontsize',labelfontsize);
title(['100-yr change in zonal mean flux at ' num2str(fluxdepth+75) ' m'],...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

% - Rel change in regional mean flux from baseline case
barwidth = 0.8;
ax = subplot(222);
solidbars = bar((1:nregs+1)-0.15,...
    100*roundn(off_fluxrm(kost_flux_zidx,:)-bl_fluxrm(kost_flux_zidx,:),roundplf)...
        ./bl_fluxrm(kost_flux_zidx,:),barwidth/3); hold on;
hatchedbars = bar((1:nregs+1)+0.15,...
    100*roundn(on_fluxrm(on_flux_zidx,:)-bl_fluxrm(kost_flux_zidx,:),roundplf)...
        ./bl_fluxrm(kost_flux_zidx,:),barwidth/3,'facecolor','white');
colornow = fcol; format_fbonoff_bars; % needs colornow,hatchedbars,solidbars:
xlim([0.5 nregs+1.5]); ylim(ylimsnow);
set(ax,'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Region','fontsize',labelfontsize);
ylabel('Relative change [%]','fontsize',labelfontsize);
title(['100-yr change in regional mean flux at ' num2str(fluxdepth+75) ' m'],...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

%---------ROW 2
% MAKE THIS INTO MAP INSTEAD!!!!
% DEPTHS DOWN TO 2000 M, INTERVAL OF 25 M
% ylimsnow = 

% - Zonal mean depth where PSR feedback goes to 0
ax = subplot(223);
plot(lat2, -fbstsd_fluxzm,...
    'color',fcol,'linestyle',fbofflinestyle,...
    'linewidth',fbofflinewidth); hold on;
xlim([-80 80]); % ylim(ylimsnow);
set(ax,'XTick',linspace(-80,80,9),...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Latitude','fontsize',labelfontsize);
ylabel('Depth [m]','fontsize',labelfontsize);
title('Zonal mean depth where PSR feedback goes to 0',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

% - Regional mean depth where PSR feedback goes to 0
barwidth = 0.8;
ax = subplot(224);
solidbars = bar(1:nregs+1, -fbstsd_fluxrm,...
                barwidth, 'facecolor', fcol);
xlim([0.5 nregs+1.5]); % ylim(ylimsnow);
set(ax,'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'TickDir','out','fontsize',tickfontsize);
xlabel('Region','fontsize',labelfontsize);
ylabel('Depth [m]','fontsize',labelfontsize);
title('Regional mean depth where PSR feedback goes to 0',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

print(f, [fig_save_path figname '_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path figname '_PSRfbpaper_final.png'], '-dpng', '-r300');
