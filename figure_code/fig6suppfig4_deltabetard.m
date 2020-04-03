close all;
clear all;
setup_figs;

%-------------------------------
% Define desired circ rate change
%-------------------------------
% - Main figure 6:
circfactor = 0.9;
figname = 'fig6';
circadj = 1;
betachylab = '\Delta\beta';
rdchdir = 'shoaling'; % remin depth change dir

% - Supplementary figure 4:
%circfactor = 1.1;
%figname = 'suppfig4';
%circadj = -1;
%betachylab = '\beta decrease';
%rdchdir = 'deepening'; % remin depth change dir

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

%----------Load beta maps
% - PSD fb on in all regions,
% new circ, 100 yr model output 
on_fname = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocksm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin1yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor)];
load(on_fname,'output');
on_beta_map = output.betamapnow(:,:,end);

% - PSD fb on - 1stdev in all regions,
% new circ, 100 yr model output
onms_fname = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdbelowm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor)];
load(onms_fname,'output');
onms_beta_map = output.betamapnow(:,:,end);

% - PSD fb on + 1stdev in all regions,
% new circ, 100 yr model output
onps_fname = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdabovem1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor)];
load(onps_fname,'output');
onps_beta_map = output.betamapnow(:,:,end);

% - Kost beta
% (applies to baseline 1*circ and
% any circ PSD fb off runs)
load([data_path 'Kost_beta_2deg.mat'],'beta_clim');
kost_beta_map = mean(beta_clim,3)';
kost_beta_map(M3d(:,:,1)==0)=NaN;

%----------Load/create remin depth maps
% - PSD fb on in all regions,
% new circ, 100 yr model output 
beta_name = ['kostfbon_100yr' sprintf('%.2f',circfactor) 'circ'];
load_rd_map=1;
if load_rd_map==1
    load([data_path 'remin_depth_map_' beta_name 'beta.mat'],...
        'remin_depth_map');
    on_rd_map = remin_depth_map;
else
    var_des = ...
        ['remin depths in [m] calc using PRiSM from 0.9 circ fb on, 100-yr beta map; model output used was' on_fname];
    on_rd_map = calc_remin_depth_map(...
        lon2,lat2,on_beta_map,beta_name,var_des);
end

% - PSD fb on - 1stdev in all regions,
% new circ, 100 yr model output
beta_name = ['kostfbonminusstdev_100yr' sprintf('%.2f',circfactor) 'circ'];
load_rd_map=1;
if load_rd_map==1
    load([data_path 'remin_depth_map_' beta_name 'beta.mat'],...
        'remin_depth_map');
    onms_rd_map = remin_depth_map;
else
    var_des = ...
        ['remin depths in [m] calc using PRiSM from 0.9 circ fb on minus 1stdev, 100-yr beta map; model output used was' onms_fname];
    onms_rd_map = calc_remin_depth_map(...
        lon2,lat2,onms_beta_map,beta_name,var_des);
end

% - PSD fb on + 1stdev in all regions,
% new circ, 100 yr model output
beta_name = ['kostfbonplusstdev_100yr' sprintf('%.2f',circfactor) 'circ'];
load_rd_map=1;
if load_rd_map==1
    load([data_path 'remin_depth_map_' beta_name 'beta.mat'],...
        'remin_depth_map');
    onps_rd_map = remin_depth_map;
else
    var_des = ...
        ['remin depths in [m] calc using PRiSM from 0.9 circ fb on plus 1stdev, 100-yr beta map; model output used was' onps_fname];
    onps_rd_map = calc_remin_depth_map(...
        lon2,lat2,onps_beta_map,beta_name,var_des);
end

% - Kost beta
% (applies to baseline 1*circ and
% any circ PSD fb off runs)
beta_name = 'kost';
load_rd_map=1;
if load_rd_map==1
    load([data_path 'remin_depth_map_' beta_name 'beta.mat'],...
        'remin_depth_map');
    kost_rd_map = remin_depth_map;
else
    var_des = ...
        'remin depths in [m] calc using PRiSM from annual mean Kost beta map';
    kost_rd_map = calc_remin_depth_map(...
        lon2,lat2,kost_beta_map,beta_name,var_des);
end

%-------------------------------
% Calculate global, zonal, and
% regional mean beta and remin depths
%-------------------------------
% - Calculate global mean beta
kost_beta_gmean = ...
    nansum(nansum(kost_beta_map.*area_ocn_only))...
    ./nansum(nansum(area_ocn_only));
on_beta_gmean = ...
    nansum(nansum(on_beta_map.*area_ocn_only))...
    ./nansum(nansum(area_ocn_only));
onms_beta_gmean = ...
    nansum(nansum(onms_beta_map.*area_ocn_only))...
    ./nansum(nansum(area_ocn_only));
onps_beta_gmean = ...
    nansum(nansum(onps_beta_map.*area_ocn_only))...
    ./nansum(nansum(area_ocn_only));

% - Calculate regional mean betas
kost_beta_rmeans = nan(1,nregs);
for ireg=1:nregs
    betanow = kost_beta_map(reg_map==ireg);
    area_ocn_onlynow = area_ocn_only(reg_map==ireg);
    kost_beta_rmeans(ireg) = nansum(nansum(betanow.*area_ocn_onlynow))...
        ./nansum(nansum(area_ocn_onlynow));
end
on_beta_rmeans = nan(1,nregs);
for ireg=1:nregs
    betanow = on_beta_map(reg_map==ireg);
    area_ocn_onlynow = area_ocn_only(reg_map==ireg);
    on_beta_rmeans(ireg) = nansum(nansum(betanow.*area_ocn_onlynow))...
        ./nansum(nansum(area_ocn_onlynow));
end
onms_beta_rmeans = nan(1,nregs);
for ireg=1:nregs
    betanow = onms_beta_map(reg_map==ireg);
    area_ocn_onlynow = area_ocn_only(reg_map==ireg);
    onms_beta_rmeans(ireg) = nansum(nansum(betanow.*area_ocn_onlynow))...
        ./nansum(nansum(area_ocn_onlynow));
end
onps_beta_rmeans = nan(1,nregs);
for ireg=1:nregs
    betanow = onps_beta_map(reg_map==ireg);
    area_ocn_onlynow = area_ocn_only(reg_map==ireg);
    onps_beta_rmeans(ireg) = nansum(nansum(betanow.*area_ocn_onlynow))...
        ./nansum(nansum(area_ocn_onlynow));
end

% - Calculate zonal mean betas
kost_beta_z = ...
    nansum(kost_beta_map.*area_ocn_only,2)...
    ./area_oo_zsum;
on_beta_z = ...
    nansum(on_beta_map.*area_ocn_only,2)...
    ./area_oo_zsum;
onms_beta_z = ...
    nansum(onms_beta_map.*area_ocn_only,2)...
    ./area_oo_zsum;
onps_beta_z = ...
    nansum(onps_beta_map.*area_ocn_only,2)...
    ./area_oo_zsum;

% - Calculate global mean rd
kost_rd_gmean = ...
    nansum(nansum(kost_rd_map.*area_ocn_only))...
    ./nansum(nansum(area_ocn_only));
on_rd_gmean = ...
    nansum(nansum(on_rd_map.*area_ocn_only))...
    ./nansum(nansum(area_ocn_only));
onms_rd_gmean = ...
    nansum(nansum(onms_rd_map.*area_ocn_only))...
    ./nansum(nansum(area_ocn_only));
onps_rd_gmean = ...
    nansum(nansum(onps_rd_map.*area_ocn_only))...
    ./nansum(nansum(area_ocn_only));

% - Calculate regional mean rds
kost_rd_rmeans = nan(1,nregs);
for ireg=1:nregs
    rdnow = kost_rd_map(reg_map==ireg);
    area_ocn_onlynow = area_ocn_only(reg_map==ireg);
    kost_rd_rmeans(ireg) = nansum(nansum(rdnow.*area_ocn_onlynow))...
        ./nansum(nansum(area_ocn_onlynow));
end
on_rd_rmeans = nan(1,nregs);
for ireg=1:nregs
    rdnow = on_rd_map(reg_map==ireg);
    area_ocn_onlynow = area_ocn_only(reg_map==ireg);
    on_rd_rmeans(ireg) = nansum(nansum(rdnow.*area_ocn_onlynow))...
        ./nansum(nansum(area_ocn_onlynow));
end
onms_rd_rmeans = nan(1,nregs);
for ireg=1:nregs
    rdnow = onms_rd_map(reg_map==ireg);
    area_ocn_onlynow = area_ocn_only(reg_map==ireg);
    onms_rd_rmeans(ireg) = nansum(nansum(rdnow.*area_ocn_onlynow))...
        ./nansum(nansum(area_ocn_onlynow));
end
onps_rd_rmeans = nan(1,nregs);
for ireg=1:nregs
    rdnow = onps_rd_map(reg_map==ireg);
    area_ocn_onlynow = area_ocn_only(reg_map==ireg);
    onps_rd_rmeans(ireg) = nansum(nansum(rdnow.*area_ocn_onlynow))...
        ./nansum(nansum(area_ocn_onlynow));
end

% - Calculate zonal mean rds
kost_rd_z = ...
    nansum(kost_rd_map.*area_ocn_only,2)...
    ./area_oo_zsum;
on_rd_z = ...
    nansum(on_rd_map.*area_ocn_only,2)...
    ./area_oo_zsum;
onms_rd_z = ...
    nansum(onms_rd_map.*area_ocn_only,2)...
    ./area_oo_zsum;
onps_rd_z = ...
    nansum(onps_rd_map.*area_ocn_only,2)...
    ./area_oo_zsum;

%-------------------------------
% Plot figure
%-------------------------------
%ploterrbars=0;
%ebwidth = 100; % 1/ebwidth is the width, technically
tickfontsize = 9;
labelfontsize = 10; labelfontwt = 'normal';
titlefontsize = 11; titlefontwt = 'bold';
barwidth = 0.3; offset = barwidth/2;
betacol = [0 0.55 0]; rdcol = [0.73 0 0.73];

f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 8.5 5]);

betamin = 3.25; betamax = 5.25;
betaticks = linspace(3,5,5);
rdmin = 0; rdmax = 1700;
rdticks = linspace(0,1600,5);
rdticklabs = compose('%g',rdticks/100);

subplot(221);
[ax,hbeta,hrd]=plotyy(...
    lat2,kost_beta_z,...
    lat2,kost_rd_z);
set(hbeta,'Color',betacol,'LineStyle','-','LineWidth',2);
set(hrd,'Color',rdcol,'LineStyle','-','Linewidth',2);
set(ax(1),'YColor',betacol,'XLim',[-80 80],...
    'XTick',linspace(-80,80,9),'box','off',...
    'YLim',[betamin betamax],'YTick',betaticks,...
    'fontsize',tickfontsize);
set(ax(2),'YColor',rdcol,'XLim',[-80 80],...
    'YTickLabels',rdticklabs,...
    'YLim',[rdmin rdmax],'YTick',rdticks,...
    'fontsize',tickfontsize);
rl=refline(ax(1),0,ax(1).YLim(2)); rl.Color='k';
ylabel(ax(1),'\beta',...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
ylabel(ax(2),'Remin depth [x100 m]',...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
set(ax,'TickDir','out');
xlabel('Latitude','fontsize',labelfontsize);
title(['Baseline_ zonal mean '...
    sprintf('%s{%f %f %f}','\color[rgb]',betacol) '\beta' ...
    '\color{black} & ' ...
    sprintf('%s{%f %f %f}','\color[rgb]',rdcol) 'remin depth'],...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

subplot(222);
[ax,hbeta,hrd]=plotyy(...
    (1:nregs+1)-offset,[kost_beta_rmeans kost_beta_gmean],...
    (1:nregs+1)+offset,[kost_rd_rmeans kost_rd_gmean],...
    @(x,y) bar(x,y,barwidth), @(x,y) bar(x,y,barwidth));
set(hbeta,'FaceColor',betacol);
set(hrd,'FaceColor',rdcol);
set(ax(1),'YColor',betacol,'YLim',[betamin betamax],'box','off',...
    'YTick',betaticks,'XLim',[0.5 nregs+1.5],...
    'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'fontsize',tickfontsize);
set(ax(2),'YColor',rdcol,'YLim',[rdmin rdmax],...
    'YTick',rdticks,'YTickLabels',rdticklabs,...
    'XLim',[0.5 nregs+1.5],'XTick',[],...
    'fontsize',tickfontsize);
rl1=refline(ax(1),0,ax(1).YLim(2)); rl1.Color='k'; % adds box top edge
rl2=refline(ax(2),0,ax(2).YLim(1)); rl2.Color='k'; % make x-axis black
ylabel(ax(1),'\beta','FontSize',labelfontsize,...
    'FontWeight',labelfontwt);
ylabel(ax(2),'Remin depth [x100 m]',...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
set(ax,'TickDir','out');
xlabel('Region','fontsize',labelfontsize);
title(['Baseline_ regional mean '...
    sprintf('%s{%f %f %f}','\color[rgb]',betacol) '\beta' ...
    '\color{black} & ' ...
    sprintf('%s{%f %f %f}','\color[rgb]',rdcol) 'remin depth'],...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

if circfactor==0.9
    betamax = 0.08; betalsn = 10;
    betanoticklabs = 1:2:9;
elseif circfactor==1.1
    betamax = 0.09; betalsn = 11;
    betanoticklabs = 4:2:10;
end 
betamin = -0.01;
betaticks = linspace(betamin,betamax,betalsn);
betaticks(2) = 0; % wasn't perfectly zero for some reason before
betaticklabs = compose('%g',betaticks);
for ilab=betanoticklabs
    betaticklabs{ilab} = ' ';
end
if circfactor==0.9
    rdmax = 80; rdlsn = 10;
    rdnoticklabs = 1:2:9;
elseif circfactor==1.1
    rdmax = 90; rdlsn = 11;
    rdnoticklabs = 4:2:10;
end
rdmin = -10;
rdticks = linspace(rdmin,rdmax,rdlsn);
rdticklabs = compose('%g',rdticks);
for ilab=rdnoticklabs
    rdticklabs{ilab} = ' ';
end

subplot(223);
[ax,hbeta,hrd]=plotyy(...
    lat2,circadj*(on_beta_z-kost_beta_z),...
    lat2,-circadj*(on_rd_z-kost_rd_z));
set(hbeta,'Color',betacol,'LineStyle','-','LineWidth',2);
set(hrd,'Color',rdcol,'LineStyle','-','Linewidth',2);
set(ax(1),'YColor',betacol,'XLim',[-80 80],...
    'YTickLabels',betaticklabs,...
    'XTick',linspace(-80,80,9),'box','off',...
    'YLim',[betamin betamax],'YTick',betaticks,...
    'fontsize',tickfontsize);
set(ax(2),'YColor',rdcol,'XLim',[-80 80],...
    'YTickLabels',rdticklabs,...
    'YLim',[rdmin rdmax],'YTick',rdticks,...
    'fontsize',tickfontsize);
rl1=refline(ax(1),0,ax(1).YLim(2)); rl1.Color='k';
rl2=refline(ax(1),0,0); rl2.Color='k';
ylabel(ax(1),betachylab,...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
ylabel(ax(2),['Remin depth ' rdchdir ' [m]'],...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
set(ax,'TickDir','out');
xlabel('Latitude','fontsize',labelfontsize);
title('100-yr_ zonal change w/ PSR feedback',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

subplot(224);
[ax,hbeta,hrd]=plotyy(...
    (1:nregs+1)-offset,circadj*([on_beta_rmeans on_beta_gmean]-[kost_beta_rmeans kost_beta_gmean]),...
    (1:nregs+1)+offset,-circadj*([on_rd_rmeans on_rd_gmean]-[kost_rd_rmeans kost_rd_gmean]),...
    @(x,y) bar(x,y,barwidth), @(x,y) bar(x,y,barwidth));
set(hbeta,'FaceColor',betacol);
set(hrd,'FaceColor',rdcol);
set(ax(1),'YColor',betacol,'YLim',[betamin betamax],'box','off',...
    'YTickLabels',betaticklabs,...
    'YTick',betaticks,'XLim',[0.5 nregs+1.5],...
    'XTick',1:nregs+1,'XTickLabel',[reg_names {'GLB'}],...
    'fontsize',tickfontsize);
set(ax(2),'YColor',rdcol,'YLim',[rdmin rdmax],...
    'YTickLabels',rdticklabs,...
    'YTick',rdticks,...
    'XLim',[0.5 nregs+1.5],'XTick',[],...
    'fontsize',tickfontsize);
rl1=refline(ax(1),0,ax(1).YLim(2)); rl1.Color='k'; % adds box top edge
rl2=refline(ax(2),0,0); rl2.Color='k'; % make x-axis black
ylabel(ax(1),betachylab,'FontSize',labelfontsize,...
    'FontWeight',labelfontwt);
ylabel(ax(2),['Remin depth ' rdchdir ' [m]'],...
    'FontSize',labelfontsize,'FontWeight',labelfontwt);
set(ax,'TickDir','out');
xlabel('Region','fontsize',labelfontsize);
title('100-yr_ regional change w/ PSR feedback',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

disp('pause - manually move lhs subplots over to the left moreby 12 or 20 left arrow clicks')
pause

print(f, [fig_save_path figname '_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path figname '_PSRfbpaper_final.png'], '-dpng', '-r300');

disp('glob mean beta, Kost:');
kost_beta_gmean
disp(['glob mean beta, 100 yr fb on - Kost, ' num2str(circfactor) ' circ:']);
on_beta_gmean

disp('glob mean rd, Kost:');
kost_rd_gmean
disp(['glob mean rd, 100 yr fb on - Kost, ' num2str(circfactor) ' circ:']);
on_rd_gmean

disp('reg mean beta, Kost:');
kost_beta_rmeans
disp(['reg mean beta, 100 yr fb on - Kost, ' num2str(circfactor) ' circ:']);
[on_beta_rmeans on_beta_gmean]-[kost_beta_rmeans kost_beta_gmean]

disp('reg mean rd, Kost:');
kost_rd_rmeans
disp(['reg mean rd, 100 yr fb on - Kost, ' num2str(circfactor) ' circ:']);
[on_rd_rmeans on_rd_gmean]-[kost_rd_rmeans kost_rd_gmean]
