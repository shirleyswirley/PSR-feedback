close all;
clear all;
setup_figs;

%-------------------------------
% Load grid variables
%-------------------------------
reg_struct = load([data_path ...
    'PRiSM_regions_2deg.mat']);
M3d = reg_struct.M3d; gridd = reg_struct.gridd;
lon2 = gridd.xt; lat2 = gridd.yt;
isurf = find(M3d(:,:,1)==1);
area_ocn_only = nan(size(gridd.Areat)); % m^2
area_ocn_only(isurf) = gridd.Areat(isurf);
tot_ocn_area = sum(gridd.Areat(isurf)); % [m^2]

%-------------------------------
% Define desired circ rates
%-------------------------------
% circulation factors we want to plot
circfactorarray = [0.9 1.1];
numcircfactors = length(circfactorarray);

%-------------------------------
% Load export model output
%-------------------------------

% ----- Baseline 1*circ
bl_fname = [model_output_path 'BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_Kostbetainit_PSDfboff_epsnolims_1.00circ_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_spiralintoeqiternum0.mat'];
load(bl_fname,'output');
bl_ssce = output.expCmapSS;
bl_sstotce = nansum(nansum(bl_ssce.*area_ocn_only));

% ----- PSR fb off, new circs
for icirc=1:numcircfactors
    off_fnames{icirc} = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfboffKostbeta_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_1ov1000dt_ts100yrs_prevNMinit.mat',circfactorarray(icirc))];
    load(off_fnames{icirc},'output');
    off_ssce(:,:,icirc) = output.expCmapnow(:,:,end);
    off_sstotce(icirc) = nansum(nansum(off_ssce(:,:,icirc).*area_ocn_only)); 
    off_otce(:,icirc) = output.totexpC; % tot carbon export over time
    off_time(:,icirc) = output.time;
end

% ----- PSR fb on in all regions, new circs
for icirc = 1:numcircfactors
    on_fnames{icirc} = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocksm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin1yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactorarray(icirc))];
    load(on_fnames{icirc},'output');
    on_ssce(:,:,icirc) = output.expCmapnow(:,:,end);
    on_sstotce(icirc) = nansum(nansum(on_ssce(:,:,icirc).*area_ocn_only));
    on_otce(:,icirc) = output.totexpC; % tot carbon export over time
    on_time(:,icirc) = output.time;
end

% ----- PSR fb on in all regions, new circs, beta vs. export - 1stdev
for icirc = 1:numcircfactors
    onms_fnames{icirc} = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdbelowm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactorarray(icirc))];
    load(onms_fnames{icirc},'output');
    onms_ssce(:,:,icirc) = output.expCmapnow(:,:,end);
    onms_sstotce(icirc) = nansum(nansum(onms_ssce(:,:,icirc).*area_ocn_only));
    onms_otce(:,icirc) = output.totexpC; % tot carbon export over time
    onms_time(:,icirc) = output.time;
end

% ----- PSR fb on in all regions, new circs, beta vs. export + 1stdev
for icirc = 1:numcircfactors
    onps_fnames{icirc} = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdabovem1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactorarray(icirc))];
    load(onps_fnames{icirc},'output');
    onps_ssce(:,:,icirc) = output.expCmapnow(:,:,end);
    onps_sstotce(icirc) = nansum(nansum(onps_ssce(:,:,icirc).*area_ocn_only));
    onps_otce(:,icirc) = output.totexpC; % tot carbon export over time
    onps_time(:,icirc) = output.time;
end

disp('pause')
pause

%-------------------------------
% Calculate SS abs and rel changes
% in global mean export from baseline case
%-------------------------------
off_acce = nan(numcircfactors,1);
for icirc=1:numcircfactors
    off_acce(icirc)=(off_sstotce(icirc)-bl_sstotce)./tot_ocn_area;
    off_rcce(icirc)=(off_sstotce(icirc)-bl_sstotce)./bl_sstotce;
end

on_acce = nan(numcircfactors,1);
for icirc=1:numcircfactors
    on_acce(icirc)=(on_sstotce(icirc)-bl_sstotce)./tot_ocn_area;
    on_rcce(icirc)=(on_sstotce(icirc)-bl_sstotce)./bl_sstotce;
end

onms_acce = nan(numcircfactors,1);
for icirc=1:numcircfactors
    onms_acce(icirc)=(onms_sstotce(icirc)-bl_sstotce)./tot_ocn_area;
    onms_rcce(icirc)=(onms_sstotce(icirc)-bl_sstotce)./bl_sstotce;
end

onps_acce = nan(numcircfactors,1);
for icirc=1:numcircfactors
    onps_acce(icirc)=(onps_sstotce(icirc)-bl_sstotce)./tot_ocn_area;
    onps_rcce(icirc)=(onps_sstotce(icirc)-bl_sstotce)./bl_sstotce;
end

%-------------------------------
% Plot figure
%-------------------------------
f=figure;
set(f,'color','white','units',...
    'inches','position',[1 1 9 5]);

fbofflinestyle='-'; fbonlinestyle='--';
linewidth=3; linecol='b';
labelfontsize = 12; labelfontwt = 'bold';

% - Abs change over time in
% global mean export from baseline case
subplot(121);
p(1)=plot([-10 0],[0 0],'k','linewidth',linewidth); hold on;
for icirc=1:numcircfactors
    % - PSR fb off runs
    plot([0 off_time(:,icirc)'],...
        ([bl_sstotce off_otce(:,icirc)']-bl_sstotce)./tot_ocn_area,...
        'color','b','linestyle',fbofflinestyle,'linewidth',linewidth);
    % - PSR fb on runs
    plot([0 on_time(:,icirc)'],...
        ([bl_sstotce on_otce(:,icirc)']-bl_sstotce)./tot_ocn_area,...
        'color','b','linestyle',fbonlinestyle,'linewidth',linewidth);
end
ylabel('Global mean export change [molC m^{-2} yr^{-1}]');
xlabel('Years after circulation change');
xlim([-10 100]); ylim([-0.3 0.3]);
set(gca,'TickLength',[0.02 0.01],'YTick',-0.3:0.1:0.3,...
        'fontsize',labelfontsize,'fontweight',labelfontwt,...
        'YTick',-0.3:0.1:0.3);
%legend(' ',' ',' ',' ',' '); % just to get the legend elems on the fig
%legend(p,{'Baseline','Feedback off','Feedback on'});

% - Barplot of SS abs change in
% global mean export from baseline case
sp = subplot(122);
% move subplot to left
pos = get(sp, 'Position');
posnew = pos; posnew(1) = posnew(1) - 0.1; set(sp,'Position',posnew);
hold on;
% for some reason, gotta run this twice to get the hatching patterns to look correct
for iweird = 1:2
    for icirc = 1:numcircfactors
        % - PSR fb off runs
        barcolor=[0 0 0.8];
        h1 = bar(1,off_acce(icirc),'facecolor',barcolor);hold on;
        % - PSR fb on runs
        barcolor='w';
        h1 = bar(1,on_acce(icirc),'facecolor',barcolor);
        hatchfill2(h1,'hatchstyle','single','hatchangle',45,...
            'hatchdensity',15,'hatchcolor',[0 0 0.8],'hatchlinewidth',2);
        eb = errorbar(1,on_acce(icirc),...
             abs(on_acce(icirc)-onms_acce(icirc)),... neg
             abs(on_acce(icirc)-onps_acce(icirc)),... pos
            'color','k','linewidth',2);
    end
end 
ylim([-0.3 0.3]);
set(gca,'XColor','None','YColor','None');
plot(xlim,[0 0],'-w'); % covers pesky black line at zero
%legend(' ',' ',' ',' ',' ',' ',' '); % just to get the legend elems on the fig

print(f, [fig_save_path 'fig5_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'fig5_PSRfbpaper_final.png'], '-dpng', '-r300');

disp('SS abs export change from baseline w/ PSR fb off:');
off_acce
disp('SS rel export change from baseline w/ PSR fb off:');
off_rcce
disp('SS abs export change from baseline w/ PSR fb on:');
on_acce
disp('SS rel export change from baseline w/ PSR fb on:');
on_rcce
disp('SS abs export change from baseline w/ PSR fb on, w/ beta vs. export - 1 stdev:');
onms_acce
disp('SS rel export change from baseline w/ PSR fb on, w/ beta vs. export - 1 stdev:');
onms_rcce
disp('SS abs export change from baseline w/ PSR fb on, w/ beta vs. export + 1 stdev:');
onps_acce
disp('SS rel export change from baseline w/ PSR fb on, w/ beta vs. export + 1 stdev:');
onps_rcce
