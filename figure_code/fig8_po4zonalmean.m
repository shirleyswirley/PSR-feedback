close all;
clear all;
setup_figs;

%-------------------------------
% Define desired circ rate change
%-------------------------------
circfactor = 0.9;
figname = 'fig8';

%-------------------------------
% Load grid variables + model output
%-------------------------------
% - Load grid variables + regions map
reg_struct = load([data_path ...
    'PRiSM_regions_2deg.mat']);
M3d = reg_struct.M3d; gridd = reg_struct.gridd;
lon2 = gridd.xt; lat2 = gridd.yt;
depth = gridd.zt; [~,deepdidx]=min(abs(depth-1500));
iocn = find(M3d==1); isurf = find(M3d(:,:,1)==1);

%----------Load po4 maps
% - Baseline 1*circ, steady-state
bl_fname = [model_output_path 'BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_Kostbetainit_PSDfboff_epsnolims_1.00circ_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_spiralintoeqiternum0.mat'];
load(bl_fname,'output');
bl_p = M3d*NaN;
bl_p(iocn) = output.po4;
bl_pzm = squeeze(nanmean(bl_p,2)); % zonal mean

% - PSR fb off,
% new circ, 100 yr model output
off_fname = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfboffKostbeta_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor)];
load(off_fname,'output');
off_p = M3d*NaN;
off_p(iocn) = output.po4end;
off_pzm = squeeze(nanmean(off_p,2)); % zonal mean

% - PSR fb on in all regions,
% new circ, 100 yr model output 
on_fname = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocksm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor)];
load(on_fname,'output');
on_p = M3d*NaN;
on_p(iocn) = output.po4end;
on_pzm = squeeze(nanmean(on_p,2)); % zonal mean

%-------------------------------
% Plot figure
%-------------------------------
seqcmap = cbrewer('seq','YlGnBu',20,'linear');
divcmap = flipud(cbrewer('div','RdBu',21,'linear'));
titlesnow = {'Baseline zonal mean PO_4',...
    '100-yr PO_4 change WITHOUT PSR feedback',...
    '100-yr PO_4 change WITH PSR feedback',...
    'PO_4 WITH feedback - WITHOUT feedback'};
unitsnow = {'[mmol m^{-3}]','[mmol m^{-3}]',...
    '[mmol m^{-3}]','[mmol m^{-3}]'};
didxsnow = 1:deepdidx;
plotsnow = nan(length(depth(didxsnow)),length(lat2),4); 
plotsnow(:,:,1) = bl_pzm(:,didxsnow)';
plotsnow(:,:,2) = (off_pzm(:,didxsnow)-bl_pzm(:,didxsnow))';
plotsnow(:,:,3) = (on_pzm(:,didxsnow)-bl_pzm(:,didxsnow))';
plotsnow(:,:,4) = (on_pzm(:,didxsnow)-off_pzm(:,didxsnow))';
f=figure; set(f,'color','white',...
    'units','inches','position',[0.5 0.5 5 11]);
for isp=1:4
    ax=subplot(4,1,isp);
    pcolor(lat2,depth(didxsnow),plotsnow(:,:,isp));
    shading flat; cb = colorbar; cb.TickLength = 0.05;
    cb.Label.String = unitsnow(isp); 
    xlim([-80 80]); xticks(-80:20:80);
    yticks(250:250:1500); yticklabels(); grid on;
    set(gca,'YDir','reverse','TickLength',...
        [0.03 0.035],'layer','top');
    xlabel('Latitude'); ylabel('Depth [m]');
    if isp==1
        colormap(ax,seqcmap);
    elseif isp==2 | isp==3
        cmin = -0.06; cmax = 0.06;
        caxis([cmin cmax]); colormap(ax,divcmap);
        mapvarnow = plotsnow(:,:,isp);
        cbnumticks = 5;
        extend_cbar_ticklabels; % needs cmin,cmax,mapvarnow,cb,cbnumticks
    elseif isp==4
        caxis([-max(max(abs(plotsnow(:,:,isp)))) ...
                max(max(abs(plotsnow(:,:,isp))))]);
        colormap(ax,divcmap);
    end
    title(titlesnow{isp});
end

print(f, [fig_save_path figname '_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path figname '_PSRfbpaper_final.png'], '-dpng', '-r300');
