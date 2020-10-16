close all;
clear all;
setup_figs;

%-------------------------------
% Define desired circ rate change
%-------------------------------
circfactor = 0.9;
circadj = 1;
figname = 'fig10';

%-------------------------------
% Load grid variables + regions map
%-------------------------------
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
    reg_struct.regnamesabbrev(5:8)]
nregs = length(reg_names);
%figure; pcolor(reg_map); shading flat; colorbar;

%-------------------------------
% Load export model output
%-------------------------------

% ----- Baseline 1*circ, steady-state
bl_fname = [model_output_path 'BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_Kostbetainit_PSDfboff_epsnolims_1.00circ_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_spiralintoeqiternum0.mat'];
load(bl_fname,'output');
bl_ce = output.expCmapSS;
% - Regional + global means
%needs mapnow,nregs,reg_map,area_ocn_only,tot_ocn_area
mapnow = bl_ce; calc_reg_global_means; bl_cerm = mapnowrm;

% ----- PSR fb off, new circ, 100 yr model output
off_fname = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfboffKostbeta_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor)];
load(off_fname,'output');
off_ce = output.expCmapnow(:,:,end);
% - Regional + global means
%needs mapnow,nregs,reg_map,area_ocn_only,tot_ocn_area
mapnow = off_ce; calc_reg_global_means; off_cerm = mapnowrm;

% ----- PSR fb on, one region at a time + GLB, new circ, 100 yr model output
for ionreg=1:nregs
    on_fnames{ionreg} = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocksm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit_%sonly.mat',circfactor,reg_names{ionreg})];
end
on_fnames{nregs+1} = [model_output_path sprintf('BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocksm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactor)];
on_ces = nan(length(lat2),length(lon2),nregs+1);
on_cerms = nan(nregs+1,nregs+1);
for ionreg=1:(nregs+1)
    load(on_fnames{ionreg},'output');
    on_ces(:,:,ionreg) = output.expCmapnow(:,:,end);
    % - Regional + global means
    %needs mapnow,nregs,reg_map,area_ocn_only,tot_ocn_area
    mapnow = on_ces(:,:,ionreg); calc_reg_global_means;
    on_cerms(ionreg,:) = mapnowrm;
end

%-------------------------------
% Calculate export-based fb strengths
%-------------------------------
% - Calculate fb strength
roundpl = -2;
fbsizemap = nan(length(lat2),length(lon2),nregs+1);
fbsizemat = nan(nregs+1,nregs+1);
for ionreg = 1:(nregs+1)
    fbsizemap(:,:,ionreg) = ...
        -roundn(on_ces(:,:,ionreg)-off_ce,roundpl)...
        ./(off_ce-bl_ce);
    mapnow = fbsizemap(:,:,ionreg);
    %needs mapnow,nregs,reg_map,area_ocn_only,tot_ocn_area
    calc_reg_global_means;
    fbsizemat(ionreg,:) = mapnowrm;
end

% - Calculate contribution to global fb strength
% (normalized fb strength)
for irow=1:size(fbsizemat,1)
    fbsizematglobctrb(irow,:) = fbsizemat(irow,:)./fbsizemat(end,:);
end
fbsizematglobctrb(fbsizemat==0) = 0;

%-------------------------------
% Plot figure
%-------------------------------
divcmap = flipud(cbrewer('div','RdYlBu',25,'linear'));
labelfontsize = 10; labelfontwt = 'bold';
titlefontsize = 11; titlefontwt = 'bold';
numfontsize = 9; numfontwt = 'normal';
cbfontsize = 12; cbfontwt = 'normal';
cbticklen = 0.02;

f=figure; set(f,'color','white',...
    'units','inches','position',[1 1 6 10]);

subplot(2,1,1);
imagesc(100*fbsizemat);
caxis([-23 23]);
colormap(divcmap);
cb=colorbar; set(cb,'fontsize',cbfontsize,'fontweight',cbfontwt);
cb.Label.String = '[%]'; cb.TickLength = cbticklen;
set(gca,'XTick',1:size(fbsizemat,2),'XTickLabel',[reg_names {'GLB'}]);
set(gca,'YTick',1:size(fbsizemat,1),'YTickLabel',[reg_names {'GLB'}]);
set(gca,'TickDir','out','FontSize',labelfontsize,'FontWeight',labelfontwt);
title('Regional feedback strengths_ ', ...
    'fontsize',titlefontsize,'fontweight',titlefontwt);
ylabel('Feedback-on region');
xlabel('Affected region');
mat = 100*fbsizemat;
matl = mat(:); textStrings = num2str(matl,'%0.2f');
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
textStrings(abs(matl)>1) = cellstr(num2str(matl(abs(matl)>1),'%0.1f'));
textStrings(abs(matl)<0.005)={'0'};
[x,y] = meshgrid(1:size(mat,2));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'fontsize',numfontsize,'fontweight',numfontwt);

subplot(2,1,2);
imagesc(100*fbsizematglobctrb);
caxis([-100 100]);
colormap(divcmap);
cb=colorbar; set(cb,'fontsize',cbfontsize,'fontweight',cbfontwt);
cb.Label.String = '[%]'; cb.TickLength = cbticklen;
set(gca,'XTick',1:size(fbsizemat,2),'XTickLabel',[reg_names {'GLB'}]);
set(gca,'YTick',1:size(fbsizemat,1),'YTickLabel',[reg_names {'GLB'}]);
set(gca,'TickDir','out','FontSize',labelfontsize,'FontWeight',labelfontwt);
title('Normalized regional feedback strengths_ ', ...
    'fontsize',titlefontsize,'fontweight',titlefontwt);
ylabel('Feedback-on region');
xlabel('Affected region');
mat = 100*fbsizematglobctrb; matl = mat(:);
textStrings = num2str(matl,'%0.2f');
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
textStrings(abs(matl)>1) = cellstr(num2str(matl(abs(matl)>1),'%0.1f'));
textStrings(abs(matl)<0.005)={'0'};
textStrings(mat(:)==100)={'100'};
[x,y] = meshgrid(1:size(mat,2));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'fontsize',numfontsize,'fontweight',numfontwt);

print(f, [fig_save_path figname '_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path figname '_PSRfbpaper_final.png'], '-dpng', '-r300');
