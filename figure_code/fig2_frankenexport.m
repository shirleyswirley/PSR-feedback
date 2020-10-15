close all;
clear all;
setup_figs;

%-------------------------------
% Load frankenstein C export map + export regions map
%-------------------------------
franken = load([data_path 'franken_export_1degdiffgrid.mat']);
lon1 = franken.grid.xt; lat1 = franken.grid.yt;
fcexp_map = franken.Cfrank;
reg_map = franken.R2d;
reg_names = franken.regnames;
reg_names = [reg_names,'Indian','Arctic','MedSeaPlus'];
reg_names_sh = {'AAZ','SAZ','STA','STP','ETA','ETP','NA','NP','IND','ARC','MED'};
[lon1grid,lat1grid] = meshgrid(lon1,lat1);
reg_map(reg_map==0 & lat1grid<30) = 9; % Indian Ocean region
reg_map(reg_map==0 & lat1grid>60) = 10; % Arctic Ocean region
reg_map(reg_map==0) = 11; % Mediterannean Sea + other seas

%-------------------------------
% Plot figure
%-------------------------------
% - Define plot params
cmap=flipud(cbrewer('seq','YlGnBu',20,'linear'));
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
linewidth = 2; % for regional contour lines
linecolor = [0 0 0]; % for regional contour lines
labelfontsize = 10;
titlefontsize = 12;
titlefontwt = 'bold';
cbticklen = 0.03;
cbnumticks = 4;
cmin = -0.5; cmax = 1;
plotregcontours=1;

% - Plot annual mean frankenstein carbon export
f=figure;
set(f,'color','white','units','inches',...
    'position',[0.5 0.5 4 2.5],'resize','off');
m_proj(mapproj,'lon',[0 360],'lat',[min(lat1) max(lat1)])
m_pcolor(lon1,lat1,log10(fcexp_map));
m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
caxis([cmin cmax]); colormap(cmap); shading flat;
cb = colorbar; cb.TickLength = cbticklen; cb.FontSize=labelfontsize;
cb.Ticks = [-0.5, 0, 0.5, 1]; cb.TickLabels = {'<10^{-0.5}', '10^0', '10^{0.5}', '>10^1'};
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
title('Export [molC m^{-2} yr^{-1}]',...
    'fontsize',titlefontsize,'fontweight',titlefontwt);

print(f, [fig_save_path 'fig3_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'fig3_PSRfbpaper_final.png'], '-dpng', '-r300');
