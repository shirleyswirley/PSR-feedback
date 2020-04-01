runloadsection=1;

if runloadsection==1

clear all;
close all;
addpath(genpath('/jetwork2/cram/matlabtools/m_map'));
addpath(genpath('/graid1/shirlleu/matlabroutines'));

% - Define the beta, eratio, and npp algs you want to use
betanamearray = {'Kost','GuidiallUVP'};
erationamearray = {'Laws2000','D2005PP','Laws2011D'};
eratioshortnamearray = {'L2000','D2005','L2011'};
NPPnamearray = {'VGPM','VGPME','CbPM'};
betan = length(betanamearray);
eration = length(erationamearray);
nppn = length(NPPnamearray);
timen = 160; % num of months in time series

%--------------------------------------------------
% Load and regrid beta vs. export
%--------------------------------------------------
% Remember that this grid is ever so slightly diff from the Kost and Guidi beta grids,
% but it's totally ok b/c both are just lower res from orig sat data + close enough
load('/graid1/shirlleu/NPPslopeproj/matfiles/globalexportmaps_spatvareratioandezdepth/expCann_D2005PPfromphysatchl1.mat','lat1','lon1');

% - Temporally regressed beta vs. normalized export maps
% Initialize beta vs. export filenames
fnamebase = ['/graid1/shirlleu/NPPslopeproj/matfiles/globalbetavsexporttempregressmaps/'];
% Initialize matrix to hold all beta vs. export maps
bvse = nan(length(lat1),length(lon1),betan,eration,nppn);
% Load beta vs. export maps
for ibeta = 1:betan
    for ieratio = 1:eration
        for inpp = 1:nppn
            load([fnamebase NPPnamearray{inpp} 'based/' betanamearray{ibeta} 'betavs' erationamearray{ieratio} 'exportn1tempregmo.mat'],'betavsetempregresscoeffmapzeros');
            % GET RID OF CRAZY OUTLIERS
            betavsetempregresscoeffmapzeros(abs(betavsetempregresscoeffmapzeros)>=5)=nan;
            length(find(abs(betavsetempregresscoeffmapzeros)>=5)) % displays how many outliers there were in each beta vs export map
            bvse(:,:,ibeta,ieratio,inpp) = betavsetempregresscoeffmapzeros;
        end
    end
end

% - Put b vs. e on PRiSM grid so we can use those approximate areas to calculate area-wted reg means
load('/graid1/shirlleu/PRiSM_GCM/sran+burial+ballast/newmatfiles/newtonmethod/BGRID_t+s+14c+cfc11_REDI_sol1_sigma0.10_pop2dop_ann_np3_PSDfboff_unifbeta_1.0circ.mat','output');
M3d = output.M3d;
iocn = find(M3d==1);
isurf = find(M3d(:,:,1)==1);
clear output;
pregs = load('/graid1/shirlleu/PRiSM_GCM/data/PRiSMregions.mat');
pgrid = pregs.gridd;
regmappgrid = pregs.R2d;
regnames = pregs.regnames; % pregs.regnames(1:9)
regnamesro = [regnames(1:4) regnames(9) regnames(5:8) regnames(10:end)];
regnamesabbrev = pregs.regnamesabbrev;
regnamesabbrevro = [regnamesabbrev(1:4) regnamesabbrev(9) regnamesabbrev(5:8) regnamesabbrev(10:end)];
Areatocnonly = nan(size(pgrid.Areat)); Areatocnonly(isurf) = pgrid.Areat(isurf);
Areatzonalsum = nansum(Areatocnonly,2);
lon1p = [lon1(181:end)' lon1(1:180)'+360];
[lon1pgrid,lat1grid] = meshgrid(lon1p,lat1);
bvsepgrid1 = [bvse(:,181:end,:,:,:) bvse(:,1:180,:,:,:)];
bvsepgrid = nan(length(pgrid.yt),length(pgrid.xt),betan,eration,nppn);
for ibeta = 1:betan
    for ieratio = 1:eration
        for inpp = 1:nppn
            bvsetemp = bvsepgrid1(:,:,ibeta,ieratio,inpp);
            bvsepgrid(:,:,ibeta,ieratio,inpp) = interp2(lon1pgrid,lat1grid,bvsetemp,pgrid.XT,pgrid.YT);
        end
    end
end

end % end runloadsection

%--------------------------------------------------
% SUPPLEMENTARY FIGURES 3 AND 4
% Plot all of the temporally-regressed beta vs. export maps
%--------------------------------------------------
plotsection=0;

if plotsection==1

redblue = flipud(cbrewer('div','RdYlBu',100,'linear'));
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
cmax = 1.4;

% - Global maps
for ibeta = 1:betan
    f=figure;
    set(f,'color','white','units','inches','position',[0.5 0.5 11 6]);
    colormap(redblue);
    isp =  1; % subplot idx
    for inpp = 1:nppn
        for ieratio = 1:eration
            subplot(3,3,isp);
            m_proj(mapproj,'lon',[0 360],'lat',[min(pgrid.yt) max(pgrid.yt)]);
            m_pcolor(pgrid.xt,pgrid.yt,bvsepgrid(:,:,ibeta,ieratio,inpp));
            m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
            caxis([-cmax cmax]);
            shading flat;
            [arrb,arrt]=computepointyColorbarlims(bvsepgrid(:,:,ibeta,ieratio,inpp),-cmax,cmax);
            pointyColorbar(arrb,arrt);
            m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
            title([betanamearray{ibeta} 'beta vs.' 10 NPPnamearray{inpp} ' NPP+' erationamearray{ieratio} ' e-ratio']);
            isp = isp+1;
        end
    end
end
 
end

%--------------------------------------------------
% Weight separate CbPM, VGPM, and VGPM-Eppley maps using
% Tom's regional weights from the PNAS transfer efficiency paper 
%--------------------------------------------------

% - Create weights matrix (5-D! lat, lon, betaalg, eratioalg, NPPalg)
%regnames =  {'Antarc','SubAnt','ST.Atl','ST.Pac','Trop.Atl','Trop.Pac','N.Atl','N.Pac','Indian'};
% Beta map weights defined at beginning
% Export map weights from Tom's PNAS transfer efficiency supplementary table 2
nonwtedwt=1/(eration*nppn); %1/9 
regwtstable = nan(length(regnames),eration,nppn);
%regwtstable(:,1,1) = [0.1139 0.3207 0.2308 0.0504 0.0656 0.0656 0.0478 0      nonwtedwt];
%regwtstable(:,2,1) = [0.1508 0.2328 0.1677 0.0300 0.0729 0.0729 0.0697 0.0026 nonwtedwt];
%regwtstable(:,3,1) = [0.0927 0.0454 0.0975 0.0208 0.0445 0.0445 0.1169 0.1855 nonwtedwt];
%regwtstable(:,1,2) = [0.1507 0.0420 0.1419 0.0663 0.1213 0.1213 0.1184 0.1197 nonwtedwt];
%regwtstable(:,2,2) = [0.1349 0.0212 0.0993 0.0435 0.1516 0.1516 0.1294 0.2379 nonwtedwt];
%regwtstable(:,3,2) = [0.0622 0.0036 0.0636 0.0292 0.1080 0.1080 0.1308 0.1211 nonwtedwt];
%regwtstable(:,1,3) = [0.0478 0.2014 0.0900 0.2688 0.1667 0.1667 0.1263 0.0107 nonwtedwt];
%regwtstable(:,2,3) = [0.1215 0.1141 0.0640 0.2695 0.1047 0.1047 0.1322 0.0978 nonwtedwt];
%regwtstable(:,3,3) = [0.1255 0.0188 0.0451 0.2216 0.1648 0.1648 0.1286 0.2247 nonwtedwt];
regwtstable(:,1,1) = [0.1139 0.3207 0.2308 0.0504 0.0656 0.0656 0.0478 0      nonwtedwt nonwtedwt nonwtedwt];
regwtstable(:,2,1) = [0.1508 0.2328 0.1677 0.0300 0.0729 0.0729 0.0697 0.0026 nonwtedwt nonwtedwt nonwtedwt];
regwtstable(:,3,1) = [0.0927 0.0454 0.0975 0.0208 0.0445 0.0445 0.1169 0.1855 nonwtedwt nonwtedwt nonwtedwt];
regwtstable(:,1,2) = [0.1507 0.0420 0.1419 0.0663 0.1213 0.1213 0.1184 0.1197 nonwtedwt nonwtedwt nonwtedwt];
regwtstable(:,2,2) = [0.1349 0.0212 0.0993 0.0435 0.1516 0.1516 0.1294 0.2379 nonwtedwt nonwtedwt nonwtedwt];
regwtstable(:,3,2) = [0.0622 0.0036 0.0636 0.0292 0.1080 0.1080 0.1308 0.1211 nonwtedwt nonwtedwt nonwtedwt];
regwtstable(:,1,3) = [0.0478 0.2014 0.0900 0.2688 0.1667 0.1667 0.1263 0.0107 nonwtedwt nonwtedwt nonwtedwt];
regwtstable(:,2,3) = [0.1215 0.1141 0.0640 0.2695 0.1047 0.1047 0.1322 0.0978 nonwtedwt nonwtedwt nonwtedwt];
regwtstable(:,3,3) = [0.1255 0.0188 0.0451 0.2216 0.1648 0.1648 0.1286 0.2247 nonwtedwt nonwtedwt nonwtedwt];

% - Check to make sure weight maps were constructed properly
% Looks good!
regwtsmap = nan(length(pgrid.yt),length(pgrid.xt),eration,nppn);
for inpp = 1:nppn
    for ieratio = 1:eration
        regwtsmaptemp = nan(size(regmappgrid));
        for ireg = 1:length(regnames)
            regwtsmaptemp(regmappgrid==ireg) = regwtstable(ireg,ieratio,inpp); 
        end
        regwtsmap(:,:,ieratio,inpp) = regwtsmaptemp;
    end
end
% f=figure;
% set(f,'color','white','units','inches','position',[0.5 0.5 10 7]);
% isp =  1; % subplot idx
% for inpp = 1:nppn
%     for ieratio = 1:eration
%         subplot(3,3,isp);
%         pcolor(pgrid.xt,pgrid.yt,regwtsmap(:,:,ieratio,inpp));
%         shading flat;colorbar;
%         caxis([0 max(max(max(max(regwtstable))))]);
%         title(['All beta vs.' 10 NPPnamearray{inpp} ' npp+' erationamearray{ieratio} ' e-ratio reg wts']);
%         isp = isp+1;
%     end
% end

%--------------------------------------------------
% CALCULATE OBS WEIGHTED MEAN BETA VS. EXPORT GLOBAL MAP
% + OBS WEIGHTED REGIONAL MEANS USING 
% METHOD 1: TEMPORALLY REGRESS BY GRID PT AND
% THEN SPATIALLY AVERAGE REGIONAL BETA VS. EXPORT SLOPES 
%--------------------------------------------------
regwtedmeanbvsemap = nan(length(pgrid.yt),length(pgrid.xt),betan);
regwtedmeanbvsez = nan(length(pgrid.yt),betan);
regwtedmean1stdbvsemap = nan(length(pgrid.yt),length(pgrid.xt),betan);
regwtedmeanovstdbvsemap = nan(length(pgrid.yt),length(pgrid.xt),betan);
nummodagmap = nan(length(pgrid.yt),length(pgrid.xt),betan);
regwtedmeanbvsermat = nan(length(regnames),betan,eration,nppn);
regwtedmeanbvser = nan(length(regnames),betan);
regwted1stdbvser = nan(length(regnames),betan);
bvseregblocksmap = nan(length(pgrid.yt),length(pgrid.xt),betan);

for ibeta = 1:betan

    bvsepgridnow = squeeze(bvsepgrid(:,:,ibeta,:,:)); 

    % - Compute weighted mean beta vs. export slope map
    regwtedmeanbvsemap(:,:,ibeta) = sum(sum(regwtsmap.*bvsepgridnow,4),3);

    % - Compute zonal mean beta vs. export slope 
    regwtedmeanbvsez(:,ibeta) = nansum(regwtedmeanbvsemap(:,:,ibeta).*Areatocnonly,2)./Areatzonalsum; 
    
    % - Compute weighted stdev above/below mean beta vs. export slope map
    regwted1stdbvsemap(:,:,ibeta) = zeros(length(pgrid.yt),length(pgrid.xt));
    for inpp = 1:nppn
        for ieratio = 1:eration
            regwted1stdbvsemap(:,:,ibeta) = regwted1stdbvsemap(:,:,ibeta) + regwtsmap(:,:,ieratio,inpp).*(bvsepgridnow(:,:,ieratio,inpp)-regwtedmeanbvsemap(:,:,ibeta)).^2;
        end
    end
    regwted1stdbvsemap(:,:,ibeta) = regwted1stdbvsemap(:,:,ibeta).^(1/2);
    
    % - Compute zonal mean stdev above/below mean beta vs. export slope map
    regwted1stdbvsez(:,ibeta) = nansum(regwted1stdbvsemap(:,:,ibeta).*Areatocnonly,2)./Areatzonalsum;
    
    % - Choose first or second plot type for significance hatching methods
    hatchtype=1;
    
    if hatchtype==1
        % - FIRST HATCH TYPE:
        % Map of weighted mean beta vs. export slope map w/ hatch
        % where abs(mean)/stdev > 1
        regwtedmeanovstdbvsemap(:,:,ibeta) = abs(regwtedmeanbvsemap(:,:,ibeta))./regwted1stdbvsemap(:,:,ibeta);
    
    elseif hatchtype==2
        % - SECOND HATCH TYPE: 
        % Map of weighted mean beta vs. export slope map w/ hatch
        % where greater than some % of models agree on the sign of the slope
        numag = 6;
        nummodpos = nan(length(pgrid.yt),length(pgrid.xt));
        nummodneg = nan(length(pgrid.yt),length(pgrid.xt));
        for ilon = 1:length(pgrid.xt)
            for ilat = 1:length(pgrid.yt)
                    bvsetemp = bvsepgridnow(ilat,ilon,:,:);
                    nummodpos(ilat,ilon) = length(find(bvsetemp>0));
                    nummodneg(ilat,ilon) = length(find(bvsetemp<0));
            end
        end
        nummodpos(isnan(bvsepgridnow(:,:,1,1)))=nan;
        nummodneg(isnan(bvsepgridnow(:,:,1,1)))=nan;
        nummodag = zeros(length(pgrid.yt),length(pgrid.xt));
        nummodag(nummodneg>6|nummodpos>6) = 1;
        nummodag(isnan(bvsepgridnow(:,:,1,1)))=nan;
        nummodagmap(:,:,ibeta) = nummodag;
    end
    
    % - Compute regional mean temporally regressed beta vs. export slope for each NPP/e-ratio/beta combo 
    for ieratio = 1:eration
        for inpp = 1:nppn
            bvsetemp = bvsepgridnow(:,:,ieratio,inpp);
            for iregion = 1:length(regnames)
                regwtedmeanbvsermat(iregion,ibeta,ieratio,inpp) = nansum(bvsetemp(regmappgrid==iregion).*pgrid.Areat(regmappgrid==iregion))/nansum(pgrid.Areat(regmappgrid==iregion));
            end
        end
    end
    regwtedmeanbvser(:,ibeta) = sum(sum(squeeze(regwtedmeanbvsermat(:,ibeta,:,:)).*regwtstable,2),3);

    % - Compute alternative regionally export weighted beta vs. export slopes    
    % CONFIRMED THAT THIS IS THE SAME AS REGWTEDMEANBVSER. YAY!!!
    for iregion = 1:length(regnames)
        regwtedmeanbvsemapnow = regwtedmeanbvsemap(:,:,ibeta);
        regwtedmeanbvser1(iregion,ibeta) = nansum(regwtedmeanbvsemapnow(regmappgrid==iregion).*pgrid.Areat(regmappgrid==iregion))./nansum(pgrid.Areat(regmappgrid==iregion));
    end
    
    % - Compute regional mean temporally regressed beta vs. export slope stdevs for each NPP/e-ratio/beta combo 
    for iregion = 1:length(regnames)
        regwted1stdbvser(iregion,ibeta) = sum(sum( ...
            squeeze(regwtstable(iregion,:,:)).* ...
            (squeeze(regwtedmeanbvsermat(iregion,ibeta,:,:))-regwtedmeanbvser(iregion,ibeta)).^2));
    end
    regwted1stdbvser(:,ibeta) = regwted1stdbvser(:,ibeta).^(1/2);

    % - Compute alternative regionally export weighted beta vs. export slope stdevs 
    % THESE NUMBERS ARE ALL LARGER THAN REGWTED1STDBVSER. I THINK THAT REGWTED1STDBVSER IS THE CORRECT WAY TO CALCULATE
    % THE REGIONAL MEAN STDEVS B/C WE WANT TO KNOW THE ERROR OF THE REG MEAN...I THINK...
    for iregion = 1:length(regnames)
        regwted1stdbvsemapnow = regwted1stdbvsemap(:,:,ibeta);
        regwted1stdbvser1(iregion,ibeta) = nansum(regwted1stdbvsemapnow(regmappgrid==iregion).*pgrid.Areat(regmappgrid==iregion))./nansum(pgrid.Areat(regmappgrid==iregion));
    end

    % - Create final regional block beta vs. export slope maps 
    bvseregblocksmaptemp = nan(size(regmappgrid));
    for iregion = 1:length(regnames)
        bvseregblocksmaptemp(regmappgrid==iregion) = regwtedmeanbvser(iregion,ibeta);  
    end
    bvseregblocksmap(:,:,ibeta) = bvseregblocksmaptemp;

end

%--------------------------------------------------
% CREATE PLOTS OF AN OLD VERSION OF FIGURE 4 AND SUPPLEMENTARY FIGURE 5
%--------------------------------------------------
plotsection = 0;

if plotsection==1

redblue = flipud(cbrewer('div','RdYlBu',100,'linear'));
pgrid.xtshift = [pgrid.xt(91:180)-360 pgrid.xt(1:90)];
linewidth=2; % for regional contour lines
linecolor=[0 0 0]; % for regional contour lines
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];

% - Plot global maps of b vs. e with chosen significance hatching
% w/ Kost and Guidi together in a 3 by 2 subplots fig
f=figure;
set(f,'color','white','units','inches','position',[0.5 0.5 11 10]);
for ibeta = 1:betan

    % - Reg wted mean beta vs. norm export slope map
    cmax = 1.06;
    % Kost min, max = -1.0541, 0.2229
    % Guidi min, max = -0.9772, 0.9740
    subplot(3,betan,ibeta);
    m_proj(mapproj,'lon',[0 360],'lat',[min(pgrid.yt) max(pgrid.yt)]);
    %m_pcolor(pgrid.xt,pgrid.yt,regwtedmeanbvsemap(:,:,ibeta));
    m_contourf(pgrid.xt,pgrid.yt,regwtedmeanbvsemap(:,:,ibeta),30);
    m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
    caxis([-cmax cmax]);
    colormap(redblue);freezeColors;
    shading flat;
    cbfreeze(colorbar);
    m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
    hold on;
    if hatchtype==1
        [c,h] = m_contourf(pgrid.xt,pgrid.yt,regwtedmeanovstdbvsemap(:,:,ibeta),[1 1],...
            'linewidth',1,'linecolor',[0.4 0.4 0.4],'linestyle','-');
    elseif hatchtype==2
        [c,h] = m_contourf(pgrid.xt,pgrid.yt,nummodagmap(:,:,ibeta),[1 1],...
            'linewidth',1,'linecolor',[0.4 0.4 0.4],'linestyle','-');
    end
    hp = findobj(h,'type','patch');
    %hatch(hp,[45 9 1],'k-'); %note: hatch(hp,[angle step width],'r--') 
    hatch(hp,45,[0.1 0.1 0.1],'-',8,1);
    if hatchtype==1
        title(['Reg wted mean ' betanamearray{ibeta} ' \beta vs. norm export slope map']);
    elseif hatchtype==2
        title(['Reg wted mean ' betanamearray{ibeta} ' \beta vs. norm export slope map; # mod agr=' num2str(numag)]);
    end

    % - Reg mean beta vs. norm export slopes bar chart 
    % end-2 is to exclude Med and Arc regions from the final plot
    % Reorder Indian Ocean
    subplot(3,betan,ibeta+betan);
    bh=bar(1:length(regnames(1:end-2)),[regwtedmeanbvser(1:4,ibeta); regwtedmeanbvser(9,ibeta); regwtedmeanbvser(5:8,ibeta)]);hold on;
    bc=get(bh,'Children');
    set(bc,'CData',[regwtedmeanbvser(1:4,ibeta); regwtedmeanbvser(9,ibeta); regwtedmeanbvser(5:8,ibeta)]);
    caxis([-max(max(abs(regwtedmeanbvser))) max(max(abs(regwtedmeanbvser)))]);
    colormap(redblue);freezeColors;
    cbfreeze(colorbar);
    errorbar(1:length(regnames(1:end-2)),[regwtedmeanbvser(1:4,ibeta); regwtedmeanbvser(9,ibeta); regwtedmeanbvser(5:8,ibeta)],...
        [regwted1stdbvser(1:4,ibeta); regwted1stdbvser(9,ibeta); regwted1stdbvser(5:8,ibeta)],'k.');
    xlim([0.5 length(regnames(1:end-2))+0.5]);
    set(gca, 'XTickLabel', regnamesabbrevro(1:end-2));
    %set(gca, 'XTick', 1:length(regnames(1:end-2)), 'XTickLabel', regnamesabbrevro(1:end-2));
    % % OR: 
    %NumTicks = length(regnames(1:end-2));
    %L = get(gca,'XLim');
    %set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    %set(gca,'XTickLabel',[' ',regnamesabbrevro(1:end-2),' ']);
    ylim([-0.6 0.2]);
    ylabel('\beta vs. export slope');
    title('Reg mean \beta vs. norm export slope');
    %grid on;

    % - Final regional blocks global map
    subplot(3,betan,ibeta+2*betan);
    m_proj(mapproj,'lon',[0 360],'lat',[min(pgrid.yt) max(pgrid.yt)]);
    m_pcolor(pgrid.xt,pgrid.yt,bvseregblocksmap(:,:,ibeta));
    m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
    caxis([-max(max(abs(regwtedmeanbvser))) max(max(abs(regwtedmeanbvser)))]);
    colormap(redblue);freezeColors;
    shading flat;cbfreeze(colorbar);
    m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
    hold on;
    nreg = max(regmappgrid(:))
    for i = 1:nreg
        RR = double(regmappgrid==i);
        RR(M3d(:,:,1)==0) = NaN;
        m_contour(pgrid.xt,pgrid.yt,RR,[.5 .5],'k','linewidth',linewidth,'linecolor',linecolor);hold on;
    end
    title('Final global \beta vs. norm export slope map');
end

end % end plotsection

%--------------------------------------------------
% FIGURE 4 AND SUPPLEMENTARY FIGURE 5
% Plot beta vs. export regionally-weighted global map, all beta vs. export
% regional means, regionally-weighted beta vs. export regional means,
% final global map of beta vs. export regionally-weighted regional means
%--------------------------------------------------
plotsection = 1;

if plotsection==1

redblue = flipud(cbrewer('div','RdYlBu',100,'linear'));
pgrid.xtshift = [pgrid.xt(91:180)-360 pgrid.xt(1:90)];
linewidth=2; % for regional contour lines
linecolor=[0 0 0]; % for regional contour lines
mapproj = 'gall-peters';
landcolor = [0.6 0.6 0.6];
labelfontsize = 14;

% - Plot global maps of b vs. e with chosen significance hatching
% and w/ all 9 regional mean b vs. e values,
% w/ Kost and Guidi in diff figs with 2 by 2 subplots
qual9 = cbrewer('qual','Set1',9);
for ibeta = 1:betan
    f=figure;
    set(f,'color','white','units','inches','position',[0.5 0.5 13 7]);

    % - Reg wted mean beta vs. norm export slope map
    cmax = 1.06;
    % Kost min, max = -1.0541, 0.2229
    % Guidi min, max = -0.9772, 0.9740
    subplot(2,2,1);
    m_proj(mapproj,'lon',[0 360],'lat',[min(pgrid.yt) max(pgrid.yt)]);
    m_contourf(pgrid.xt,pgrid.yt,regwtedmeanbvsemap(:,:,ibeta),30);
    %m_pcolor(pgrid.xt,pgrid.yt,regwtedmeanbvsemap(:,:,ibeta));
    m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
    caxis([-cmax cmax]);
    colormap(redblue);freezeColors;
    shading flat;
    h=cbfreeze(colorbar);
    set(h,'fontsize',labelfontsize);
    m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
    hold on;
    if hatchtype==1
        [c,h] = m_contourf(pgrid.xt,pgrid.yt,regwtedmeanovstdbvsemap(:,:,ibeta),[1 1],...
            'linewidth',1,'linecolor',[0.4 0.4 0.4],'linestyle','-');
    elseif hatchtype==2
        [c,h] = m_contourf(pgrid.xt,pgrid.yt,nummodagmap(:,:,ibeta),[1 1],...
            'linewidth',1,'linecolor',[0.4 0.4 0.4],'linestyle','-');
    end
    hp = findobj(h,'type','patch');
    %hatch(hp,[45 9 1],'k-'); %note: hatch(hp,[angle step width],'r--') 
    hatch(hp,45,[0.1 0.1 0.1],'-',8,1);
    if hatchtype==1
        title(['Reg wted mean ' betanamearray{ibeta} ' \beta vs. norm export slope map']);
    elseif hatchtype==2
        title(['Reg wted mean ' betanamearray{ibeta} ' \beta vs. norm export slope map; # mod agr=' num2str(numag)]);
    end

    % - Reg mean beta vs. norm export slopes for all 9 combos of NPP + e-ratio
    subplot(2,2,2);
    regwtedmeanbvsermattemp1 = squeeze(regwtedmeanbvsermat(:,ibeta,:,:));
    % end-2 is to exclude Med and Arc region from the final plot
    regwtedmeanbvsermattemp = [regwtedmeanbvsermattemp1(1:end-2,:,1) regwtedmeanbvsermattemp1(1:end-2,:,2) regwtedmeanbvsermattemp1(1:end-2,:,3)];
    % Reorder Indian Ocean
    regwtedmeanbvsermattemp = [regwtedmeanbvsermattemp(1:4,:); regwtedmeanbvsermattemp(9,:); regwtedmeanbvsermattemp(5:8,:)]; 
    bh=bar(1:length(regnames(1:end-2)),regwtedmeanbvsermattemp,1);
    set(bh,'edgecolor',[0.35 0.35 0.35]);
    colormap(qual9);freezeColors;
    exportlabels = {[NPPnamearray{1} '+' eratioshortnamearray{1}],[NPPnamearray{2} '+' eratioshortnamearray{1}],[NPPnamearray{3} '+' eratioshortnamearray{1}],...
        [NPPnamearray{1} '+' eratioshortnamearray{2}],[NPPnamearray{2} '+' eratioshortnamearray{2}],[NPPnamearray{3} '+' eratioshortnamearray{2}],...
        [NPPnamearray{1} '+' eratioshortnamearray{3}],[NPPnamearray{2} '+' eratioshortnamearray{3}],[NPPnamearray{3} '+' eratioshortnamearray{3}]};
    caxis([0 9]);
    h=cbfreeze(colorbar('YTick',0.5:length(exportlabels)+0.5,'YTickLabels',exportlabels));
    %cbfreeze(colorbar('YTickLabels',exportlabels));
    set(h,'fontsize',12,'fontweight','bold');
    xlim([0.5 length(regnames(1:end-2))+0.5]);
    set(gca,'XTickLabel',regnamesabbrevro(1:end-2),'fontsize',labelfontsize-2,'fontweight','bold');
    ylim([-0.65 0.2]);
    ylabel('\beta vs. export slope');
    title('Reg mean \beta vs. norm export slope');

     % - Reg mean beta vs. norm export slopes bar chart 
     % end-2 is to exclude Med and Arc regions from the final plot
     % Reorder Indian Ocean
     subplot(2,2,3);
     bh=bar(1:length(regnames(1:end-2)),[regwtedmeanbvser(1:4,ibeta); regwtedmeanbvser(9,ibeta); regwtedmeanbvser(5:8,ibeta)]);hold on;
     bc=get(bh,'Children');
     set(bc,'CData',[regwtedmeanbvser(1:4,ibeta); regwtedmeanbvser(9,ibeta); regwtedmeanbvser(5:8,ibeta)]);
     caxis([-max(max(abs(regwtedmeanbvser))) max(max(abs(regwtedmeanbvser)))]);
     colormap(redblue);freezeColors;
     h=cbfreeze(colorbar);
     set(h,'fontsize',labelfontsize);
     errorbar(1:length(regnames(1:end-2)),[regwtedmeanbvser(1:4,ibeta); regwtedmeanbvser(9,ibeta); regwtedmeanbvser(5:8,ibeta)],...
         [regwted1stdbvser(1:4,ibeta); regwted1stdbvser(9,ibeta); regwted1stdbvser(5:8,ibeta)],'k.');
     xlim([0.5 length(regnames(1:end-2))+0.5]);
     set(gca,'XTickLabel', regnamesabbrevro(1:end-2),'fontsize',labelfontsize-2,'fontweight','bold');
     %set(gca, 'XTick', 1:length(regnames(1:end-2)), 'XTickLabel', regnamesabbrevro(1:end-2));
     % % OR: 
     %NumTicks = length(regnames(1:end-2));
     %L = get(gca,'XLim');
     %set(gca,'XTick',linspace(L(1),L(2),NumTicks))
     %set(gca,'XTickLabel',[' ',regnamesabbrevro(1:end-2),' ']);
     ylim([-0.6 0.2]);
     ylabel('\beta vs. export slope');
     title('Reg mean \beta vs. norm export slope');
     %grid on;
 
     % - Final regional blocks global map
     subplot(2,2,4);
     m_proj(mapproj,'lon',[0 360],'lat',[min(pgrid.yt) max(pgrid.yt)]);
     m_pcolor(pgrid.xt,pgrid.yt,bvseregblocksmap(:,:,ibeta));
     m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
     caxis([-max(max(abs(regwtedmeanbvser))) max(max(abs(regwtedmeanbvser)))]);
     colormap(redblue);freezeColors;
     shading flat;h=cbfreeze(colorbar);
     set(h,'fontsize',labelfontsize);
     m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
     hold on;
     nreg = max(regmappgrid(:))
     for i = 1:nreg
         RR = double(regmappgrid==i);
         RR(M3d(:,:,1)==0) = NaN;
         m_contour(pgrid.xt,pgrid.yt,RR,[.5 .5],'k','linewidth',linewidth,'linecolor',linecolor);hold on;
     end
     title('Final global \beta vs. norm export slope map');
end

end % end plotsection

%--------------------------------------------------
% CREATE SUPPLEMENTARY PLOTS!!!!!
%--------------------------------------------------

plotsection=0;

if plotsection==1

f=figure;
set(f,'color','white','units','inches','position',[0.5 0.5 11 10]);

for ibeta = 1:betan

    % - Reg wted beta vs. norm export stdev map
    subplot(3,betan,ibeta);
    m_proj(mapproj,'lon',[0 360],'lat',[min(pgrid.yt) max(pgrid.yt)]);
    m_pcolor(pgrid.xt,pgrid.yt,regwted1stdbvsemap(:,:,ibeta));
    m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
    caxis([-max(max(abs(regwted1stdbvsemap(:,:,ibeta)))) max(max(abs(regwted1stdbvsemap(:,:,ibeta))))]);
    colormap(redblue);freezeColors;
    shading flat;cbfreeze(colorbar);
    m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
    hold on;
    nreg = max(regmappgrid(:))
    for i = 1:nreg
        RR = double(regmappgrid==i);
        RR(M3d(:,:,1)==0) = NaN;
        m_contour(pgrid.xt,pgrid.yt,RR,[.5 .5],'k','linewidth',linewidth,'linecolor',linecolor);hold on;
    end
    title(['Reg wted mean' betanamearray{ibeta} '\beta vs. norm export slope 1-stdev map']);

    % - Reg mean beta vs. norm export slopes for all 9 combos of NPP + e-ratio
    subplot(3,betan,ibeta+betan);
    regwtedmeanbvsermattemp1 = squeeze(regwtedmeanbvsermat(:,ibeta,:,:));
    % end-1 is to exclude Med region from the final plot
    regwtedmeanbvsermattemp = [regwtedmeanbvsermattemp1(1:end-1,:,1) regwtedmeanbvsermattemp1(1:end-1,:,2) regwtedmeanbvsermattemp1(1:end-1,:,3)];
    bar(1:length(regnames(1:end-1)),regwtedmeanbvsermattemp);
    colormap(jet);
    xlim([0.5 length(regnames(1:end-1))+0.5]);
    set(gca, 'XTickLabel', regnamesabbrev(1:end-1));
    ylim([-0.65 0.2]);
    ylabel('\beta vs. export slope');
    title('Reg mean \beta vs. norm export slope');

end

end % end plotsection

% %--------------------------------------------------
% % OLD FIRST ATTEMPT AT PLOTTING
% %--------------------------------------------------
% 
% % - Calculate max and min lats for each region
% minreglats = nan(1,length(regnames));
% maxreglats = nan(1,length(regnames));
% for ireg = 1:length(regnames)
%     [latidx,lonidx]=find(regmappgrid==ireg);
%     minreglats(ireg) = pgrid.yt(min(latidx));
%     maxreglats(ireg) = pgrid.yt(max(latidx));
% end
% 
% % - Create plot!!!
% f=figure;
% set(f,'color','white','units','inches','position',[0.5 0.5 5 7]);
% colors = ...
%     [    0         0    1.0000;
%     1.0000         0         0;
%          0    1.0000         0;
%     1.0000    0.1034    0.7241;
%     1.0000    0.8276         0;
%          0    0.3448         0;
%     0.5172    0.5172    1.0000;
%     0.6207    0.3103    0.2759;
%          0    1.0000    0.7586;
%          0         0    0.1724;]
% 
% % - Plot global maps of b vs. e with chosen significance hatching
% subplot(211);
% pcolor(pgrid.xt,pgrid.yt,regwtedmeanbvsemap);hold on;
% colormap(bipolar);
% caxis([-1 1]);
% shading flat;colorbar;
% hold on;
% if hatchtype==1
%     [c,h] = contourf(pgrid.xt,pgrid.yt,regwtedmeanovstdbvsemap,[1 1],'linewidth',2);
% elseif hatchtype==2
%     [c,h] = contourf(pgrid.xt,pgrid.yt,nummodag,[1 1],'linewidth',2);
% end
% hp = findobj(h,'type','patch');
% hatch(hp,[45 7 1],'k-'); %note: hatch(hp,[angle step width],'r--') 
% if hatchtype==1
%     title(['Regionally wted mean ' betaname ' \beta vs. normalized export slope map']);
% elseif hatchtype==2
%     title(['Regionally wted mean ' betaname ' \beta vs. normalized export slope map; # mod agr=' num2str(numag)]);
% end
% 
% % - Plot the final regional mean temp regressed beta vs. export slopes w/ error bars
% subplot(212);
% plot(pgrid.yt,regwtedmeanbvsez,'k-','linewidth',2); hold on;
% %plot(pgrid.yt,regwtedmeanbvsez+regwted1stdbvsez,'k--','linewidth',1);
% %plot(pgrid.yt,regwtedmeanbvsez-regwted1stdbvsez,'k--','linewidth',1);
% for ireg = 1:length(regnames)
%     p(ireg)=plot([minreglats(ireg) maxreglats(ireg)],[regwtedmeanbvser(ireg) regwtedmeanbvser(ireg)],'color',colors(ireg,:),'linewidth',2);
% end
% xlim([-80 80]);
% xlabel('Lat');
% ylabel('\beta vs. export slope');
% title('Zonal or regional mean \beta vs. normalized export slope');
% grid on;
% 
% disp('pause')
% pause
% 
% errorbar(regwtedmeanbvse_tblocks1,regwted1stdbvse_tblocks1,'k');hold on;
% plot([1 11],[0 0],'k--');
% NumTicks = length(regnames)+2;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% set(gca,'XTickLabel',[' ',regnames,' ']);
% ylabel('\beta vs. export slope');
% suplabel(['monthly temp regress by grid pt then spat avg - ' betaname ' beta'],'t');
% 
% % - Plot all of the regional mean temp regressed beta vs. export slopes for each NPP/e-ratio/beta combo (18 total) w/ the final slopes bolded
% f=figure; set(f,'color','white','units','inches','position',[0.5 0.5 10 10]);
% subplot(121);
% plotcols = distinguishable_colors(numel(regwtedmeanbvse_tblocks1mat));
% icol=1;
% for ieratio = 1:eration
%     for inpp = 1:nppn
%         plot(regwtedmeanbvse_tblocks1mat(:,ieratio,inpp),'color',plotcols(icol,:));hold on;
%         icol = icol+1;
%     end
% end
% plot(regwtedmeanbvse_tblocks1,'k','linewidth',4);
% plot([1 11],[0 0],'k--');
% NumTicks = length(regnames);
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% set(gca,'XTickLabel',regnames);
% ylabel('\beta vs. export slope');
% subplot(122);
