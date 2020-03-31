clear all;
close all;
addpath(genpath('/graid1/shirlleu/matlabroutines'));
addpath(genpath('/jetwork2/cram/matlabtools/m_map'));

nummos = 160; % length of export and beta time series

%--------------------------------------------------
% Load and plot all export maps
% (Everything on satellite grid since this is how everything was computed)
%--------------------------------------------------
% - Initialize export filenames
opts.NPPalgarray = {'VGPM','VGPME','CbPM'};
% choose: VGPM, CbPM, VGPME, physatchl1
opts.eratioezdepthstr = 'spatvareratioandezdepth';
% choose: spatvareratioandezdepth, spatvareratioconstezdepth,
% consteratioandezdepth, consteratiospatvarezdepth
opts.fnamebaseexpC = ['/graid1/shirlleu/NPPslopeproj/matfiles/globalexportmaps_' opts.eratioezdepthstr '/'];
opts.eratioalgarray = {'Laws2000','D2005PP','Laws2011D'};
% choose: 'EP1979','Laws2000','D2005PP','D2005chl','Lutz2007','H2011','Laws2011N','Laws2011D'
eration = length(opts.eratioalgarray);
nppn = length(opts.NPPalgarray);

% - Initialize matrix to hold all annual export maps
load([opts.fnamebaseexpC 'expCann_D2005PPfromphysatchl1.mat'],'lat1','lon1','PP1ann');
lon1 = [lon1(181:end)' lon1(1:180)'+360];
expCalgs = nan(length(lat1),length(lon1),nummos,length(opts.NPPalgarray),length(opts.eratioalgarray)); % num of months, num of diff NPP maps, num of diff eratio maps 
PPalgs = nan(length(lat1),length(lon1),nummos,length(opts.NPPalgarray),length(opts.eratioalgarray)); % num of months, num of diff NPP maps, num of diff eratio maps

% - Load export maps
for iNPPalg = 1:length(opts.NPPalgarray)
    for ieratioalg = 1:length(opts.eratioalgarray)
        load([opts.fnamebaseexpC 'expCmo_' opts.eratioalgarray{ieratioalg} 'from' opts.NPPalgarray{iNPPalg} '.mat'],'expCmo1','PPmo1','exptime');
        expCalgs(:,:,:,iNPPalg,ieratioalg) = [expCmo1(:,181:end,:) expCmo1(:,1:180,:)];
        PPalgs(:,:,:,iNPPalg,ieratioalg) = [PPmo1(:,181:end,:) PPmo1(:,1:180,:)];
    end
end
time = exptime;

% SUPPLEMENTARY FIGURE 2
% - Plot all export maps
plotsection=1;

if plotsection==1
  cmap=cbrewer('seq','YlGnBu',100,'linear');
  mapproj = 'gall-peters';
  landcolor = [0.6 0.6 0.6];
  f=figure;
  set(f,'color','white','units','inches','position',[0.5 0.5 15 7]);
  isubplot=1;
  for iNPPalg = 1:length(opts.NPPalgarray)
      for ieratioalg = 1:length(opts.eratioalgarray)
          expCtemp = nanmean(expCalgs(:,:,:,iNPPalg,ieratioalg),3);
          m_proj(mapproj,'lon',[0 360],'lat',[min(lat1) max(lat1)])
          subplot(3,3,isubplot);
          %m_contourf(lon1,lat1,expCtemp,20);
          m_pcolor(lon1,lat1,expCtemp);
          m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
          caxis([0 8]);
          colormap(flipud(cmap));
          shading flat;
          [arrb,arrt]=computepointyColorbarlims(expCtemp,0,8);
          pointyColorbar(arrb,arrt);
          m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
          title([opts.NPPalgarray{iNPPalg} ' NPP+' opts.eratioalgarray{ieratioalg} ' e-ratio']);
          isubplot=isubplot+1;
      end
  end
  %suplabel(['Export [molC/m^2/yr] - ' opts.eratioezdepthstr],'t');
end % end plotsection

%--------------------------------------------------
% Weight separate CbPM, VGPM, and VGPM-Eppley maps using
% Tom's regional weights from the PNAS transfer efficiency paper 
% (Put everything on satellite grid)
%--------------------------------------------------

% - Load region map
franken = load('/graid1/shirlleu/NPPslopeproj/matfiles/franken_export.mat');
fgrid = franken.grid;
regmap = franken.R2d;
regmap = interp2([fgrid.xt(end)-360 fgrid.xt],fgrid.yt,[regmap(:,end) regmap],lon1,lat1,'nearest');
regnames = franken.regnames;
regnames = [regnames,'Indian','Arctic','MedSeaPlus'];
regnamesabbrev = {'AAZ','SAZ','STA','STP','ETA','ETP','NA','NP','IND','ARC','MED'};
[lon1grid,lat1grid] = meshgrid(lon1,lat1);
regmap(regmap==0 & lat1grid<30) = 9; % Indian Ocean region
regmap(regmap==0 & lat1grid>60) = 10; % Arctic Ocean region
regmap(regmap==0) = 11; % Mediterannean Sea + other seas

% - Create weights matrix (5-D! lat, lon, betaalg, eratioalg, NPPalg)
%regnames =  {'Antarc','SubAnt','ST.Atl','ST.Pac','Trop.Atl','Trop.Pac','N.Atl','N.Pac','Indian'};
% Beta map weights defined at beginning
% Export map weights from Tom's PNAS transfer efficiency supplementary table 2
nonwtedwt=1/(eration*nppn); %1/9 
regwtstable = nan(length(regnames),eration,nppn);
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
regwtsmap = nan(length(lat1),length(lon1),eration,nppn);
for inpp = 1:nppn
    for ieratio = 1:eration
        regwtsmaptemp = nan(size(regmap));
        for ireg = 1:length(regnames)
            regwtsmaptemp(regmap==ireg) = regwtstable(ireg,ieratio,inpp);
        end
        regwtsmap(:,:,ieratio,inpp) = regwtsmaptemp;
    end
end
plotsection=0;
if plotsection==1
  f=figure;
  set(f,'color','white','units','inches','position',[0.5 0.5 10 7]);
  isp =  1; % subplot idx
  for inpp = 1:nppn
      for ieratio = 1:eration
          subplot(3,3,isp);
          pcolor(lon1,lat1,regwtsmap(:,:,ieratio,inpp));
          shading flat;colorbar;
          caxis([0 max(max(max(max(regwtstable))))]);
          title(['All beta vs.' 10 opts.NPPalgarray{inpp} ' npp+' opts.eratioalgarray{ieratio} ' e-ratio reg wts']);
          isp = isp+1;
      end
  end
end % end plotsection

%--------------------------------------------------
% Create franken export monthly time series
%--------------------------------------------------
expCfrankents = zeros(length(lat1),length(lon1),nummos); 
for itime = 1:nummos 
    for iNPPalg = 1:length(opts.NPPalgarray)
        for ieratioalg = 1:length(opts.eratioalgarray)
            expCfrankents(:,:,itime) = expCfrankents(:,:,itime) + expCalgs(:,:,itime,iNPPalg,ieratioalg).*regwtsmap(:,:,ieratioalg,iNPPalg);
        end
    end
end

%--------------------------------------------------
% Load and plot all beta maps
% (Everything on satellite grid since this is how everything was computed)
%--------------------------------------------------
% - Initialize beta file/map names
opts.fnamebasebeta = '/graid1/shirlleu/NPPslopeproj/matfiles/betamaps/';
opts.betaalgarray = {'Kost','GuidiallUVP'};

% - Initialize matrix to hold all monthly mean beta maps
% for gridding purposes (use physat chl+derivatives grid, which
% is ALMOST EXACTLY (shifted by ~0.1 degrees) the same
% as the Kostadinov beta grid)
betaalgs = nan(length(lat1),length(lon1),nummos,length(opts.betaalgarray)); 
 
% - Load monthly Kost beta map
kost = load([opts.fnamebasebeta 'Xi1mol.mat'],'Ximo1');
betaalgs(:,:,:,1) = [kost.Ximo1(:,181:end,:) kost.Ximo1(:,1:180,:)];
    
% - Load monthly Guidi map
% (map has the same lon and lat as Kost beta map already)
guidi = load('/graid1/shirlleu/NPPslopeproj/matfiles/betamaps/guidialluvp_beta1mo.mat');
betaalgs(:,:,:,2) = [guidi.betamo1(:,181:end,:) guidi.betamo1(:,1:180,:)];
 
% - Plot all beta maps
plotsection=0

if plotsection==1
  cmap=cbrewer('seq','YlGnBu',100,'linear');
  mapproj = 'gall-peters';
  landcolor = [0.6 0.6 0.6];
  f=figure;
  set(f,'color','white','units','inches','position',[0.5 0.5 10 3.5]);
  isubplot=1;
  for ialg=1:size(betaalgs,4)
      betatemp = nanmean(betaalgs(:,:,:,ialg),3);
      m_proj(mapproj,'lon',[0 360],'lat',[min(lat1) max(lat1)])
      subplot(1,2,isubplot);
      m_pcolor(lon1,lat1,betatemp);
      m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
      caxis([2 6]);
      colormap(cmap);
      shading flat;
      freezeColors;cbfreeze(colorbar);
      m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
      title(opts.betaalgarray{ialg});
      isubplot=isubplot+1;
  end
  %suplabel('beta','t');
end % end plotsection

%--------------------------------------------------
% Load and plot Guidi UVP locations
% (Everything on satellite grid since this is how everything was computed)
%--------------------------------------------------
twguidipath = ['/jetwork2/tweber/Data/particles/' ...
                 'base_Weber_2_4_5_slope_7_15_new.mat']
load(twguidipath); % base_all
ns = length(base_all);
glat = nan(1,ns);
glon = nan(1,ns);
for i = 1:ns
    glat(i) = base_all(i).latitude;
    glon(i) = base_all(i).longitude;
    gdate(i,1) = base_all(i).datem;
    mPSDt10(i) = nanmean(base_all(i).slope(find(base_all(i).depth<=10)));
    mPSDt20(i) = nanmean(base_all(i).slope(find(base_all(i).depth<=20)));
    mPSDt40(i) = nanmean(base_all(i).slope(find(base_all(i).depth<=40)));
    mPSDt50(i) = nanmean(base_all(i).slope(find(base_all(i).depth<=50)));
end 
glon(glon<0)=glon(glon<0)+360;

%--------------------------------------------------
% For plotting Guidi UVP measurements on top of Guidi beta
% global map, decide which depth of depth-avged beta to use
%--------------------------------------------------
mPSDmat = [-mPSDt10' -mPSDt20' -mPSDt40' -mPSDt50'];

plotsection=0;
if plotsection==1
  cmap=cbrewer('seq','YlGnBu',100,'linear');
  mapproj = 'gall-peters';
  landcolor = [0.6 0.6 0.6];
  f=figure;
  set(f,'color','white','units','inches','position',[0.5 0.5 10 7]);
  for idepth=1:4
      subplot(2,2,idepth);
      mPSDnow = mPSDmat(:,idepth);
      m_proj(mapproj,'lon',[0 360],'lat',[min(lat1) max(lat1)])
      m_scatter(glon,glat,40,mPSDnow);
      m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',10);
      caxis([2 6]);
      colormap(cmap);
      shading flat;
      freezeColors; cbfreeze(colorbar);
      m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
  end
  
  f=figure;
  set(f,'color','white','units','inches','position',[0.5 0.5 10 7]);
  filledgh = inpaint_nans(betatemp);
  ghPSD_coresp = interp2(lon1grid,lat1grid,filledgh,glon,glat)';
  for idepth=1:4
      subplot(2,2,idepth);
      mPSDnow = mPSDmat(:,idepth);
      ikeep = ~(isnan(ghPSD_coresp) | isnan(mPSDnow));
      coeffs = polyfit(ghPSD_coresp(ikeep), mPSDnow(ikeep), 1);
      fittedX = linspace(min(ghPSD_coresp), max(ghPSD_coresp), length(mPSDnow(ikeep)));
      fittedY = polyval(coeffs, fittedX);
      [jnk1,jnk2,jnk3,jnk4,stats]=regress(mPSDnow(ikeep),[ones(size(ghPSD_coresp(ikeep))) ghPSD_coresp(ikeep)]);
      scatter(ghPSD_coresp, mPSDnow); hold on;
      plot(fittedX,fittedX,'k');
      plot(fittedX,fittedY,'r');
      ylim([0 8]);
      xlim([2 6]);
      xlabel('hirata+guidi \beta');
      ylabel('UVP \beta');
      title(['R^2=' num2str(stats(1))]);
  end

end % end plotsection

% Based on the plots above, let's use the avg beta over the top 50 m

%--------------------------------------------------
% Choose lats and lons defining boxes over which I'll avg
% to get export and beta time series/scatter plots
%--------------------------------------------------
% - SUBTROPICAL EXAMPLE
%latn = -33;
%lats = -34; 
latn = -26;
lats = -34;
%lonl = -130.5;
%lonr = -129.5;
%lonl = 229.5;
%lonr = 230.5;
lonl = 225.5;
lonr = 234.5;
latidxst=find(lat1<latn&lat1>lats);
lonidxst=find(lon1<lonr&lon1>lonl);

% - TROPICAL EXAMPLE
%latn = 0.5;
latn = 3;
%lats = -0.5;
lats = -3;
%lonl = 175.5;
lonl = 260;
%lonr = 176.5;
lonr = 270;
latidxt=find(lat1<latn&lat1>lats);
lonidxt=find(lon1<lonr&lon1>lonl);

%--------------------------------------------------
% MAIN TEXT FIGURE 3A AND 3B
% Plot franken export map next to Guidi beta map.
% Guidi beta map w/ in situ obs locations;
% Franken export map w/ tropical and subtropical
% example regions boxed.
% Both only the export map with regional contours.
%--------------------------------------------------
plotsection=1;

if plotsection==1
  cmap=cbrewer('seq','YlGnBu',100,'linear');
  mapproj = 'gall-peters';
  landcolor = [0.6 0.6 0.6];
  linewidth=2; % for regional contour lines
  linecolor=[0 0 0]; % for regional contour lines
  labelfontsize = 14;

  f=figure;
  set(f,'color','white','units','inches','position',[0.5 0.5 10 3.5],'resize','off');
  
  plotregcontours=0;
  subplot(121);
  ialg=2; % 2 = hirata+guidi
  betatemp = nanmean(betaalgs(:,:,:,ialg),3);
  m_proj(mapproj,'lon',[0 360],'lat',[min(lat1) max(lat1)])
  m_contourf(lon1,lat1,betatemp,25);
  %m_pcolor(lon1,lat1,betatemp);
  m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
  caxis([2 6]);
  colormap(cmap);
  shading flat;
  freezeColors; h=cbfreeze(colorbar);
  set(h,'fontsize',labelfontsize);
  m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
  hold on;
  colormap([1 1 1; 1 1 1; 1 1 1]);
  m_scatter(glon,glat,15,-mPSDt50','o');
  freezeColors; %cbfreeze(colorbar);
  if plotregcontours==1
    %nreg = max(regmap(:))
    nreg = 8;
    for i = 1:nreg
        RR = double(regmap==i);
        RR(isnan(regmap)) = NaN;
        m_contour(lon1,lat1,RR,[.5 .5],'k','linewidth',linewidth,'linecolor',linecolor);hold on;
    end
  end
  title(opts.betaalgarray{ialg});
  
  plotregcontours=1;
  expCfrank = interp2([fgrid.xt(end)-360 fgrid.xt],fgrid.yt,[franken.Cfrank(:,end) franken.Cfrank],lon1,lat1);
  % Next few lines are to fill in the small nan islands east of Australia
  latn = -10; lats = -30; lonl = 160; lonr = 184;
  latidxnan=find(lat1<latn&lat1>lats);
  lonidxnan=find(lon1<lonr&lon1>lonl);
  expCfrank1 = inpaint_nans(expCfrank(latidxnan,lonidxnan)); % interpolate over small nan islands
  expCfrank(latidxnan,lonidxnan)=expCfrank1; % fill in the small nan islands
  %figure; pcolor(flipud(expCfrank)); shading flat; colorbar;
  % Next few lines are to fill in the small nan islands south east of Australia
  latn = -35; lats = -48; lonl = 160; lonr = 184;
  latidxnan=find(lat1<latn&lat1>lats);
  lonidxnan=find(lon1<lonr&lon1>lonl);
  expCfrank1 = inpaint_nans(expCfrank(latidxnan,lonidxnan)); % interpolate over small nan islands
  expCfrank(latidxnan,lonidxnan)=expCfrank1; % fill in the small nan islands
  %figure; pcolor(flipud(expCfrank)); shading flat; colorbar;
  % Next few lines are to fill in the small nan islands at the NA/STA boundary
  latn = 40; lats = 35; lonl = 333; lonr = 337;
  latidxnan=find(lat1<latn&lat1>lats);
  lonidxnan=find(lon1<lonr&lon1>lonl);
  expCfrank1 = inpaint_nans(expCfrank(latidxnan,lonidxnan)); % interpolate over small nan islands
  expCfrank(latidxnan,lonidxnan)=expCfrank1; % fill in the small nan islands
  %figure; pcolor(flipud(expCfrank)); shading flat; colorbar;
  % Next few lines are to fill in the small nan islands in the Trop Pac
  latn = 23; lats = 18; lonl = 200; lonr = 210;
  latidxnan=find(lat1<latn&lat1>lats);
  lonidxnan=find(lon1<lonr&lon1>lonl);
  expCfrank1 = inpaint_nans(expCfrank(latidxnan,lonidxnan)); % interpolate over small nan islands
  expCfrank(latidxnan,lonidxnan)=expCfrank1; % fill in the small nan islands
  %figure; pcolor(flipud(expCfrank)); shading flat; colorbar;
  subplot(122);
  m_proj(mapproj,'lon',[0 360],'lat',[min(lat1) max(lat1)])
  m_contourf(lon1,lat1,log10(expCfrank),30);
  %m_pcolor(lon1,lat1,log10(expCfrank));
  m_grid('xtick',[0 90 180 270 360],'xlabeldir','middle','fontsize',labelfontsize);
  caxis([-0.5 1]);
  %caxis([0 8]);
  colormap(flipud(cmap));
  shading flat;
  [arrb,arrt]=computepointyColorbarlims(log10(expCfrank),-0.5,1);
  pointyColorbar(arrb,arrt);
  h=pointyColorbar(arrb,arrt);
  set(h,'fontsize',labelfontsize);
  m_coast('patch',landcolor,'edgecolor','k'); % must go after shading flat
  if plotregcontours==1
    hold on;
    %nreg = max(regmap(:))
    nreg = 8;
    for i = 1:nreg
        RR = double(regmap==i);
        RR(isnan(regmap)) = NaN;
        m_contour(lon1,lat1,RR,[.5 .5],'k','linewidth',linewidth,'linecolor',linecolor);hold on;
    end
  end
  hold on;
  tropex = zeros(size(regmap));
  tropex(latidxt,lonidxt) = 1;
  tropex(isnan(regmap)) = NaN;
  m_contour(lon1,lat1,tropex,[.5 .5],'linewidth',linewidth,'linecolor',[1 1 1]);hold on;
  hold on;
  subtropex = zeros(size(regmap));
  subtropex(latidxst,lonidxst) = 1;
  subtropex(isnan(regmap)) = NaN;
  m_contour(lon1,lat1,subtropex,[.5 .5],'linewidth',linewidth,'linecolor',[1 1 1]);hold on;
  title('franken export');
end % end plotsection

%--------------------------------------------------
% MAIN TEXT FIGURE 3C AND 3D
% Guidi beta versus franken export scatterplots 
% and time series
%--------------------------------------------------
plotsection=1;

if plotsection==1

f=figure;
set(f,'color','white','units','inches','position',[0.5 0.5 10 7]);
for ireg = 1:2

    if ireg==1 % tropics
        lonidx = lonidxt;
        latidx = latidxt;
        expCmin = 0;
        expCmax = 2;
        betamin = 3.25; 
        betamax = 4.25;
        betaticks = [betamin:0.25:betamax];
        iNPPalgme = 3;
        ieratioalgme = 3;
        iNPPalgle = 1;
        ieratioalgle = 1;
    elseif ireg==2 % subtropics
        lonidx = lonidxst;
        latidx = latidxst;
        expCmin = 0;
        expCmax = 2.5;
        betamin = 4; 
        betamax = 6.5;
        betaticks = [betamin:0.5:betamax];
        iNPPalgme = 2;
        ieratioalgme = 3;
        iNPPalgle = 1;
        ieratioalgle = 2;
    end

    % - Select the chosen beta map time-series
    ibetaalg=2;
    %betanow = squeeze(betaalgs(latidx,lonidx,:,ibetaalg));
    betanow = squeeze(nanmean(nanmean(squeeze(betaalgs(latidx,lonidx,:,ibetaalg)),2),1));
    
    % - Locate desired export time-series and regress against beta
    %expCfrankentsnorm1pt = squeeze(expCfrankents(latidx,lonidx,:)./nanmean(expCfrankents(latidx,lonidx,:))); 
    expCfrankentsnorm1pt = squeeze(nanmean(nanmean(squeeze(expCfrankents(latidx,lonidx,:)),2),1))./nanmean(nanmean(nanmean(expCfrankents(latidx,lonidx,:))));
    [bvseccfranken,bvsercfranken,bvserifranken,bvsenumtimepts,bvsepvalfranken]=tempcorrnanwithintercept(expCfrankentsnorm1pt,betanow);
    % me = most extreme
    %expCmetsnorm1pt = squeeze(expCalgs(latidx,lonidx,:,iNPPalgme,ieratioalgme))./nanmean(expCalgs(latidx,lonidx,:,iNPPalgme,ieratioalgme));
    expCmetsnorm1pt = squeeze(nanmean(nanmean(squeeze(expCalgs(latidx,lonidx,:,iNPPalgme,ieratioalgme)),2),1))./nanmean(nanmean(nanmean(expCalgs(latidx,lonidx,:,iNPPalgme,ieratioalgme))));
    [bvseccme,bvsercme,bvserime,bvsenumtimepts,bvsepvalme]=tempcorrnanwithintercept(expCmetsnorm1pt,betanow);
    % le = least extreme
    %expCletsnorm1pt = squeeze(expCalgs(latidx,lonidx,:,iNPPalgle,ieratioalgle))./nanmean(expCalgs(latidx,lonidx,:,iNPPalgle,ieratioalgle));
    expCletsnorm1pt = squeeze(nanmean(nanmean(squeeze(expCalgs(latidx,lonidx,:,iNPPalgle,ieratioalgle)),2),1))./nanmean(nanmean(nanmean(expCalgs(latidx,lonidx,:,iNPPalgle,ieratioalgle))));
    [bvseccle,bvsercle,bvserile,bvsenumtilepts,bvsepvalle]=tempcorrnanwithintercept(expCletsnorm1pt,betanow);
    
    subplot(2,2,2*ireg-1);
    tsbegin = 4; % 28 = Dec99 
    tsend = 77;
    alljans = 5:12:length(time(1:tsend));
    [ax,h1,h2] = plotyy(time(tsbegin:tsend),expCfrankentsnorm1pt(tsbegin:tsend),time(tsbegin:tsend),betanow(tsbegin:tsend));
    % remember that left and right y axis ticks must align 
    set(h1,'linewidth',2);set(h2,'linewidth',2);
    set(ax(1),'xlim',[time(tsbegin) time(tsend)],'xtick',time(alljans),'ylim',[expCmin expCmax],'ytick',[expCmin:0.5:expCmax]);
    set(ax(2),'xlim',[time(tsbegin) time(tsend)],'xtick',time(alljans),'ylim',[betamin betamax],'ytick',betaticks);
    datetick(ax(1),'x',12,'keeplimits','keepticks');
    datetick(ax(2),'x',12,'keeplimits','keepticks');
    ylabel('Norm export');
    %ylabel('\beta');
    %title([opts.betaalgarray{ibetaalg} 'beta vs. franken export;' 10 '(' num2str(lat1(latidx)) ',' num2str(lon1(lonidx)) '),' num2str(bvsercfranken)]);
    title([opts.betaalgarray{ibetaalg} 'beta vs. franken export;' 10 'btwn ' num2str(min(lat1(latidx))) ', ' num2str(max(lat1(latidx))) ' lat;' 10 num2str(min(lon1(lonidx))) ', ' num2str(max(lon1(lonidx))) ' lon' 10 num2str(bvsercfranken)]);
    hold on;
    %set(gca,'FontName','Arial','FontSize',20,'FontWeight','Bold');
    set(gca,'xgrid','on')
    %grid on;
    subplot(2,2,2*ireg);
    plot(expCfrankentsnorm1pt,betanow,'.','markersize',20,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor',[0.5 0.5 0.5]);
    set(gca,'ylim',[betamin betamax],'ytick',betaticks);
    xlabel('Normalized export');
    ylabel('\beta');
    xlim([expCmin expCmax]);
    ylim([betamin betamax]);
    hold on;
    plot([expCmin expCmax],bvsercfranken*[expCmin expCmax]+bvserifranken,'linewidth',2,'color','k');
    plot([expCmin expCmax],bvsercme*[expCmin expCmax]+bvserime,'linewidth',2,'color','k','linestyle','--');
    plot([expCmin expCmax],bvsercle*[expCmin expCmax]+bvserile,'linewidth',2,'color','k','linestyle','--');
    %title(['(' num2str(lat1(latidx)) ',' num2str(lon1(lonidx)) ');' num2str(bvsercme) ';' num2str(bvsercfranken) ';' num2str(bvsercle)]);
    title(['btwn ' num2str(min(lat1(latidx))) ', ' num2str(max(lat1(latidx))) ' lat;' 10 num2str(min(lon1(lonidx))) ', ' num2str(max(lon1(lonidx))) ' lon' 10 num2str(bvsercme) ';' num2str(bvsercfranken) ';' num2str(bvsercle)]);
    %set(gca,'FontName','Arial','FontSize',20,'FontWeight','Bold');

end % end tropical/subtropical plotting

end % end plot section

%--------------------------------------------------
% Kost and Guidi beta scatterplots versus one export at a time
% to find most and least extreme beta vs. export slopes
%--------------------------------------------------
plotsection=0;

if plotsection==1

% define either subtropical or tropical region here
latidx = latidxst;
lonidx = lonidxst;

for ibetaalg = 1:length(opts.betaalgarray)

    % - Select the chosen beta map time-series
    %betanow = squeeze(betaalgs(latidx,lonidx,:,ibetaalg));
    betanow = squeeze(nanmean(nanmean(squeeze(betaalgs(latidx,lonidx,:,ibetaalg)),2),1));

    f1=figure;
    set(f1,'color','white','units','inches','position',[0.5 0.5 12 10]);
    f2=figure;
    set(f2,'color','white','units','inches','position',[0.5 0.5 12 10]);
    isubplot=1;

    for iNPPalg = 1:length(opts.NPPalgarray)
        for ieratioalg = 1:length(opts.eratioalgarray)
    
            % - Select the chosen expC map time-series
            %expCnow = squeeze(expCalgs(latidx,lonidx,:,iNPPalg,ieratioalg))./nanmean(expCalgs(latidx,lonidx,:,iNPPalg,ieratioalg));
            expCnow = squeeze(nanmean(nanmean(squeeze(expCalgs(latidx,lonidx,:,iNPPalg,ieratioalg)),2),1))./nanmean(nanmean(nanmean(expCalgs(latidx,lonidx,:,iNPPalg,ieratioalg))));
        
            % - Generate temp regress coeff at chosen spatial grid pt
            [bvsetempcorrelcoeff,bvsetempregresscoeff,bvsetempregressint,bvsenumtimepts,bvsepval]=tempcorrnanwithintercept(expCnow,betanow);
            
            % - Looking at scatter plots at chosen spatial grid pt
            figure(f1);
            subplot(3,3,isubplot);
            plot(expCnow,betanow,'.','markersize',20,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor',[0.5 0.5 0.5]);
            %xlabel('Export (molC m^{-2} yr^{-1})');
            xlabel('Normalized export');
            ylabel('\beta');
            %xlim([0 2.5]);
            %ylim([3.5 6.5]);
            %title([opts.betaalgarray{ibetaalg} 'beta vs. ' opts.NPPalgarray{iNPPalg} 'NPP+' opts.eratioalgarray{ieratioalg} 'e-ratio;' 10 '(' num2str(lat1(latidx)) ',' num2str(lon1(lonidx)) '),' num2str(bvsetempregresscoeff)]);
            title([opts.betaalgarray{ibetaalg} 'beta vs. ' opts.NPPalgarray{iNPPalg} 'NPP+' opts.eratioalgarray{ieratioalg} 'e-ratio;' 10 'btwn ' num2str(min(lat1(latidx))) ', ' num2str(max(lat1(latidx))) ' lat;' 10 num2str(min(lon1(lonidx))) ', ' num2str(max(lon1(lonidx))) ' lon; ' num2str(bvsetempregresscoeff)]);
            hold on;
            plot([min(expCnow) max(expCnow)],bvsetempregresscoeff*[min(expCnow) max(expCnow)]+bvsetempregressint,...
                'linewidth',2,'color','k');
            set(gca,'FontName','Arial','FontSize',20,'FontWeight','Bold');
            %grid on;
            %f=figure;
            %set(f,'color','white','units','inches','position',[1 1 10 6]);
            %plotyy(1:160,squeeze(expC(latidx,lonidx,:)),1:160,squeeze(beta(latidx,lonidx,:)));

            % - Looking at time series at chosen spatial grid pt
            figure(f2);
            subplot(3,3,isubplot);
            plotyy(1:100,expCnow(1:100),1:100,betanow(1:100));
            %xlabel('Export (molC m^{-2} yr^{-1})');
            ylabel('Norm export');
            %ylabel('\beta');
            %xlim([0 2.5]);
            %ylim([3.5 6.5]);
            %title([opts.betaalgarray{ibetaalg} 'beta vs. ' opts.NPPalgarray{iNPPalg} 'NPP+' opts.eratioalgarray{ieratioalg} 'e-ratio;' 10 '(' num2str(lat1(latidx)) ',' num2str(lon1(lonidx)) '),' num2str(bvsetempregresscoeff)]);
            title([opts.betaalgarray{ibetaalg} 'beta vs. ' opts.NPPalgarray{iNPPalg} 'NPP+' opts.eratioalgarray{ieratioalg} 'e-ratio;' 10 'btwn ' num2str(min(lat1(latidx))) ', ' num2str(max(lat1(latidx))) ' lat;' 10 num2str(min(lon1(lonidx))) ', ' num2str(max(lon1(lonidx))) ' lon; ' num2str(bvsetempregresscoeff)]);
            hold on;
            %set(gca,'FontName','Arial','FontSize',20,'FontWeight','Bold');
            %grid on;

            isubplot = isubplot+1; 

        end % end eratio alg for loop
    end % end NPP alg for loop

end % end beta for loop

end % end plot section
