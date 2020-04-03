runsection=1;

if runsection==1

clear all;
close all;
addpath(genpath('/graid1/shirlleu/matlabroutines'));
addpath(genpath('/graid1/shirlleu/Pcycle_3D'));
fignum=0;

%%%%%%%%%%%%%%%%%
% MODES TO EDIT:
%%%%%%%%%%%%%%%%%

% circulation factors we want to plot
circfactorarray = [0.9 1.0 1.1];
circ1idx=2; % circfactorarray idx where circ factor is 1.0
numcircfactors = length(circfactorarray);

% names for legends
runnames = {'fboffK','allonKost'}; 
numruntypes = length(runnames); 

% Type of runs you are comparing - OCMIP, PRiSM NM, PRiSM TS
modeltype = {'PRiSM_TS','PRiSM_TS'};

% Feedback status of runs you are comparing
fbstatus = {'Pfboff','Pfbon'};
baserunidxK = 1; % cell idx of the baseline PSD fb off run in cell array fbstatus (Kost)
baserunidxG = 3; % cell idx of the baseline PSD fb off run in cell array fbstatus (Guidi)

% beta map for each run
betatype = {'Kost','Kost'};

% filename stem for saving pdfs
namestem = 'PRiSMTS_progprodtype_regionalblocks_PSDoffvson'

%%%%%%%%%%%%%%%%%
% OUTPUT FROM DIFFERENT RUNS TO EDIT:
%%%%%%%%%%%%%%%%%

for icirc=1:length(circfactorarray)

    %%%%%%%%%%%%%%%%%
    % KOSTADINOV
    %%%%%%%%%%%%%%%%%
    if circfactorarray(icirc)==1
    %---------------------
    % Baseline 1*circ
    %---------------------
        fname{icirc,1} = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/newtonmethod/spiralintoeq/prog/PSDfboff/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_Kostbetainit_PSDfboff_epsnolims_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_spiralintoeqiternum0.mat',circfactorarray(icirc));
    else
    %---------------------
    % PSD fb off
    %---------------------
        fname{icirc,1} = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/PSDfboff/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfboffKostbeta_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_1ov1000dt_ts500yrs_prevNMinit.mat',circfactorarray(icirc));
    end % end psd fb off circfactor if else

    %---------------------
    % Baseline 1*circ
    %---------------------
    if circfactorarray(icirc)==1
        fname{icirc,2} = fname{icirc,1};
    else
        %---------------------
        % PSD fb on, all regions turned on
        %---------------------
        fname{icirc,2} = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocksm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts500yrs_prevNMinit.mat',circfactorarray(icirc));
    end % end psd fb on circfactor if else

end

icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors] % don't want 1*circ

    %%%%%%%%%%%%%%%%%
    % KOSTADINOV
    %%%%%%%%%%%%%%%%%

    %---------------------
    % PSD fb on, all regions turned on, beta vs. export - 1stdev
    %---------------------
    fnamestd{icircstore,1} = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdbelowm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts300yrs_prevNMinit.mat',circfactorarray(icirc));

    %---------------------
    % PSD fb on, all regions turned on, beta vs. export + 1stdev
    %---------------------
    fnamestd{icircstore,2} = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdabovem1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts300yrs_prevNMinit.mat',circfactorarray(icirc));

    icircstore = icircstore + 1;

end

%%%%%%%%%%%%%%%%%
% load grid and observations - should be the same for all runs
%%%%%%%%%%%%%%%%%

% unpack output - should be the same for all runs
load('/graid1/shirlleu/PRiSM_GCM/sran+burial+ballast/newmatfiles/newtonmethod/BGRID_t+s+14c+cfc11_REDI_sol1_sigma0.10_pop2dop_ann_np3_PSDfboff_unifbeta_1.0circ.mat','output');
gridd = output.grid;
pobs = output.po4star;
M3d = output.M3d;
iocn = find(M3d==1);
isurf = find(M3d(:,:,1)==1);
m = length(iocn);
ns = length(isurf);
clear output;

% grid output - should be the same for all runs
Areat = gridd.Areat;
lat = gridd.yt;
totocnarea = sum(Areat(isurf)); % [m^2]
Areatocnonly = nan(size(Areat)); Areatocnonly(isurf) = Areat(isurf); 
Areatzonalsum = nansum(Areatocnonly,2); 

%%%%%%%%%%%%%%%%%
% prepare variables to store diff run output in
%%%%%%%%%%%%%%%%%

[ylen,xlen,zlen]=size(M3d);
SSexpCmap = zeros(ylen,xlen,numcircfactors,numruntypes);
SSexpCmapstd = zeros(ylen,xlen,numcircfactors-1,2*length(unique(betatype))); % 1-stdev above and below for each beta

%%%%%%%%%%%%%%%%%
% load output for all runs (column vectors)
%%%%%%%%%%%%%%%%%

for icirc=1:length(circfactorarray)
    for irun=1:numruntypes
        load(fname{icirc,irun},'output');
        if icirc==circ1idx % use PRiSM_NM ref circ baseline run
            SSexpCmap(:,:,icirc,irun) = output.expCmapSS;
            if irun==baserunidxK
                SStotexpCbaselineK = nansum(nansum(SSexpCmap(:,:,icirc,irun).*Areatocnonly));
            elseif irun==baserunidxG
                SStotexpCbaselineG = nansum(nansum(SSexpCmap(:,:,icirc,irun).*Areatocnonly));
            end
        else
            SSexpCmap(:,:,icirc,irun) = output.expCmapnow(:,:,end);
            totexpCovtime(:,icirc,irun) = output.totexpC;
            time(:,icirc,irun) = output.time;
        end
        clear output;
    end
end

for icircstore=1:length(circfactorarray)-1
    for irun=1:2*length(unique(betatype)) % 1-stdev above and below for each beta
        load(fnamestd{icircstore,irun},'output');
        SSexpCmapstd(:,:,icircstore,irun) = output.expCmapnow(:,:,end);
        totexpCovtimestd(:,icircstore,irun) = output.totexpC;
        timestd(:,icircstore,irun) = output.time;
        clear output;
    end
end

%%%%%%%%%%%%%%%%%
% calculate variables needed for plotting
%%%%%%%%%%%%%%%%%

% - Abs change in global mean export from baseline case
absch_expC = nan(numcircfactors-1,numruntypes);
itercirc=1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun=1:numruntypes
        if strcmp(betatype(irun),'Kost')
            absch_expC(itercirc,irun)=totexpCovtime(end,icirc,irun)-SStotexpCbaselineK;
        elseif strcmp(betatype(irun),'Guidi')
            absch_expC(itercirc,irun)=totexpCovtime(end,icirc,irun)-SStotexpCbaselineG;
        end
    end
    itercirc=itercirc+1;
end
absch_expC = absch_expC/totocnarea;

% - Rel change in global mean export from baseline case
relch_expC = nan(numcircfactors-1,numruntypes);
itercirc=1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun=1:numruntypes
        if strcmp(betatype(irun),'Kost')
            relch_expC(itercirc,irun)=(totexpCovtime(end,icirc,irun)-SStotexpCbaselineK)/SStotexpCbaselineK;
        elseif strcmp(betatype(irun),'Guidi')
            relch_expC(itercirc,irun)=(totexpCovtime(end,icirc,irun)-SStotexpCbaselineG)/SStotexpCbaselineG;
        end
    end
    itercirc=itercirc+1;
end

% - Abs change in global mean export from baseline case for stdev +/- runs
absch_expCstd = nan(numcircfactors-1,numruntypes);
itercirc=1;
for icirc=1:numcircfactors-1
    for irun=1:2*length(unique(betatype)) % 1-stdev above and below for each beta
        if strcmp(betatype(irun),'Kost')
            absch_expCstd(itercirc,irun)=totexpCovtimestd(end,icirc,irun)-SStotexpCbaselineK;
        elseif strcmp(betatype(irun),'Guidi')
            absch_expCstd(itercirc,irun)=totexpCovtimestd(end,icirc,irun)-SStotexpCbaselineG;
        end
    end
    itercirc=itercirc+1;
end
absch_expCstd = absch_expCstd/totocnarea;

% - Rel change in global mean export from baseline case for stdev +/- runs
relch_expCstd = nan(numcircfactors-1,numruntypes);
itercirc=1;
for icirc=1:numcircfactors-1
    for irun=1:2*length(unique(betatype))
        if strcmp(betatype(irun),'Kost')
            relch_expCstd(itercirc,irun)=(totexpCovtimestd(end,icirc,irun)-SStotexpCbaselineK)/SStotexpCbaselineK;
        elseif strcmp(betatype(irun),'Guidi')
            relch_expCstd(itercirc,irun)=(totexpCovtimestd(end,icirc,irun)-SStotexpCbaselineG)/SStotexpCbaselineG;
        end
    end
    itercirc=itercirc+1;
end

end % end runsection==1

%%%%%%%%%%%%%%%%%
% Create figures!!!
%%%%%%%%%%%%%%%%%

plotsection=1;

if plotsection==1

% - Time step toward SS w/ fb on/off fig
f=figure;
%set(f,'color','white','units','inches','position',[1 1 13 5]);
set(f,'color','white','units','inches','position',[1 1 9 5]);

fbofflinestyle='-';
fbonlinestyle='--';
linewidth=3;

subplot(121);
p(1)=plot([-10 0],[0 0],'k','linewidth',linewidth); hold on;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun=1:numruntypes
        if strcmp(betatype(irun),'Kost')
            h1=plot([0 time(:,icirc,irun)'],([SStotexpCbaselineK totexpCovtime(:,icirc,irun)']-SStotexpCbaselineK)/totocnarea);hold on;
            set(h1,'Color','b');
        elseif strcmp(betatype(irun),'Guidi')
            h1=plot([0 time(:,icirc,irun)'],([SStotexpCbaselineG totexpCovtime(:,icirc,irun)']-SStotexpCbaselineG)/totocnarea);hold on;
            set(h1,'Color','r');
        end
        if strcmp(fbstatus(irun),'Pfboff')
            set(h1,'LineStyle',fbofflinestyle,'LineWidth',linewidth);
        elseif strcmp(fbstatus(irun),'Pfbon')
            set(h1,'LineStyle',fbonlinestyle,'LineWidth',linewidth);
        end
    end
end
%xlim([-10 500]);
%ylim([-0.3 0.3]);
ylabel('\DeltaGlobal mean export [molC/m^2/yr]');
xlabel('Time before or after perturbation (years)');
set(gca,'fontsize',12,'fontweight','bold','YTick',-0.3:0.1:0.3);
%legend(' ',' ',' ',' ',' '); % just to get the legend elems on the fig
%legend(p,{'Baseline','Feedback off','Feedback on'});
%grid on;
    
for iweird = 1:2 % for some reason, gotta run this twice to get the hatching pattersn to look correct

    % - Abs change in global mean export from baseline case
    subplot(122);
    for icircstore = 1:numcircfactors-1 % circfactors excluding 1*circ
        for irun = 1:numruntypes
            if strcmp(betatype(irun),'Kost')
                if strcmp(fbstatus(irun),'Pfbon') % hatched fb on
                    barcolor='w';
                else
                    barcolor=[0 0 0.8];
                end
                ibeta = 1;
            elseif strcmp(betatype(irun),'Guidi')
                barcolor='r';
                ibeta = 2;
            end
            h1 = bar(ibeta,absch_expC(icircstore,irun),'facecolor',barcolor);hold on;
            if strcmp(fbstatus(irun),'Pfbon') % hatched fb on
                hatchfill2(h1,'hatchstyle','single','hatchangle',45,'hatchdensity',15,'hatchcolor',[0 0 0.8],'hatchlinewidth',2);
            end
            if strcmp(fbstatus(irun),'Pfbon') % fb on runs
                errors = abs(absch_expC(icircstore,irun) - absch_expCstd(icircstore,strcmp(betatype,betatype(irun))));
                eb = errorbar(ibeta,absch_expC(icircstore,irun),errors(1),errors(2),'color','k','linewidth',2); % errorbar(X,Y,L,U)
                errorbar_tick(eb,6);
            end
        end
    end
    %ylim([-0.3 0.3]);
    set(gca, 'XTick', 1, 'XTickLabel', {'Kost'},'fontsize',12,'fontweight','bold', 'YTick' ,-0.3:0.1:0.3);
    title('Abs \Delta in mean export from baseline','fontsize',10);
    %legend(' ',' ',' ',' ',' ',' ',' '); % just to get the legend elems on the fig
    %grid on;
    
    % - Rel change in global mean export from baseline case
    plotsp3=0;
    if plotsp3==1
        subplot(133);
        for icircstore = 1:numcircfactors-1 % circfactors excluding 1*circ
            for irun = 1:numruntypes
                if strcmp(betatype(irun),'Kost')
                    barcolor='b';
                    ibeta = 1;
                elseif strcmp(betatype(irun),'Guidi')
                    barcolor='r';
                    ibeta = 2;
                end
                h1 = bar(ibeta,relch_expC(icircstore,irun),'facecolor',barcolor);hold on;
                %if strcmp(fbstatus(irun),'Pfboff') % hatched fb off
                if strcmp(fbstatus(irun),'Pfbon') % hatched fb on
                    h2 = findobj(h1,'type','patch');
                    h3 = hatchfill(h2,'single',45,15,[1 1 1]); % angle, spacing, bgcolor
                    set(h3,'color',barcolor,'linewidth',2);
                end
                if strcmp(fbstatus(irun),'Pfbon') % fb on runs
                    errors = abs(relch_expC(icircstore,irun) - relch_expCstd(icircstore,strcmp(betatype,betatype(irun))));
                    eb = errorbar(ibeta,relch_expC(icircstore,irun),errors(1),errors(2),'color','k','linewidth',2); % errorbar(X,Y,L,U)
                    errorbar_tick(eb,6);
                end
            end
        end
        set(gca, 'XTick', 1, 'XTickLabel', {'Kost'},'fontsize',12,'fontweight','bold');
        %ylim([-0.3 0.3]);
        title('Rel change in global mean export from baseline case');
        %grid on;
    end % end if plotsp3

end % end iweird

end % end plotsection==1

print(f, ['yr500but300stdeverrors_fig5timestepglobalfb_Kbeta_PSDfbpaper.pdf'], '-dpdf', '-r300');
print(f, ['yr500but300stdeverrors_fig5timestepglobalfb_Kbeta_PSDfbpaper.png'], '-dpng', '-r300');
