addpath(genpath('/graid1/shirlleu/matlabroutines'));

%%%%%%%%%%%%%%%%%
% MODES TO EDIT:
%%%%%%%%%%%%%%%%%

% beta type I want to look at
betatype = 'Kost'; % choose Guidi or Kost

% circulation factors we want to plot
circfactorarray = [0.9 1.0];
circ09idx=1; % circfactorarray idx where circ factor is 1.0
circ1idx=2; % circfactorarray idx where circ factor is 1.0
numcircfactors = length(circfactorarray);

% define regions/number of regions 
pregs = load('/graid1/shirlleu/PRiSM_GCM/data/PRiSMregions.mat');
regmappgrid = pregs.R2d;
regnames = pregs.regnames(1:9); % pregs.regnames(1:9)
regnamesabbrev = pregs.regnamesabbrev(1:9);
% reorder Indian Ocean
regnamesabbrev = [regnamesabbrev(1:4) regnamesabbrev(9) regnamesabbrev(5:8)];
%f=figure; set(f,'color','white'); pcolor(regmappgrid); shading flat; colorbar;
numregions=length(regnames);
 
% names for legends
runnames = {'fboff','allon'}; 
numruntypes = length(runnames); 

% Type of runs you are comparing - OCMIP, PRiSM NM, PRiSM TS
modeltype = {'PRiSM_TS','PRiSM_TS'};

% Feedback status of runs you are comparing
fbstatus = {'Pfboff','Pfbon'};
baserunidx = 1; % cell idx of the baseline PSD fb off run in cell array fbstatus

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
        fname{icirc,1} = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/PSDfboff/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfboffKostbeta_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_1ov1000dt_ts100yrs_prevNMinit.mat',circfactorarray(icirc));
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
        fname{icirc,2} = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocksm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactorarray(icirc));
    end % end psd fb on circfactor if else
end % end circfactor for loop

icircstore=1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors] % don't want 1*circ

    %%%%%%%%%%%%%%%%%
    % KOSTADINOV
    %%%%%%%%%%%%%%%%%

    %---------------------
    % PSD fb on, all regions turned on, beta vs. export - 1stdev
    %---------------------
    fnamestd{icircstore,1} = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdbelowm1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactorarray(icirc));

    %---------------------
    % PSD fb on, all regions turned on, beta vs. export + 1stdev
    %---------------------
    fnamestd{icircstore,2} = sprintf('/graid1/shirlleu/PRiSM_GCM/flexible/newmatfiles/timestep/prog/slopeonlynormexpC_betachangescheme/regionalblocks/bvse_Gbetawt0Kbetawt1/BGRID_t+s+14c+cfc11_RAYLEIGH_sol1_sigma0.10_pop2dop_ann_np3_PSDfbonbvsen1_regblocks1stdabovem1slopeonlynexpC_Kostbetainit_epsbtwn2and6p5_%.2fcirc_BUR0_BAL0_slowsink0_varCP1_progprodtype_0.5mumax_recalcremin10yrs_1ov1000dt_ts100yrs_prevNMinit.mat',circfactorarray(icirc));

end % end circfactor for loop

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

expC = zeros(ylen,xlen,numcircfactors,numruntypes);
expCz = zeros(ylen,numcircfactors,numruntypes);
expCregs = zeros(numcircfactors,numruntypes,numregions+1); % +1 for global
expCstd = zeros(ylen,xlen,numcircfactors-1,2);
expCzstd = zeros(ylen,numcircfactors-1,2);
expCregsstd = zeros(numcircfactors-1,2,numregions+1); % +1 for global

po4 = zeros(length(iocn),numcircfactors,numruntypes);
Pmat = zeros(ylen,xlen,zlen,numcircfactors,numruntypes);
po4tc = zeros(ylen,xlen,numcircfactors,numruntypes);
po4tcz = zeros(ylen,numcircfactors,numruntypes);
po4tcregs = zeros(numcircfactors,numruntypes,numregions+1); % +1 for global

po4std = zeros(length(iocn),numcircfactors,numruntypes);
Pmatstd = zeros(ylen,xlen,zlen,numcircfactors,numruntypes);
po4tcstd = zeros(ylen,xlen,numcircfactors-1,2);
po4tczstd = zeros(ylen,numcircfactors-1,2);
po4tcregsstd = zeros(numcircfactors-1,2,numregions+1); % +1 for global

beta = zeros(ylen,xlen,numcircfactors,numruntypes);
betaz = zeros(ylen,numcircfactors,numruntypes);
betaregs = zeros(numcircfactors,numruntypes,numregions+1); % +1 for global

betastd = zeros(ylen,xlen,numcircfactors-1,2);
betazstd = zeros(ylen,numcircfactors-1,2);
betaregsstd = zeros(numcircfactors-1,2,numregions+1); % +1 for global

%%%%%%%%%%%%%%%%%
% load output for all runs (column vectors)
%%%%%%%%%%%%%%%%%

for icirc=1:numcircfactors
    for irun=1:numruntypes
        load(fname{icirc,irun},'output');
        if icirc==circ1idx % use PRiSM_NM ref circ baseline run
            expC(:,:,icirc,irun) = output.expCmapSS;
            po4(:,icirc,irun) = output.po4;
            load(fname{icirc,irun},'opts');
            if isfield(opts,'betamapinit')
                beta(:,:,icirc,irun) = opts.betamapinit;
            else
                betatemp = M3d(:,:,1)*NaN; betatemp(isurf) = opts.beta;
                beta(:,:,icirc,irun) = betatemp;
            end
        else
            expC(:,:,icirc,irun) = output.expCmapnow(:,:,end);
            po4(:,icirc,irun) = output.po4end;
            if isfield(output,'betamapnow') % PSD fb on
               beta(:,:,icirc,irun) = output.betamapnow(:,:,end);
            elseif isfield(output.opts,'spiraledintoeqfname') % PSD fb off, spiraled into eq beta map
               sp = load(output.opts.spiraledintoeqfname,'opts');
               beta(:,:,icirc,irun) = sp.opts.betamapinit;
            elseif isfield(output.opts,'initcondfname') % PSD fb off, init cond beta map
               ic = load(output.opts.initcondfname,'opts');
               beta(:,:,icirc,irun) = ic.opts.betamapinit;
            else
               beta(:,:,icirc,irun) = nan;
            end
        end
        clear output;
    end
end

for icircstore=1:numcircfactors-1
    for irun=1:2 % 1-stdev above and below for each beta
        load(fnamestd{icircstore,irun},'output');
        expCstd(:,:,icircstore,irun) = output.expCmapnow(:,:,end);
        po4std(:,icircstore,irun) = output.po4end;
        if isfield(output,'betamapnow') % PSD fb on
           betastd(:,:,icirc,irun) = output.betamapnow(:,:,end);
        elseif isfield(output.opts,'spiraledintoeqfname') % PSD fb off, spiraled into eq beta map
           sp = load(output.opts.spiraledintoeqfname,'opts');
           betastd(:,:,icirc,irun) = sp.opts.betamapinit;
        elseif isfield(output.opts,'initcondfname') % PSD fb off, init cond beta map
           ic = load(output.opts.initcondfname,'opts');
           betastd(:,:,icirc,irun) = ic.opts.betamapinit;
        else
           betastd(:,:,icirc,irun) = nan; 
        end
        clear output;
    end
end

%%%%%%%%%%%%%%%%%
% grid po4 and extract chosen po4 depth for all runs
%%%%%%%%%%%%%%%%%

[depth,depthidx]=min(abs(gridd.zt-200));

for icirc=1:length(circfactorarray)
    for irun=1:numruntypes
        Pmattemp = M3d*NaN; Pmattemp(iocn) = po4(:,icirc,irun);
        Pmat(:,:,:,icirc,irun) = Pmattemp;
    end
end
po4tc = squeeze(Pmat(:,:,depthidx,:,:));

for icircstore=1:numcircfactors-1
    for irun=1:2 % 1-stdev above and below for each beta
        Pmattemp = M3d*NaN; Pmattemp(iocn) = po4std(:,icirc,irun);
        Pmatstd(:,:,:,icircstore,irun) = Pmattemp;
    end
end
po4tcstd = squeeze(Pmatstd(:,:,depthidx,:,:));

%-------------------------------------------------------------------
% CALCULATE MEANS AND CHANGES IN EXPORT

calcexportsection=1;

if calcexportsection==1

%%%%%%%%%%%%%%%%%
% Calculate global, zonal, regional mean export for all runs
%%%%%%%%%%%%%%%%%
roundple = -3;

% - Calculate global, zonal, regional mean export values
for icirc=1:numcircfactors
    for irun=1:numruntypes
        expCz(:,icirc,irun) = nansum(expC(:,:,icirc,irun).*Areatocnonly,2)./Areatzonalsum;
        totexpC(icirc,irun) = nansum(nansum(expC(:,:,icirc,irun).*Areatocnonly)); 
        expCtemp = expC(:,:,icirc,irun);  
        for iregion=1:4
            expCregs(icirc,irun,iregion) = nansum(expCtemp(regmappgrid==iregion).*Areatocnonly(regmappgrid==iregion))./nansum(Areatocnonly(regmappgrid==iregion)); % regional exports
        end
        % Reorder Indian Ocean
        for iregion=5
            expCregs(icirc,irun,iregion) = nansum(expCtemp(regmappgrid==9).*Areatocnonly(regmappgrid==9))./nansum(Areatocnonly(regmappgrid==9)); % regional exports
        end
        for iregion=6:9
            expCregs(icirc,irun,iregion) = nansum(expCtemp(regmappgrid==iregion-1).*Areatocnonly(regmappgrid==iregion-1))./nansum(Areatocnonly(regmappgrid==iregion-1)); % regional exports
        end
        expCregs(icirc,irun,end) = totexpC(icirc,irun)./totocnarea; % global export
    end
end

% - Calculate global, zonal, regional mean export values for stdev +/- runs
for icircstore=1:numcircfactors-1
    for irun=1:2 % 1-stdev above and below for each beta
        expCzstd(:,icircstore,irun) = nansum(expCstd(:,:,icircstore,irun).*Areatocnonly,2)./Areatzonalsum;
        totexpCstd(icircstore,irun) = nansum(nansum(expCstd(:,:,icircstore,irun).*Areatocnonly)); 
        expCtemp = expCstd(:,:,icircstore,irun);  
        for iregion=1:4
            expCregsstd(icircstore,irun,iregion) = nansum(expCtemp(regmappgrid==iregion).*Areatocnonly(regmappgrid==iregion))./nansum(Areatocnonly(regmappgrid==iregion)); % regional exports
        end
        % Reorder Indian Ocean
        for iregion=5
            expCregsstd(icircstore,irun,iregion) = nansum(expCtemp(regmappgrid==9).*Areatocnonly(regmappgrid==9))./nansum(Areatocnonly(regmappgrid==9)); % regional exports
        end
        for iregion=6:9
            expCregsstd(icircstore,irun,iregion) = nansum(expCtemp(regmappgrid==iregion-1).*Areatocnonly(regmappgrid==iregion-1))./nansum(Areatocnonly(regmappgrid==iregion-1)); % regional exports
        end
        expCregsstd(icircstore,irun,end) = totexpCstd(icircstore,irun)./totocnarea; % global export
    end
end

% - Calculate absolute regional/global mean export changes from baseline
absch_expCregs = nan(numcircfactors-1,numruntypes,numregions+1); 
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun=1:numruntypes
        absch_expCregs(icircstore,irun,:)=squeeze(expCregs(icirc,irun,:)-expCregs(circ1idx,baserunidx,:));
    end
    icircstore = icircstore+1;
end
absch_expCregs = squeeze(absch_expCregs);

% - Calculate absolute regional/global mean export changes from baseline for stdev +/- runs
absch_expCregsstd = nan(numcircfactors-1,2,numregions+1);
for icircstore=1:numcircfactors-1
    for irun=1:2 % 1-stdev above and below for each beta
        absch_expCregsstd(icircstore,irun,:)=squeeze(expCregsstd(icircstore,irun,:)-expCregs(circ1idx,baserunidx,:));
    end
end
absch_expCregsstd = squeeze(absch_expCregsstd);

% - Calculate relative regional/global mean export changes from baseline
relch_expCregs = nan(numcircfactors-1,numruntypes,numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun=1:numruntypes
        relch_expCregs(icircstore,irun,:)=squeeze(roundn(expCregs(icirc,irun,:)-expCregs(circ1idx,baserunidx,:),roundple)./expCregs(circ1idx,baserunidx,:));
    end
    icircstore = icircstore+1;
end
relch_expCregs = squeeze(relch_expCregs);

% - Calculate relative regional/global mean export changes from baseline for stdev +/- runs
relch_expCregsstd = nan(numcircfactors-1,2,numregions+1);
for icircstore=1:numcircfactors-1
    for irun=1:2
        relch_expCregsstd(icircstore,irun,:)=squeeze(roundn(expCregsstd(icircstore,irun,:)-expCregs(circ1idx,baserunidx,:),roundple)./expCregs(circ1idx,baserunidx,:));
    end
end
relch_expCregsstd = squeeze(relch_expCregsstd);

% - Calculate regional/global mean feedback strength
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
fbst_expCregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1); 
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbst_expCregs(icircstore,irunstore,:)=squeeze(100*roundn(expCregs(icirc,baserunidx,:)-expCregs(icirc,ifbonrun,:),roundple)./(expCregs(icirc,baserunidx,:)-expCregs(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbst_expCregs = squeeze(fbst_expCregs);

% - Calculate regional/global mean feedback strength for stdev +/- runs
fbst_expCregsstd = nan(numcircfactors-1,2,numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun = 1:2
        fbst_expCregsstd(icircstore,irun,:)=squeeze(100*roundn(expCregs(icirc,baserunidx,:)-expCregsstd(icircstore,irun,:),roundple)./(expCregs(icirc,baserunidx,:)-expCregs(circ1idx,baserunidx,:)));
    end
    icircstore = icircstore+1;
end
fbst_expCregsstd = squeeze(fbst_expCregsstd);

end % end calcexportsection

%-------------------------------------------------------------------
% CALC MEANS AND CHANGES IN PO4 AT THERMOCLINE (200 M DEPTH)

calcpo4tcsection=1;

if calcpo4tcsection==1

%%%%%%%%%%%%%%%%%
% Calculate global, zonal, regional mean po4tc for all runs
%%%%%%%%%%%%%%%%%
roundplp = -3;

% - Calculate global, zonal, regional mean po4tc values
for icirc=1:numcircfactors
    for irun=1:numruntypes
        po4tcz(:,icirc,irun) = nansum(po4tc(:,:,icirc,irun).*Areatocnonly,2)./Areatzonalsum;
        totpo4tc(icirc,irun) = nansum(nansum(po4tc(:,:,icirc,irun).*Areatocnonly));
        po4tctemp = po4tc(:,:,icirc,irun);
        for iregion=1:4
            po4tcregs(icirc,irun,iregion) = nansum(po4tctemp(regmappgrid==iregion).*Areatocnonly(regmappgrid==iregion))./nansum(Areatocnonly(regmappgrid==iregion)); % regional po4tcs
        end
        % Reorder Indian Ocean
        for iregion=5
            po4tcregs(icirc,irun,iregion) = nansum(po4tctemp(regmappgrid==9).*Areatocnonly(regmappgrid==9))./nansum(Areatocnonly(regmappgrid==9)); % regional po4tcs
        end
        for iregion=6:9
            po4tcregs(icirc,irun,iregion) = nansum(po4tctemp(regmappgrid==iregion-1).*Areatocnonly(regmappgrid==iregion-1))./nansum(Areatocnonly(regmappgrid==iregion-1)); % regional po4tcs
        end
        po4tcregs(icirc,irun,end) = totpo4tc(icirc,irun)./totocnarea; % global po4tc
    end
end

% - Calculate global, zonal, regional mean po4tc values for stdev +/- runs
for icircstore=1:numcircfactors-1
    for irun=1:2 % 1-stdev above and below for each beta
        po4tczstd(:,icircstore,irun) = nansum(po4tcstd(:,:,icircstore,irun).*Areatocnonly,2)./Areatzonalsum;
        totpo4tcstd(icircstore,irun) = nansum(nansum(po4tcstd(:,:,icircstore,irun).*Areatocnonly));
        po4tctemp = po4tcstd(:,:,icircstore,irun);
        for iregion=1:4
            po4tcregsstd(icircstore,irun,iregion) = nansum(po4tctemp(regmappgrid==iregion).*Areatocnonly(regmappgrid==iregion))./nansum(Areatocnonly(regmappgrid==iregion)); % regional po4tcs
        end
        % Reorder Indian Ocean
        for iregion=5
            po4tcregsstd(icircstore,irun,iregion) = nansum(po4tctemp(regmappgrid==9).*Areatocnonly(regmappgrid==9))./nansum(Areatocnonly(regmappgrid==9)); % regional po4tcs
        end
        for iregion=6:9
            po4tcregsstd(icircstore,irun,iregion) = nansum(po4tctemp(regmappgrid==iregion-1).*Areatocnonly(regmappgrid==iregion-1))./nansum(Areatocnonly(regmappgrid==iregion-1)); % regional po4tcs
        end
        po4tcregsstd(icircstore,irun,end) = totpo4tcstd(icircstore,irun)./totocnarea; % global po4tc
    end
end

% - Calculate absolute regional/global mean po4tc changes from baseline
absch_po4tcregs = nan(numcircfactors-1,numruntypes,numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun=1:numruntypes
        absch_po4tcregs(icircstore,irun,:)=squeeze(po4tcregs(icirc,irun,:)-po4tcregs(circ1idx,baserunidx,:));
    end
    icircstore = icircstore+1;
end
absch_po4tcregs = squeeze(absch_po4tcregs);

% - Calculate absolute regional/global mean po4tc changes from baseline for stdev +/- runs
absch_po4tcregsstd = nan(numcircfactors-1,2,numregions+1);
for icircstore=1:numcircfactors-1
    for irun=1:2 % 1-stdev above and below for each beta
        absch_po4tcregsstd(icircstore,irun,:)=squeeze(po4tcregsstd(icircstore,irun,:)-po4tcregs(circ1idx,baserunidx,:));
    end
end
absch_po4tcregsstd = squeeze(absch_po4tcregsstd);

% - Calculate relative regional/global mean po4tc changes from baseline
relch_po4tcregs = nan(numcircfactors-1,numruntypes,numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun=1:numruntypes
        relch_po4tcregs(icircstore,irun,:)=squeeze(roundn(po4tcregs(icirc,irun,:)-po4tcregs(circ1idx,baserunidx,:),roundplp)./po4tcregs(circ1idx,baserunidx,:));
    end
    icircstore = icircstore+1;
end
relch_po4tcregs = squeeze(relch_po4tcregs);

% - Calculate relative regional/global mean po4tc changes from baseline for stdev +/- runs
relch_po4tcregsstd = nan(numcircfactors-1,2,numregions+1);
for icircstore=1:numcircfactors-1
    for irun=1:2
        relch_po4tcregsstd(icircstore,irun,:)=squeeze(roundn(po4tcregsstd(icircstore,irun,:)-po4tcregs(circ1idx,baserunidx,:),roundplp)./po4tcregs(circ1idx,baserunidx,:));
    end
end
relch_po4tcregsstd = squeeze(relch_po4tcregsstd);

% - Calculate regional/global mean feedback strength
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
fbst_po4tcregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbst_po4tcregs(icircstore,irunstore,:)=squeeze(100*roundn(po4tcregs(icirc,baserunidx,:)-po4tcregs(icirc,ifbonrun,:),roundplp)./(po4tcregs(icirc,baserunidx,:)-po4tcregs(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbst_po4tcregs = squeeze(fbst_po4tcregs);

% - Calculate regional/global mean feedback strength for stdev +/- runs
fbst_po4tcregsstd = nan(numcircfactors-1,2,numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun = 1:2
        fbst_po4tcregsstd(icircstore,irun,:)=squeeze(100*roundn(po4tcregs(icirc,baserunidx,:)-po4tcregsstd(icircstore,irun,:),roundplp)./(po4tcregs(icirc,baserunidx,:)-po4tcregs(circ1idx,baserunidx,:)));
    end
    icircstore = icircstore+1;
end
fbst_po4tcregsstd = squeeze(fbst_po4tcregsstd);

end % end calcpo4tc section

%%%%%%%%%%%%%%%%%
% Create figures!!!
%%%%%%%%%%%%%%%%%

plotsection=1;

if plotsection==1
%--------------------------------------------------
% FIGURE 6:
% Row 1: mean export and po4tc in the baseline case, zonal then regional
% Row 2: rel change in mean export from baseline, zonal then regional
% Row 3: rel change in mean po4tc from baseline, zonal then regional 
%--------------------------------------------------
fbofflinestyle='-';
fbonlinestyle='--';
fbofflinewidth=2;
fbonlinewidth=2;
%expCcol = 'r';
expCcol = [0 0 0.8];
%po4col = 'm';
po4col = [1 0.4 0];
fontwt = 'bold';
fontsize = 12;

fignum=1;
f=figure(fignum); set(f,'Color','White','Units','Inches','Position',[1 1 10 9.5]);

%---------ROW 1
if strcmp(betatype,'Kost'); po4max = 3.2; expCmax = 10; end
if strcmp(betatype,'Guidi'); po4max = 2.4; expCmax = 8; end

% - Zonally avged export and po4tc for baseline 1*circ case
subplot(421);
[AX,H1,H2]=plotyy(lat,expCz(:,circ1idx,baserunidx),lat,po4tcz(:,circ1idx,baserunidx));
set(H1,'Color',expCcol,'LineStyle','-','LineWidth',2);
set(H2,'Color',po4col,'LineStyle','-','Linewidth',2);
set(AX(1),'YColor',expCcol,'XLim',[-80 80],'YLim',[0 expCmax],'YTick',[0:expCmax/4:expCmax]);
set(AX(2),'YColor',po4col,'XLim',[-80 80],'YLim',[0 po4max],'YTick',[0:po4max/4:expCmax]);
ylabel(AX(1),'[molC m^2 yr^{-1}]','FontSize',fontsize,'FontWeight',fontwt);
ylabel(AX(2),'[mmolP m^{-3}]','FontSize',fontsize,'FontWeight',fontwt);
set(AX,'FontSize',fontsize,'FontWeight',fontwt,'FontName','Arial','TickDir','out');
xlabel('Latitude');
title('Zonal mean export, po4tc baseline case');

% - Regionally avged export and po4tc for baseline 1*circ case
barwidth=0.3;
offset = barwidth/2;
subplot(422);
[AX,H1,H2]=plotyy((1:numregions+1)-offset,squeeze(expCregs(circ1idx,baserunidx,:)),...
    (1:numregions+1)+offset,squeeze(po4tcregs(circ1idx,baserunidx,:)),...
    @(x,y) bar(x,y,barwidth), @(x,y) bar(x,y,barwidth));
%[AX,H1,H2] =plotyy((1:numregions+1)-offset,squeeze(expCregs(circ1idx,baserunidx,:)),...
%    (1:numregions+1)+offset,squeeze(po4tcregs(circ1idx,baserunidx,:)),'bar','bar');
set(H1,'FaceColor',expCcol);
set(H2,'FaceColor',po4col);
set(AX(1),'YColor',expCcol,'YLim',[0 expCmax],'YTick',[0:expCmax/4:expCmax],'XLim',[0.5 numregions+1.5],'XTick',1:numregions+1,'XTickLabel',[regnamesabbrev {'glob'}]);
set(AX(2),'YColor',po4col,'YLim',[0 po4max],'YTick',[0:po4max/4:po4max],'XLim',[0.5 numregions+1.5],'XTick',[]);
set(AX,'FontSize',fontsize,'FontWeight',fontwt,'FontName','Arial','TickDir','out');
ylabel(AX(1),'[molC m^2 yr^{-1}]','FontSize',fontsize,'FontWeight',fontwt);
ylabel(AX(2),'[mmolP m^{-3}]','FontSize',fontsize,'FontWeight',fontwt);
title('Reg mean po4tc baseline case');

%---------ROW 3
% - Relative change in zonal mean export from baseline case
subplot(425);
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun=1:numruntypes
        h1=plot(lat,100*roundn(expCz(:,icirc,irun)-expCz(:,circ1idx,baserunidx),roundple)...
            ./expCz(:,circ1idx,baserunidx),'color',expCcol); hold on;
        if strcmp(fbstatus(irun),'Pfboff')
            set(h1,'LineStyle',fbofflinestyle,'LineWidth',fbofflinewidth);
        elseif strcmp(fbstatus(irun),'Pfbon')
            set(h1,'LineStyle',fbonlinestyle,'LineWidth',fbonlinewidth);
        end
    end
end
xlim([-80 80]);
ylim([-17 1]);
%ylim([-18 0.5]);
xlabel('Lat');
ylabel('[%]');
title('Rel change in zonal mean export from baseline case');
%set(gca,'xgrid','on');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Rel change in regionally avged export from baseline case
barwidth = 0.8;
ebwidth = 100; % 1/ebwidth is the width, technically
for iweird = 1:2
    subplot(426);
    solidbars = bar(1:numregions+1,100*relch_expCregs(baserunidx,:)',barwidth);hold on;
    hatchedbars = bar(1:numregions+1,100*relch_expCregs(baserunidx+1,:)',barwidth,'facecolor','white');
    colornow = expCcol; formatbars_newfig6;
    set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
    ylabel('[%]');
    ylim([-17 1]);
    xlim([0.5 numregions+1.5]);
    title('Rel change in reg mean export from baseline case');
end
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

%---------ROW 2
% - Relative change in zonal mean po4tc from baseline case
subplot(423);
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun=1:numruntypes
        h1=plot(lat,100*roundn(po4tcz(:,icirc,irun)-po4tcz(:,circ1idx,baserunidx),roundplp)...
            ./po4tcz(:,circ1idx,baserunidx),'color',po4col); hold on;
        if strcmp(fbstatus(irun),'Pfboff')
            set(h1,'LineStyle',fbofflinestyle,'LineWidth',fbofflinewidth);
        elseif strcmp(fbstatus(irun),'Pfbon')
            set(h1,'LineStyle',fbonlinestyle,'LineWidth',fbonlinewidth);
        end
    end
end
plot(lat,zeros(size(lat)),'k');
xlim([-80 80]);
ylim([-4.2 2]);
xlabel('Lat');
ylabel('[%]');
title('Rel change in zonal mean po4tc from baseline case');
%set(gca,'xgrid','on');
set(gca,'FontName','Arial','FontSize',12,'FontWeight',fontwt,'TickDir','out');

% - Rel change in regionally avged po4tc from baseline case
barwidth = 0.8;
ebwidth = 100; % 1/ebwidth is the width, technically
for iweird = 1:2
    subplot(424);
    solidbars = bar(1:numregions+1,100*relch_po4tcregs(baserunidx,:)',barwidth);hold on;
    hatchedbars = bar(1:numregions+1,100*relch_po4tcregs(baserunidx+1,:)',barwidth,'facecolor','white');
    colornow = po4col; formatbars_newfig6;
    set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
    ylabel('[%]');
    ylim([-4.2 2]);
    xlim([0.5 numregions+1.5]);
    title('Rel change in reg mean po4tc from baseline case');
end
set(gca,'FontName','Arial','FontSize',12,'FontWeight',fontwt,'TickDir','out');

%---------ROW 4
fbofflinestyle='-';
fbonlinestyle='--';
fbofflinewidth=2;
fbonlinewidth=2;
expCcol = [0 0 0.8];
po4col = [1 0.4 0];
fontwt = 'bold';
fontsize = 12;
deltacircovcirc=-0.1;

% - Zonal mean fb strength: (fboff-fbon)/(fboff-baseline)
subplot(427);
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCz(:,icirc,baserunidx)-expCz(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(expCz(:,icirc,baserunidx)-expCz(:,circ1idx,baserunidx)),...
            'Color',expCcol,'LineStyle','-','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            (roundn(po4tcz(:,icirc,baserunidx)-po4tcz(:,icirc,fbonrunidx(ifbonrun)),roundplp)...
            ./po4tcz(:,circ1idx,baserunidx))./ ...
            (deltacircovcirc+roundn(po4tcz(:,icirc,baserunidx)-po4tcz(:,circ1idx,baserunidx),roundplp)...
            ./po4tcz(:,circ1idx,baserunidx)), ...
            'Color',po4col,'LineStyle','-','LineWidth',2); hold on;
    end
end
xlabel('Lat');
ylabel('[%]');
title('[(po4 fboff-fbon)/baseline]/[fboff deltaw/w+deltaP/P] [%]');
xlim([-80 80]);
ylim([0 28]);
%legend('expC','po4tc');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Regional mean fb strength: (fboff-fbon)/(fboff-baseline)
ploterrbars=0;
barwidth=1;
offset = 0.125;
ebwidth = 100; % 1/ebwidth is the width, technically
subplot(428);
fbst_po4tcregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbst_po4tcregs(icircstore,irunstore,:)=squeeze(100*(roundn(po4tcregs(icirc,baserunidx,:)-po4tcregs(icirc,ifbonrun,:),roundplp)./po4tcregs(circ1idx,baserunidx,:))./ ...
            (deltacircovcirc+roundn(po4tcregs(icirc,baserunidx,:)-po4tcregs(circ1idx,baserunidx,:),roundplp)./po4tcregs(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbst_po4tcregs = squeeze(fbst_po4tcregs);
b=bar(1:numregions+1,[fbst_expCregs fbst_po4tcregs],barwidth);hold on;
set(b(1),'facecolor',expCcol);
set(b(2),'facecolor',po4col);
if ploterrbars==1
    fbsterrors = abs(bsxfun(@minus,fbst_expCregs',fbst_expCregsstd));
    eb2 = errorbar((1:numregions+1)-offset,fbst_expCregs',fbsterrors(1,:)',fbsterrors(2,:)','linestyle','none','color','k','linewidth',2); % errorbar(X,Y,L,U)
    errorbar_tick(eb2,ebwidth);
end
set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
ylabel('[%]');
ylim([0 28]);
xlim([0.5 numregions+1.5]);
title('(fboff-fbon)/(fboff-baseline) [%]');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

saveas(f,'newfig6.png')
end % end plotsection

disp('pause')
pause

%--------------------------------------------------
% FIGURE 7:
% Row 1:
% Row 2:
% Row 3:
%--------------------------------------------------
deltacircovcirc=-0.1;

plotsection=1;

if plotsection==1

fignum=fignum+1;
f=figure(fignum); set(f,'Color','White','Units','Inches','Position',[1 1 10 9.5]);
fbofflinestyle='-';
fbonlinestyle='--';
fbofflinewidth=2;
fbonlinewidth=2;
expCcol = [0 0 0.8];
po4col = [1 0.4 0];
fontwt = 'bold';
fontsize = 12;

% - Zonal mean rel change diff: (fboff-fbon)/baseline
subplot(421);
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCz(:,icirc,baserunidx)-expCz(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./expCz(:,circ1idx,baserunidx),...
            'Color',expCcol,'LineStyle','-','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(po4tcz(:,icirc,baserunidx)-po4tcz(:,icirc,fbonrunidx(ifbonrun)),roundplp)...
            ./po4tcz(:,circ1idx,baserunidx),...
            'Color',po4col,'LineStyle','-','LineWidth',2); hold on;
    end
end
xlabel('Lat');
ylabel('[%]');
title('(fboff-fbon)/baseline [%]');
xlim([-80 80]);
ylim([-3.5 0]);
legend('expC','po4tc');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Regional mean rel change diff: (fboff-fbon)/baseline
subplot(422);
barwidth = 1;
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
fbonminusoffrel_po4tcregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbonminusoffrel_po4tcregs(icircstore,irunstore,:)=squeeze(100*roundn(po4tcregs(icirc,baserunidx,:)-po4tcregs(icirc,ifbonrun,:),roundplp)./po4tcregs(circ1idx,baserunidx,:));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbonminusoffrel_po4tcregs = squeeze(fbonminusoffrel_po4tcregs);
fbonminusoffrel_expCregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbonminusoffrel_expCregs(icircstore,irunstore,:)=squeeze(100*roundn(expCregs(icirc,baserunidx,:)-expCregs(icirc,ifbonrun,:),roundple)./expCregs(circ1idx,baserunidx,:));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbonminusoffrel_expCregs = squeeze(fbonminusoffrel_expCregs);
b=bar(1:numregions+1,[fbonminusoffrel_expCregs fbonminusoffrel_po4tcregs],barwidth);hold on;
set(b(1),'facecolor',expCcol);
set(b(2),'facecolor',po4col);
set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
ylabel('[%]');
ylim([-3.5 0]);
xlim([0.5 numregions+1.5]);
title('(fboff-fbon)/(baseline) [%]');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Zonal mean fb strength: (fboff-fbon)/(fboff-baseline)
subplot(425);
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCz(:,icirc,baserunidx)-expCz(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(expCz(:,icirc,baserunidx)-expCz(:,circ1idx,baserunidx)),...
            'Color',expCcol,'LineStyle','-','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat', 100*...
            (roundn(po4tcz(:,icirc,baserunidx)-po4tcz(:,icirc,fbonrunidx(ifbonrun)),roundplp)...
            ./po4tcz(:,circ1idx,baserunidx))./ ...
            (roundn(expCz(:,icirc,baserunidx)-expCz(:,circ1idx,baserunidx),roundple)...
            ./expCz(:,circ1idx,baserunidx)), ...
            'Color',po4col,'LineStyle','-','LineWidth',2); hold on;
    end
end
xlabel('Lat');
ylabel('[%]');
title('((po4 fboff-fbon)/baseline) / ((expC fboff-baseline)/baseline) [%]');
xlim([-80 80]);
ylim([0 28]);
%legend('expC','po4tc');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Regional mean fb strength: (fboff-fbon)/(fboff-baseline)
ploterrbars=0;
barwidth=1;
offset = 0.125;
ebwidth = 100; % 1/ebwidth is the width, technically
subplot(426);
fbst_po4tcregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbst_po4tcregs(icircstore,irunstore,:)=squeeze(100*(roundn(po4tcregs(icirc,baserunidx,:)-po4tcregs(icirc,ifbonrun,:),roundplp)./po4tcregs(circ1idx,baserunidx,:))./ ...
            (roundn(expCregs(icirc,baserunidx,:)-expCregs(circ1idx,baserunidx,:),roundple)./expCregs(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbst_po4tcregs = squeeze(fbst_po4tcregs);
b=bar(1:numregions+1,[fbst_expCregs fbst_po4tcregs],barwidth);hold on;
set(b(1),'facecolor',expCcol);
set(b(2),'facecolor',po4col);
if ploterrbars==1
    fbsterrors = abs(bsxfun(@minus,fbst_expCregs',fbst_expCregsstd));
    eb2 = errorbar((1:numregions+1)-offset,fbst_expCregs',fbsterrors(1,:)',fbsterrors(2,:)','linestyle','none','color','k','linewidth',2); % errorbar(X,Y,L,U)
    errorbar_tick(eb2,ebwidth);
end
set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
ylabel('[%]');
ylim([0 28]);
xlim([0.5 numregions+1.5]);
title('(fboff-fbon)/(fboff-baseline) [%]');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Zonal mean fb strength: (fboff-fbon)/(fboff-baseline)
subplot(423);
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCz(:,icirc,baserunidx)-expCz(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(expCz(:,icirc,baserunidx)-expCz(:,circ1idx,baserunidx)),...
            'Color',expCcol,'LineStyle','-','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            (roundn(po4tcz(:,icirc,baserunidx)-po4tcz(:,icirc,fbonrunidx(ifbonrun)),roundplp)...
            ./po4tcz(:,circ1idx,baserunidx))./ ...
            (deltacircovcirc+roundn(po4tcz(:,icirc,baserunidx)-po4tcz(:,circ1idx,baserunidx),roundplp)...
            ./po4tcz(:,circ1idx,baserunidx)), ...
            'Color',po4col,'LineStyle','-','LineWidth',2); hold on;
    end
end
xlabel('Lat');
ylabel('[%]');
title('[(po4 fboff-fbon)/baseline]/[fboff deltaw/w+deltaP/P] [%]');
xlim([-80 80]);
ylim([0 28]);
%legend('expC','po4tc');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Regional mean fb strength: (fboff-fbon)/(fboff-baseline)
ploterrbars=0;
barwidth=1;
offset = 0.125;
ebwidth = 100; % 1/ebwidth is the width, technically
subplot(424);
fbst_po4tcregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbst_po4tcregs(icircstore,irunstore,:)=squeeze(100*(roundn(po4tcregs(icirc,baserunidx,:)-po4tcregs(icirc,ifbonrun,:),roundplp)./po4tcregs(circ1idx,baserunidx,:))./ ...
            (deltacircovcirc+roundn(po4tcregs(icirc,baserunidx,:)-po4tcregs(circ1idx,baserunidx,:),roundplp)./po4tcregs(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbst_po4tcregs = squeeze(fbst_po4tcregs);
b=bar(1:numregions+1,[fbst_expCregs fbst_po4tcregs],barwidth);hold on;
set(b(1),'facecolor',expCcol);
set(b(2),'facecolor',po4col);
if ploterrbars==1
    fbsterrors = abs(bsxfun(@minus,fbst_expCregs',fbst_expCregsstd));
    eb2 = errorbar((1:numregions+1)-offset,fbst_expCregs',fbsterrors(1,:)',fbsterrors(2,:)','linestyle','none','color','k','linewidth',2); % errorbar(X,Y,L,U)
    errorbar_tick(eb2,ebwidth);
end
set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
ylabel('[%]');
ylim([0 28]);
xlim([0.5 numregions+1.5]);
title('(fboff-fbon)/(fboff-baseline) [%]');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

end % end plotsection

%--------------------------------------------------
% Figuring out noisy points #1
%--------------------------------------------------
expCfboffreldiff=100*roundn(expC(:,:,1,baserunidx)-expC(:,:,circ1idx,baserunidx),roundple)./expC(:,:,circ1idx,baserunidx);
po4tcfboffonreldiff = 100*roundn(po4tc(:,:,1,baserunidx)-po4tc(:,:,1,2),roundplp)...
    ./po4tc(:,:,circ1idx,baserunidx); 
po4fbst=100*po4tcfboffonreldiff./expCfboffreldiff;

runsection=0;

if runsection==1

if strcmp(betatype,'Guidi') 
    threshold=0.2;
    slim=-40;
    nlim=55;
elseif strcmp(betatype,'Kost')
    threshold=0.4; 
    slim=-80;
    nlim=55;
end

expCfboffreldiffbi=nan(size(expCfboffreldiff));
expCfboffreldiffbi(abs(expCfboffreldiff)<threshold&(gridd.YT>nlim|gridd.YT<slim))=1;
expCfboffreldiffbi(abs(expCfboffreldiff)>=threshold)=0;

po4fbstbi=po4fbst;
po4fbstbi(abs(expCfboffreldiff)<threshold&(gridd.YT>nlim|gridd.YT<slim))=nan;

rdylbu = flipud(cbrewer('div','RdYlBu',100,'linear'));

% % - Plot expCfboffreldiff and some
% fignum=fignum+1;
% f=figure(fignum);
% set(f,'color','white','units','inches','position',[1 1 10 7.5]);
% colormap(rdylbu);
% subplot(321);
% pcolor(gridd.xt,gridd.yt,expC(:,:,1,baserunidx));
% shading flat; colorbar;
% title('expC 0.9*circ fb off [molC/m^2/yr]');
% subplot(322);
% pcolor(gridd.xt,gridd.yt,expC(:,:,circ1idx,baserunidx));
% shading flat; colorbar;
% title('expC 1*circ [molC/m^2/yr]');
% subplot(323);
% pcolor(gridd.xt,gridd.yt,expC(:,:,1,baserunidx)-expC(:,:,circ1idx,baserunidx));
% caxis([-0.1 0.1]);
% shading flat; colorbar;
% title('expC (0.9*circ fboff - 1*circ) [molC/m^2/yr]');
% subplot(324);
% pcolor(gridd.xt,gridd.yt,expCfboffreldiff);
% caxis([-10 10]);
% shading flat; colorbar;
% title('expC (0.9*circ fboff - 1*circ)/(1*circ) [%]');
% subplot(325);
% pcolor(gridd.xt,gridd.yt,expCfboffreldiffbi);
% %caxis([-0.1 0.1]);
% shading flat; colorbar;
% title(['1 where abs((0.9*circ fboff - 1*circ)/1*circ)<' num2str(threshold) '%']);
% 
% % - Plot po4tcfboffonreldiff and some
% fignum=fignum+1;
% f=figure(fignum);
% set(f,'color','white','units','inches','position',[1 1 10 7.5]);
% colormap(rdylbu);
% subplot(321);
% pcolor(gridd.xt,gridd.yt,po4tc(:,:,1,baserunidx));
% shading flat; colorbar;
% title('po4tc 0.9*circ fb off [mmol/m^3]');
% subplot(322);
% pcolor(gridd.xt,gridd.yt,po4tc(:,:,circ1idx,baserunidx));
% shading flat; colorbar;
% title('po4tc 1*circ [mmol/m^3]');
% subplot(323);
% pcolor(gridd.xt,gridd.yt,po4tc(:,:,1,baserunidx)-po4tc(:,:,1,2));
% caxis([-0.05 0.05]);
% shading flat; colorbar;
% title('po4tc (0.9*circ fboff - fbon) [mmol/m^3]');
% subplot(324);
% pcolor(gridd.xt,gridd.yt,po4tcfboffonreldiff);
% caxis([-2 2]);
% shading flat; colorbar;
% title('po4tc (0.9*circ fboff - fbon)/(1*circ) [%]');

% - Plot expCfboffreldiff, po4tcfboffonreldiff,
% their ratio, where points are masked out,
% the new zonal mean values excluding the masked out points (with the old zonal mean values),
% and the new regional mean values excluding the masked out points
fignum=fignum+1;
f=figure(fignum);
set(f,'color','white','units','inches','position',[1 1 10 9.5]);
colormap(rdylbu);
subplot(421);
pcolor(gridd.xt,gridd.yt,expCfboffreldiff);
caxis([-20 20]);
shading flat; colorbar;
title('expC (0.9*circ fboff - 1*circ)/(1*circ) [%]');
subplot(422);
pcolor(gridd.xt,gridd.yt,po4tcfboffonreldiff);
caxis([-5 5]);
shading flat; colorbar;
title('po4tc (0.9*circ fboff - fbon)/(1*circ) [%]');
subplot(423);
pcolor(gridd.xt,gridd.yt,expCfboffreldiffbi);
caxis([-1 1]);
shading flat; colorbar;
title(['1 where abs((expC 0.9*circ fboff - 1*circ)/1*circ)<' num2str(threshold) '%']);
subplot(424);
pcolor(gridd.xt,gridd.yt,po4fbst);
caxis([-500 500]);
shading flat; colorbar;
title('[(po4 fboff-fbon)/baseline]/[(expC fboff-baseline)/baseline] [%]');
subplot(425);
pcolor(gridd.xt,gridd.yt,po4fbstbi);
caxis([-500 500]);
shading flat; colorbar;
title('[(po4 fboff-fbon)/baseline]/[(expC fboff-baseline)/baseline] [%]');

%----------- BEGIN RECALCULATIONS #1
% Need recalculating of:
% po4tcz, po4tcregs, expCz, expCregs, fbst_expCregs, fbst_expCregsstd 

% - Recalculate po4tcz, po4tcregs
po4tczalt = nan(size(po4tcz));
po4tcregsalt = nan(size(po4tcregs));
totpo4tcalt = nan(size(totpo4tc));
for icirc=1:numcircfactors
    for irun=1:numruntypes
        po4tctemp = po4tc(:,:,icirc,irun);
        po4tctemp(abs(expCfboffreldiff)<threshold&(gridd.YT>nlim|gridd.YT<slim)) = nan; 
        totpo4tcalt(icirc,irun) = nansum(nansum(po4tctemp.*Areatocnonly));
        po4tczalt(:,icirc,irun) = nansum(po4tctemp.*Areatocnonly,2)./Areatzonalsum;
        for iregion=1:4
            po4tcregsalt(icirc,irun,iregion) = nansum(po4tctemp(regmappgrid==iregion).*Areatocnonly(regmappgrid==iregion))./nansum(Areatocnonly(regmappgrid==iregion)); % regional po4tcs
        end
        % Reorder Indian Ocean
        for iregion=5
            po4tcregsalt(icirc,irun,iregion) = nansum(po4tctemp(regmappgrid==9).*Areatocnonly(regmappgrid==9))./nansum(Areatocnonly(regmappgrid==9)); % regional po4tcs
        end
        for iregion=6:9
            po4tcregsalt(icirc,irun,iregion) = nansum(po4tctemp(regmappgrid==iregion-1).*Areatocnonly(regmappgrid==iregion-1))./nansum(Areatocnonly(regmappgrid==iregion-1)); % regional po4tcs
        end
        po4tcregsalt(icirc,irun,end) = totpo4tcalt(icirc,irun)./totocnarea; % global po4tc
    end
end

% - Recalculate expCz, expCregs
expCzalt = nan(size(expCz));
expCregsalt = nan(size(expCregs));
totexpCalt = nan(size(totexpC));
for icirc=1:numcircfactors
    for irun=1:numruntypes
        expCtemp = expC(:,:,icirc,irun);
        expCtemp(abs(expCfboffreldiff)<threshold&(gridd.YT>nlim|gridd.YT<slim)) = nan;
        totexpCalt(icirc,irun) = nansum(nansum(expCtemp.*Areatocnonly));
        expCzalt(:,icirc,irun) = nansum(expCtemp.*Areatocnonly,2)./Areatzonalsum;
        for iregion=1:4
            expCregsalt(icirc,irun,iregion) = nansum(expCtemp(regmappgrid==iregion).*Areatocnonly(regmappgrid==iregion))./nansum(Areatocnonly(regmappgrid==iregion)); % regional exports
        end
        % Reorder Indian Ocean
        for iregion=5
            expCregsalt(icirc,irun,iregion) = nansum(expCtemp(regmappgrid==9).*Areatocnonly(regmappgrid==9))./nansum(Areatocnonly(regmappgrid==9)); % regional exports
        end
        for iregion=6:9
            expCregsalt(icirc,irun,iregion) = nansum(expCtemp(regmappgrid==iregion-1).*Areatocnonly(regmappgrid==iregion-1))./nansum(Areatocnonly(regmappgrid==iregion-1)); % regional exports
        end
        expCregsalt(icirc,irun,end) = totexpCalt(icirc,irun)./totocnarea; % global export
    end
end

% - Recalculate fbst_expCregs, fbst_expCregsstd 
expCzstdalt = nan(size(expCzstd));
totexpCstdalt = nan(size(totexpCstd));
expCregsstdalt = nan(size(expCregsstd));
for icircstore=1:numcircfactors-1
    for irun=1:2 % 1-stdev above and below for each beta
        expCstdtemp = expCstd(:,:,icircstore,irun);
        expCstdtemp(abs(expCfboffreldiff)<threshold&(gridd.YT>nlim|gridd.YT<slim)) = nan;
        expCzstdalt(:,icircstore,irun) = nansum(expCstdtemp.*Areatocnonly,2)./Areatzonalsum;
        totexpCstdalt(icircstore,irun) = nansum(nansum(expCstdtemp.*Areatocnonly));
        for iregion=1:4
            expCregsstdalt(icircstore,irun,iregion) = nansum(expCstdtemp(regmappgrid==iregion).*Areatocnonly(regmappgrid==iregion))./nansum(Areatocnonly(regmappgrid==iregion)); % regional exports
        end
        % Reorder Indian Ocean
        for iregion=5
            expCregsstdalt(icircstore,irun,iregion) = nansum(expCstdtemp(regmappgrid==9).*Areatocnonly(regmappgrid==9))./nansum(Areatocnonly(regmappgrid==9)); % regional exports
        end
        for iregion=6:9
            expCregsstdalt(icircstore,irun,iregion) = nansum(expCstdtemp(regmappgrid==iregion-1).*Areatocnonly(regmappgrid==iregion-1))./nansum(Areatocnonly(regmappgrid==iregion-1)); % regional exports
        end
        expCregsstdalt(icircstore,irun,end) = totexpCstdalt(icircstore,irun)./totocnarea; % global export
    end
end

fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
fbst_expCregsalt = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbst_expCregsalt(icircstore,irunstore,:)=squeeze(100*roundn(expCregsalt(icirc,baserunidx,:)-expCregsalt(icirc,ifbonrun,:),roundple)./(expCregsalt(icirc,baserunidx,:)-expCregsalt(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbst_expCregsalt = squeeze(fbst_expCregsalt);

fbst_expCregsstdalt = nan(numcircfactors-1,2,numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun = 1:2
        fbst_expCregsstdalt(icircstore,irun,:)=squeeze(100*roundn(expCregsalt(icirc,baserunidx,:)-expCregsstdalt(icircstore,irun,:),roundple)./(expCregsalt(icirc,baserunidx,:)-expCregsalt(circ1idx,baserunidx,:)));
    end
    icircstore = icircstore+1;
end
fbst_expCregsstdalt = squeeze(fbst_expCregsstdalt);
%----------- FINISH RECALCULATIONS #1

% - Continue plotting expCfboffreldiff, po4tcfboffonreldiff,
% their ratio, where points are masked out,
% the new zonal mean values excluding the masked out points (with the old zonal mean values),
% and the new regional mean values excluding the masked out points

% Zonal mean fb strength: (fboff-fbon)/(fboff-baseline)
%expCcol = 'r';
expCcol = [0 0 0.8];
%po4col = 'm';
po4col = [1 0.4 0];
fontwt = 'normal';
subplot(427);
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCzalt(:,icirc,baserunidx)-expCzalt(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(expCzalt(:,icirc,baserunidx)-expCzalt(:,circ1idx,baserunidx)),...
            'Color',expCcol,'LineStyle','-','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat', 100*...
            (roundn(po4tczalt(:,icirc,baserunidx)-po4tczalt(:,icirc,fbonrunidx(ifbonrun)),roundplp)...
            ./(po4tczalt(:,circ1idx,baserunidx)))./ ...
            (roundn(expCzalt(:,icirc,baserunidx)-expCzalt(:,circ1idx,baserunidx),roundple)...
            ./expCzalt(:,circ1idx,baserunidx)), ...
            po4col,'LineStyle','-','LineWidth',2); hold on;
    end
end
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCz(:,icirc,baserunidx)-expCz(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(expCz(:,icirc,baserunidx)-expCz(:,circ1idx,baserunidx)),...
            'Color',expCcol,'LineStyle','--','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat', 100*...
            (roundn(po4tcz(:,icirc,baserunidx)-po4tcz(:,icirc,fbonrunidx(ifbonrun)),roundplp)...
            ./(po4tcz(:,circ1idx,baserunidx)))./ ...
            (roundn(expCz(:,icirc,baserunidx)-expCz(:,circ1idx,baserunidx),roundple)...
            ./expCz(:,circ1idx,baserunidx)), ...
            po4col,'LineStyle','--','LineWidth',2); hold on;
    end
end
xlabel('Lat');
ylabel('[%]');
title('((po4 fboff-fbon)/baseline) / ((expC fboff-baseline)/baseline) [%]');
xlim([-80 80]);
ylim([0 28]);
legend('expCmasked','po4tcmasked','expCorig','po4tcorig');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% Regional mean fb strength: (fboff-fbon)/(fboff-baseline)
barwidth=1;
offset = 0.125;
ebwidth = 100; % 1/ebwidth is the width, technically
subplot(428);
fbst_po4tcregsalt = nan(numcircfactors-1,length(fbonrunidx),numregions+1); 
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbst_po4tcregsalt(icircstore,irunstore,:)=squeeze(100*(roundn(po4tcregsalt(icirc,baserunidx,:)-po4tcregsalt(icirc,ifbonrun,:),roundplp)./po4tcregsalt(circ1idx,baserunidx,:))./ ...
            (roundn(expCregsalt(icirc,baserunidx,:)-expCregsalt(circ1idx,baserunidx,:),roundple)./expCregsalt(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbst_po4tcregsalt = squeeze(fbst_po4tcregsalt);
b=bar(1:numregions+1,[fbst_expCregsalt fbst_po4tcregsalt],barwidth);hold on;
set(b(1),'facecolor',expCcol);
set(b(2),'facecolor',po4col);
fbsterrorsalt = abs(bsxfun(@minus,fbst_expCregsalt',fbst_expCregsstdalt));
eb2 = errorbar((1:numregions+1)-offset,fbst_expCregsalt',fbsterrorsalt(1,:)',fbsterrorsalt(2,:)','linestyle','none','color','k','linewidth',2); % errorbar(X,Y,L,U)
errorbar_tick(eb2,ebwidth);
set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
ylabel('[%]');
ylim([0 28]);
xlim([0.5 numregions+1.5]);
title('expC (fboff-fbon)/(fboff-baseline), pts masked out [%]');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

end % end runsection

%--------------------------------------------------
% Figuring out noisy points #2
%--------------------------------------------------
relcircplusrelpo4tcfboff = 100*(deltacircovcirc+roundn(po4tc(:,:,1,baserunidx)-po4tc(:,:,circ1idx,baserunidx),roundplp)./po4tc(:,:,circ1idx,baserunidx));

runsection=0;

if runsection==1

if strcmp(betatype,'Guidi')
    threshold1=2;
    slim=-40;
    nlim=55;
elseif strcmp(betatype,'Kost')
    threshold1=0.4;
    slim=-80;
    nlim=55;
end

relcircplusrelpo4tcfboffbi = nan(size(relcircplusrelpo4tcfboff));
relcircplusrelpo4tcfboffbi(abs(relcircplusrelpo4tcfboff)<threshold1)=1;
relcircplusrelpo4tcfboffbi(abs(relcircplusrelpo4tcfboff)>=threshold1)=0;
%expCfboffreldiffbi(abs(expCfboffreldiff)<threshold&(gridd.YT>nlim|gridd.YT<slim))=1;
%expCfboffreldiffbi(abs(expCfboffreldiff)>=threshold)=0;

% po4tcfboffonreldiff already defined above
% po4tcfboffonreldiff = 100*roundn(po4tc(:,:,1,baserunidx)-po4tc(:,:,1,2),-3)...
%     ./po4tc(:,:,circ1idx,baserunidx);

po4fbst1=100*po4tcfboffonreldiff./relcircplusrelpo4tcfboff;
po4fbstbi1=po4fbst1;
po4fbstbi1(abs(relcircplusrelpo4tcfboff)<threshold1)=nan;
%po4fbstbi1(abs(expCfboffreldiff)<threshold&(gridd.YT>nlim|gridd.YT<slim))=nan;

rdylbu = flipud(cbrewer('div','RdYlBu',100,'linear'));

% - Plot relcircplusrelpo4tcfboff, po4tcfboffonreldiff,
% their ratio, where points are masked out,
% the new zonal mean values excluding the masked out points (with the old zonal mean values),
% and the new regional mean values excluding the masked out points
fignum=fignum+1;
f=figure(fignum);
set(f,'color','white','units','inches','position',[1 1 10 9.5]);
colormap(rdylbu);
subplot(421);
pcolor(gridd.xt,gridd.yt,relcircplusrelpo4tcfboff);
caxis([-20 20]);
shading flat; colorbar;
title('0.9*circ fb off deltacirc/circ + deltaP/P [%]');
subplot(422);
pcolor(gridd.xt,gridd.yt,po4tcfboffonreldiff);
caxis([-5 5]);
shading flat; colorbar;
title('po4tc (0.9*circ fboff - fbon)/(1*circ) [%]');
subplot(423);
pcolor(gridd.xt,gridd.yt,relcircplusrelpo4tcfboffbi);
caxis([-1 1]);
shading flat; colorbar;
title(['1 where abs(0.9*circ fb off deltacirc/circ + deltaP/P)<' num2str(threshold1) '%']);
subplot(424);
pcolor(gridd.xt,gridd.yt,po4fbst1);
caxis([-500 500]);
shading flat; colorbar;
title('[(po4 fboff-fbon)/baseline]/[fboff deltaw/w + deltaP/P] [%]');
subplot(425);
pcolor(gridd.xt,gridd.yt,po4fbstbi1);
caxis([-500 500]);
shading flat; colorbar;
title('[(po4 fboff-fbon)/baseline]/[fboff deltaw/w + deltaP/P] [%]');

%disp('pause')
%pause

%----------- BEGIN RECALCULATIONS #2
% Need recalculating of:
% po4tcz, po4tcregs, expCz, expCregs, fbst_expCregs, fbst_expCregsstd 

% - Recalculate po4tcz, po4tcregs
po4tczalt1 = nan(size(po4tcz));
po4tcregsalt1 = nan(size(po4tcregs));
totpo4tcalt1 = nan(size(totpo4tc));
for icirc=1:numcircfactors
    for irun=1:numruntypes
        po4tctemp = po4tc(:,:,icirc,irun);
        po4tctemp(abs(relcircplusrelpo4tcfboff)<threshold1) = nan;
        totpo4tcalt1(icirc,irun) = nansum(nansum(po4tctemp.*Areatocnonly));
        po4tczalt1(:,icirc,irun) = nansum(po4tctemp.*Areatocnonly,2)./Areatzonalsum;
        for iregion=1:4
            po4tcregsalt1(icirc,irun,iregion) = nansum(po4tctemp(regmappgrid==iregion).*Areatocnonly(regmappgrid==iregion))./nansum(Areatocnonly(regmappgrid==iregion)); % regional po4tcs
        end
        % Reorder Indian Ocean
        for iregion=5
            po4tcregsalt1(icirc,irun,iregion) = nansum(po4tctemp(regmappgrid==9).*Areatocnonly(regmappgrid==9))./nansum(Areatocnonly(regmappgrid==9)); % regional po4tcs
        end
        for iregion=6:9
            po4tcregsalt1(icirc,irun,iregion) = nansum(po4tctemp(regmappgrid==iregion-1).*Areatocnonly(regmappgrid==iregion-1))./nansum(Areatocnonly(regmappgrid==iregion-1)); % regional po4tcs
        end
        po4tcregsalt1(icirc,irun,end) = totpo4tcalt1(icirc,irun)./totocnarea; % global po4tc
    end
end

% - Recalculate expCz, expCregs
expCzalt1 = nan(size(expCz));
expCregsalt1 = nan(size(expCregs));
totexpCalt1 = nan(size(totexpC));
for icirc=1:numcircfactors
    for irun=1:numruntypes
        expCtemp = expC(:,:,icirc,irun);
        expCtemp(abs(relcircplusrelpo4tcfboff)<threshold1) = nan;
        totexpCalt1(icirc,irun) = nansum(nansum(expCtemp.*Areatocnonly));
        expCzalt1(:,icirc,irun) = nansum(expCtemp.*Areatocnonly,2)./Areatzonalsum;
        for iregion=1:4
            expCregsalt1(icirc,irun,iregion) = nansum(expCtemp(regmappgrid==iregion).*Areatocnonly(regmappgrid==iregion))./nansum(Areatocnonly(regmappgrid==iregion)); % regional exports
        end
        % Reorder Indian Ocean
        for iregion=5
            expCregsalt1(icirc,irun,iregion) = nansum(expCtemp(regmappgrid==9).*Areatocnonly(regmappgrid==9))./nansum(Areatocnonly(regmappgrid==9)); % regional exports
        end
        for iregion=6:9
            expCregsalt1(icirc,irun,iregion) = nansum(expCtemp(regmappgrid==iregion-1).*Areatocnonly(regmappgrid==iregion-1))./nansum(Areatocnonly(regmappgrid==iregion-1)); % regional exports
        end
        expCregsalt1(icirc,irun,end) = totexpCalt1(icirc,irun)./totocnarea; % global export
    end
end

% - Recalculate fbst_expCregs, fbst_expCregsstd 
expCzstdalt1 = nan(size(expCzstd));
totexpCstdalt1 = nan(size(totexpCstd));
expCregsstdalt1 = nan(size(expCregsstd));
for icircstore=1:numcircfactors-1
    for irun=1:2 % 1-stdev above and below for each beta
        expCstdtemp = expCstd(:,:,icircstore,irun);
        expCstdtemp(abs(relcircplusrelpo4tcfboff)<threshold1) = nan;
        expCzstdalt1(:,icircstore,irun) = nansum(expCstdtemp.*Areatocnonly,2)./Areatzonalsum;
        totexpCstdalt1(icircstore,irun) = nansum(nansum(expCstdtemp.*Areatocnonly));
        for iregion=1:4
            expCregsstdalt1(icircstore,irun,iregion) = nansum(expCstdtemp(regmappgrid==iregion).*Areatocnonly(regmappgrid==iregion))./nansum(Areatocnonly(regmappgrid==iregion)); % regional exports
        end
        % Reorder Indian Ocean
        for iregion=5
            expCregsstdalt1(icircstore,irun,iregion) = nansum(expCstdtemp(regmappgrid==9).*Areatocnonly(regmappgrid==9))./nansum(Areatocnonly(regmappgrid==9)); % regional exports
        end
        for iregion=6:9
            expCregsstdalt1(icircstore,irun,iregion) = nansum(expCstdtemp(regmappgrid==iregion-1).*Areatocnonly(regmappgrid==iregion-1))./nansum(Areatocnonly(regmappgrid==iregion-1)); % regional exports
        end
        expCregsstdalt1(icircstore,irun,end) = totexpCstdalt1(icircstore,irun)./totocnarea; % global export
    end
end

fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
fbst_expCregsalt1 = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbst_expCregsalt1(icircstore,irunstore,:)=squeeze(100*roundn(expCregsalt1(icirc,baserunidx,:)-expCregsalt1(icirc,ifbonrun,:),roundple)./(expCregsalt1(icirc,baserunidx,:)-expCregsalt1(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbst_expCregsalt1 = squeeze(fbst_expCregsalt1);

fbst_expCregsstdalt1 = nan(numcircfactors-1,2,numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for irun = 1:2
        fbst_expCregsstdalt1(icircstore,irun,:)=squeeze(100*roundn(expCregsalt1(icirc,baserunidx,:)-expCregsstdalt1(icircstore,irun,:),roundple)./(expCregsalt1(icirc,baserunidx,:)-expCregsalt1(circ1idx,baserunidx,:)));
    end
    icircstore = icircstore+1;
end
fbst_expCregsstdalt1 = squeeze(fbst_expCregsstdalt1);
%----------- FINISH RECALCULATIONS #2

% - Continue plotting expCfboffreldiff, po4tcfboffonreldiff,
% their ratio, where points are masked out,
% the new zonal mean values excluding the masked out points (with the old zonal mean values),
% and the new regional mean values excluding the masked out points

% Zonal mean fb strength: (fboff-fbon)/(fboff-baseline)
%expCcol = 'r';
expCcol = [0 0 0.8];
%po4col = 'm';
po4col = [1 0.4 0];
fontwt = 'normal';
subplot(427);
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCzalt1(:,icirc,baserunidx)-expCzalt1(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(expCzalt1(:,icirc,baserunidx)-expCzalt1(:,circ1idx,baserunidx)),...
            'Color',expCcol,'LineStyle','-','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            (roundn(po4tczalt1(:,icirc,baserunidx)-po4tczalt1(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(po4tczalt1(:,circ1idx,baserunidx)))./ ...
            (deltacircovcirc+roundn(po4tczalt1(:,icirc,baserunidx)-po4tczalt1(:,circ1idx,baserunidx),roundplp)...
            ./po4tczalt1(:,circ1idx,baserunidx)), ...
            po4col,'LineStyle','-','LineWidth',2); hold on;
    end
end
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCz(:,icirc,baserunidx)-expCz(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(expCz(:,icirc,baserunidx)-expCz(:,circ1idx,baserunidx)),...
            'Color',expCcol,'LineStyle','--','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            (roundn(po4tcz(:,icirc,baserunidx)-po4tcz(:,icirc,fbonrunidx(ifbonrun)),roundplp)...
            ./(po4tcz(:,circ1idx,baserunidx)))./ ...
            (deltacircovcirc+roundn(po4tcz(:,icirc,baserunidx)-po4tcz(:,circ1idx,baserunidx),roundplp)...
            ./po4tcz(:,circ1idx,baserunidx)), ...
            po4col,'LineStyle','--','LineWidth',2); hold on;
    end
end
xlabel('Lat');
ylabel('[%]');
title('[(po4 fboff-fbon)/baseline]/[fboff deltaw/w + deltaP/P] [%]');
xlim([-80 80]);
ylim([0 28]);
legend('expCmasked','po4tcmasked','expCorig','po4tcorig');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% Regional mean fb strength: (fboff-fbon)/(fboff-baseline)
barwidth=1;
offset = 0.125;
ebwidth = 100; % 1/ebwidth is the width, technically
subplot(428);
fbst_po4tcregsalt1 = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbst_po4tcregsalt1(icircstore,irunstore,:)=squeeze(100*(roundn(po4tcregsalt1(icirc,baserunidx,:)-po4tcregsalt1(icirc,ifbonrun,:),roundplp)./po4tcregsalt1(circ1idx,baserunidx,:))./ ...
            (deltacircovcirc+roundn(po4tcregsalt1(icirc,baserunidx,:)-po4tcregsalt1(circ1idx,baserunidx,:),roundplp)./po4tcregsalt1(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbst_po4tcregsalt1 = squeeze(fbst_po4tcregsalt1);
b=bar(1:numregions+1,[fbst_expCregsalt1 fbst_po4tcregsalt1],barwidth);hold on;
set(b(1),'facecolor',expCcol);
set(b(2),'facecolor',po4col);
fbsterrorsalt1 = abs(bsxfun(@minus,fbst_expCregsalt1',fbst_expCregsstdalt1));
eb2 = errorbar((1:numregions+1)-offset,fbst_expCregsalt1',fbsterrorsalt1(1,:)',fbsterrorsalt1(2,:)','linestyle','none','color','k','linewidth',2); % errorbar(X,Y,L,U)
errorbar_tick(eb2,ebwidth);
set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
ylabel('[%]');
ylim([0 28]);
xlim([0.5 numregions+1.5]);
title('expC (fboff-fbon)/(fboff-baseline), pts masked out [%]');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

end % end runsection

%---------------------------------
% Where does deltaw/w+deltaP/P not explain deltaE/E?
%---------------------------------
plotsection=1;

if plotsection==1

expCfbonreldiff=100*roundn(expC(:,:,1,fbonrunidx)-expC(:,:,circ1idx,baserunidx),roundple)./expC(:,:,circ1idx,baserunidx);
relcircplusrelpo4tcfbon = 100*(deltacircovcirc+roundn(po4tc(:,:,1,fbonrunidx)-po4tc(:,:,circ1idx,baserunidx),roundplp)./po4tc(:,:,circ1idx,baserunidx));

fignum=fignum+1;
f=figure(fignum);
set(f,'color','white','units','inches','position',[1 1 12.5 8]);
rdylbu = flipud(cbrewer('div','RdYlBu',100,'linear'));
colormap(rdylbu);

subplot(331);
pcolor(gridd.xt,gridd.yt,expCfboffreldiff-relcircplusrelpo4tcfboff);
shading flat;colorbar;
caxis([-20 20]);
title('0.9circ fboff deltaE/E - (deltaw/w+deltaP/P) [%]');
subplot(332);
pcolor(gridd.xt,gridd.yt,expCfboffreldiff);
shading flat;colorbar;
caxis([-20 20]);
title('0.9circ fboff deltaE/E [%]');
subplot(333);
pcolor(gridd.xt,gridd.yt,relcircplusrelpo4tcfboff);
shading flat;colorbar;
caxis([-20 20]);
title('0.9circ fboff (deltaw/w+deltaP/P) [%]');

subplot(334);
pcolor(gridd.xt,gridd.yt,expCfbonreldiff-relcircplusrelpo4tcfbon);
shading flat;colorbar;
caxis([-20 20]);
title('0.9circ fbon deltaE/E - (deltaw/w+deltaP/P) [%]');
subplot(335);
pcolor(gridd.xt,gridd.yt,expCfbonreldiff);
shading flat;colorbar;
caxis([-20 20]);
title('0.9circ fbon deltaE/E [%]');
subplot(336);
pcolor(gridd.xt,gridd.yt,relcircplusrelpo4tcfbon);
shading flat;colorbar;
caxis([-20 20]);
title('0.9circ fbon (deltaw/w+deltaP/P) [%]');

subplot(337);
pcolor(gridd.xt,gridd.yt,100*[(po4tcfboffonreldiff./expCfboffreldiff)-(po4tcfboffonreldiff./relcircplusrelpo4tcfboff)]);
shading flat;colorbar;
caxis([-20 20]);
title('[(po4 fboff-fbon)/baseline]/[fboff deltaE/E] - []/[deltaw/w + deltaP/P] [%]');
subplot(338);
pcolor(gridd.xt,gridd.yt,100*po4tcfboffonreldiff./expCfboffreldiff);
shading flat;colorbar;
caxis([-40 40]);
title('[(po4 fboff-fbon)/baseline]/[fboff deltaE/E] [%]');
subplot(339);
pcolor(gridd.xt,gridd.yt,100*po4tcfboffonreldiff./relcircplusrelpo4tcfboff);
shading flat;colorbar;
caxis([-40 40]);
title('[(po4 fboff-fbon)/baseline]/[fboff deltaw/w + deltaP/P] [%]');

expCfbst = 100*roundn(expC(:,:,1,baserunidx)-expC(:,:,1,fbonrunidx),roundple)./(expC(:,:,1,baserunidx)-expC(:,:,circ1idx,baserunidx));
fignum=fignum+1;
f=figure(fignum);
set(f,'color','white','units','inches','position',[1 1 12.5 8]);
colormap(rdylbu);
subplot(331);
pcolor(expCfbst);
shading flat;colorbar;
caxis([-40 40]);
title(['(fboff-fbon)/(fboff-baseline)' 10 '[%]']);
subplot(332);
pcolor(expCfbst-100*po4tcfboffonreldiff./expCfboffreldiff);
shading flat;colorbar;
caxis([-20 20]);
title(['(fboff-fbon)/(fboff-baseline)-' 10 '[(po4 fboff-fbon)/baseline]/[fboff deltaE/E] [%]']);
subplot(333);
pcolor(expCfbst-100*po4tcfboffonreldiff./relcircplusrelpo4tcfboff);
shading flat;colorbar;
caxis([-20 20]);
title(['(fboff-fbon)/(fboff-baseline)-' 10 '[(po4 fboff-fbon)/baseline]/[fboff deltaw/w + deltaP/P] [%]']);

end % end plotsection

%---------------------------------
% FIGURE 7 REMADE #1 (same as FIGURE 7 above but with deltaexpCoff threshold points removed,
% though only for po4 related variables, not expC related variables):
% Row 1: 
% Row 2: 
% Row 3: 
%---------------------------------
plotsection=0;

if plotsection==1

fignum=fignum+1;
f=figure(fignum); set(f,'Color','White','Units','Inches','Position',[1 1 10 9.5]);
fbofflinestyle='-';
fbonlinestyle='--';
fbofflinewidth=2;
fbonlinewidth=2;
%expCcol = 'r';
expCcol = [0 0 0.8];
%po4col = 'm';
po4col = [1 0.4 0];
fontwt = 'normal';

% - Zonal mean rel change diff: (fboff-fbon)/baseline
subplot(421);
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCz(:,icirc,baserunidx)-expCz(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(expCz(:,circ1idx,baserunidx)),...
            'Color',expCcol,'LineStyle','-','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(po4tczalt(:,icirc,baserunidx)-po4tczalt(:,icirc,fbonrunidx(ifbonrun)),roundplp)...
            ./(po4tczalt(:,circ1idx,baserunidx)),...
            po4col,'LineStyle','-','LineWidth',2); hold on;
    end
end
xlabel('Lat');
ylabel('[%]');
title('(fboff-fbon)/baseline [%]');
xlim([-80 80]);
ylim([-3.5 0]);
legend('expC','po4tc');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Regional mean rel change diff: (fboff-fbon)/baseline
subplot(422);
barwidth = 1;
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
fbonminusoffrel_po4tcregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbonminusoffrel_po4tcregs(icircstore,irunstore,:)=squeeze(100*roundn(po4tcregsalt(icirc,baserunidx,:)-po4tcregsalt(icirc,ifbonrun,:),roundplp)./(po4tcregsalt(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbonminusoffrel_po4tcregs = squeeze(fbonminusoffrel_po4tcregs);
fbonminusoffrel_expCregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbonminusoffrel_expCregs(icircstore,irunstore,:)=squeeze(100*roundn(expCregs(icirc,baserunidx,:)-expCregs(icirc,ifbonrun,:),roundple)./(expCregs(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbonminusoffrel_expCregs = squeeze(fbonminusoffrel_expCregs);
b=bar(1:numregions+1,[fbonminusoffrel_expCregs fbonminusoffrel_po4tcregs],barwidth);hold on;
set(b(1),'facecolor',expCcol);
set(b(2),'facecolor',po4col);
set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
ylabel('[%]');
ylim([-3.5 0]);
xlim([0.5 numregions+1.5]);
title('(fboff-fbon)/(baseline) [%]');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Zonal mean fb strength: (fboff-fbon)/(fboff-baseline)
subplot(423);
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCz(:,icirc,baserunidx)-expCz(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(expCz(:,icirc,baserunidx)-expCz(:,circ1idx,baserunidx)),...
            'Color',expCcol,'LineStyle','-','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat', 100*...
            (roundn(po4tczalt(:,icirc,baserunidx)-po4tczalt(:,icirc,fbonrunidx(ifbonrun)),roundplp)...
            ./(po4tczalt(:,circ1idx,baserunidx)))./ ...
            (roundn(expCzalt(:,icirc,baserunidx)-expCzalt(:,circ1idx,baserunidx),roundple)...
            ./expCzalt(:,circ1idx,baserunidx)), ...
            po4col,'LineStyle','-','LineWidth',2); hold on;
    end
end
xlabel('Lat');
ylabel('[%]');
title('((po4 fboff-fbon)/baseline) / ((expC fboff-baseline)/baseline) [%]');
xlim([-80 80]);
ylim([0 28]);
%legend('expC','po4tc');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Regional mean fb strength: (fboff-fbon)/(fboff-baseline)
barwidth=1;
offset = 0.125;
ebwidth = 100; % 1/ebwidth is the width, technically
subplot(424);
fbst_po4tcregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1); 
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbst_po4tcregs(icircstore,irunstore,:)=squeeze(100*(roundn(po4tcregsalt(icirc,baserunidx,:)-po4tcregsalt(icirc,ifbonrun,:),roundplp)./po4tcregsalt(circ1idx,baserunidx,:))./ ...
            (roundn(expCregsalt(icirc,baserunidx,:)-expCregsalt(circ1idx,baserunidx,:),roundple)./expCregsalt(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbst_po4tcregs = squeeze(fbst_po4tcregs);
b=bar(1:numregions+1,[fbst_expCregs fbst_po4tcregs],barwidth);hold on;
set(b(1),'facecolor',expCcol);
set(b(2),'facecolor',po4col);
fbsterrors = abs(bsxfun(@minus,fbst_expCregs',fbst_expCregsstd));
eb2 = errorbar((1:numregions+1)-offset,fbst_expCregs',fbsterrors(1,:)',fbsterrors(2,:)','linestyle','none','color','k','linewidth',2); % errorbar(X,Y,L,U)
errorbar_tick(eb2,ebwidth);
set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
ylabel('[%]');
ylim([0 28]);
xlim([0.5 numregions+1.5]);
title('expC (fboff-fbon)/(fboff-baseline) [%]');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

end % end plotsection

%---------------------------------
% FIGURE 7 REMADE #2 (same as FIGURE 7 above but with (rel po4 on - off)/(deltaw/w + deltap/p)
% threshold points removed,
% though only for po4 related variables, not expC related variables):
% Row 1: 
% Row 2: 
% Row 3: 
%---------------------------------
plotsection=0;

if plotsection==1

fignum=fignum+1;
f=figure(fignum); set(f,'Color','White','Units','Inches','Position',[1 1 10 9.5]);
fbofflinestyle='-';
fbonlinestyle='--';
fbofflinewidth=2;
fbonlinewidth=2;
%expCcol = 'r';
expCcol = [0 0 0.8];
%po4col = 'm';
po4col = [1 0.4 0];
fontwt = 'normal';

% - Zonal mean rel change diff: (fboff-fbon)/baseline
subplot(421);
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCz(:,icirc,baserunidx)-expCz(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(expCz(:,circ1idx,baserunidx)),...
            'Color',expCcol,'LineStyle','-','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(po4tczalt1(:,icirc,baserunidx)-po4tczalt1(:,icirc,fbonrunidx(ifbonrun)),roundplp)...
            ./(po4tczalt1(:,circ1idx,baserunidx)),...
            po4col,'LineStyle','-','LineWidth',2); hold on;
    end
end
xlabel('Lat');
ylabel('[%]');
title('(fboff-fbon)/baseline [%]');
xlim([-80 80]);
ylim([-3.5 0]);
legend('expC','po4tc');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Regional mean rel change diff: (fboff-fbon)/baseline
subplot(422);
barwidth = 1;
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
fbonminusoffrel_po4tcregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbonminusoffrel_po4tcregs(icircstore,irunstore,:)=squeeze(100*roundn(po4tcregsalt1(icirc,baserunidx,:)-po4tcregsalt1(icirc,ifbonrun,:),roundplp)./(po4tcregsalt1(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbonminusoffrel_po4tcregs = squeeze(fbonminusoffrel_po4tcregs);
fbonminusoffrel_expCregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1);
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbonminusoffrel_expCregs(icircstore,irunstore,:)=squeeze(100*roundn(expCregs(icirc,baserunidx,:)-expCregs(icirc,ifbonrun,:),roundple)./(expCregs(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbonminusoffrel_expCregs = squeeze(fbonminusoffrel_expCregs);
b=bar(1:numregions+1,[fbonminusoffrel_expCregs fbonminusoffrel_po4tcregs],barwidth);hold on;
set(b(1),'facecolor',expCcol);
set(b(2),'facecolor',po4col);
set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
ylabel('[%]');
ylim([-3.5 0]);
xlim([0.5 numregions+1.5]);
title('(fboff-fbon)/(baseline) [%]');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Zonal mean fb strength: (fboff-fbon)/(fboff-baseline)
subplot(423);
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat',100*...
            roundn(expCz(:,icirc,baserunidx)-expCz(:,icirc,fbonrunidx(ifbonrun)),roundple)...
            ./(expCz(:,icirc,baserunidx)-expCz(:,circ1idx,baserunidx)),...
            'Color',expCcol,'LineStyle','-','LineWidth',2); hold on;
    end
end
fbonrunidx=find(strcmp(fbstatus,'Pfbon'));
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    for ifbonrun = 1:length(fbonrunidx)
        plot(lat', 100*...
            (roundn(po4tczalt1(:,icirc,baserunidx)-po4tczalt1(:,icirc,fbonrunidx(ifbonrun)),roundplp)...
            ./(po4tczalt1(:,circ1idx,baserunidx)))./ ...
            (deltacircovcirc+roundn(po4tczalt1(:,icirc,baserunidx)-po4tczalt1(:,circ1idx,baserunidx),roundplp)...
            ./po4tczalt1(:,circ1idx,baserunidx)), ...
            po4col,'LineStyle','-','LineWidth',2); hold on;
    end
end
xlabel('Lat');
ylabel('[%]');
title('[(po4 fboff-fbon)/baseline]/[fboff deltaw/w + deltaP/P] [%]');
xlim([-80 80]);
ylim([0 30]);
%legend('expC','po4tc');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

% - Regional mean fb strength: (fboff-fbon)/(fboff-baseline)
barwidth=1;
offset = 0.125;
ebwidth = 100; % 1/ebwidth is the width, technically
subplot(424);
fbst_po4tcregs = nan(numcircfactors-1,length(fbonrunidx),numregions+1); 
icircstore = 1;
for icirc=[1:circ1idx-1 circ1idx+1:numcircfactors]
    irunstore = 1;
    for ifbonrun = fbonrunidx
        fbst_po4tcregs(icircstore,irunstore,:)=squeeze(100*(roundn(po4tcregsalt1(icirc,baserunidx,:)-po4tcregsalt1(icirc,ifbonrun,:),roundplp)./po4tcregsalt1(circ1idx,baserunidx,:))./ ...
            (deltacircovcirc+roundn(po4tcregsalt1(icirc,baserunidx,:)-po4tcregsalt1(circ1idx,baserunidx,:),roundplp)./po4tcregsalt1(circ1idx,baserunidx,:)));
        irunstore = irunstore+1;
    end
    icircstore = icircstore+1;
end
fbst_po4tcregs = squeeze(fbst_po4tcregs);
b=bar(1:numregions+1,[fbst_expCregs fbst_po4tcregs],barwidth);hold on;
set(b(1),'facecolor',expCcol);
set(b(2),'facecolor',po4col);
fbsterrors = abs(bsxfun(@minus,fbst_expCregs',fbst_expCregsstd));
eb2 = errorbar((1:numregions+1)-offset,fbst_expCregs',fbsterrors(1,:)',fbsterrors(2,:)','linestyle','none','color','k','linewidth',2); % errorbar(X,Y,L,U)
errorbar_tick(eb2,ebwidth);
set(gca, 'XTick', 1:numregions+1, 'XTickLabel', [regnamesabbrev {'glob'}]);
ylabel('[%]');
ylim([0 30]);
xlim([0.5 numregions+1.5]);
title('expC (fboff-fbon)/(fboff-baseline) [%]');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight',fontwt,'TickDir','out');

end % end plotsection
