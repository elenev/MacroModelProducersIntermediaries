if ~exist('batchmode','var')
	if usejava('desktop')
		batchmode=true;
	else
		batchmode=false;
	end
end
if ~batchmode
   clear; 
end
close all;

%--------------------------------------------------------------------------
% simulation setup
%------------------------------------ --------------------------------------

% file with model
respath='./';
if ~exist('resfile','var') 
    disp('resfile not defined. Opening default file instead.');
    resfile='res_20200312_benchriskysafesigI0_s130'; 
%    resfile='res_20190808_bench_s19c46'; 
end
load([respath,resfile,'.mat']);

if exist('distmat','var')
	computeVFconvergence = true;
else
	computeVFconvergence = false;
end

% Update params
augmentParams=false;
expdef='experdef_20190830.m';
if augmentParams
    run(expdef);
    mobj=mobj.augmentParams(allexpers.bench.params);
end


% number of periods and burn-in
NT_sim=10000;
NT_ini=100;

% compute Euler equation error?
compEEErr=1;
% make graphs for bankruptcy episodes
makeBRGraphs=0;

% Winsorize
winsorize=0;
winsfile='';
cutoff=99.5;


% Force the creation of a sim_res file
force_output = 1;

% output table file
outpath='./Results/';
outstats_exog=[outpath,'statsexog_',resfile,'.xls'];
outstats_exog_nolowWI=[outpath,'statsexog_nolowWI_',resfile,'.xls'];
outstats_endog=[outpath,'statsendog_',resfile,'.xls'];
errstats=[outpath,'errstats_',resfile,'.xls'];
expdata=0;
outdata=[outpath,'series_',resfile,'.csv'];
       
%--------------------------------------------------------------------------
% start simulation
%--------------------------------------------------------------------------

% CBS model?
if ~isfield(mobj.Params,'CBS')
    CBS=false;
else
    if ~mobj.Params.CBS
        CBS=false;
    else
        CBS=true;
    end
end

% set starting point
start_ex=5;
startpt=struct;
startpt.KB=stv{1}.State.KB;
startpt.LB=stv{1}.State.LB;
if ~CBS
    startpt.WI=stv{1}.State.WI;
end
startpt.BG=stv{1}.State.BG;
startpt=orderfields(startpt,mobj.En_names);
startpt_vec=model.DSGEModel.structToVec(startpt)';
startpt_vec=[start_ex,startpt_vec];

% simulate
[simseries,varnames,errmat,Wshtrans,SDFmat]=mobj.simulate(NT_sim+1,NT_ini,startpt_vec,compEEErr);
simseries_orig=simseries;
varnames_orig=varnames;
statevec = simseries(:,1);

% Saver SDFs
exnpt=mobj.Exogenv.exnpt;
if ~CBS
    SDFRmat = SDFmat(:,1:exnpt);
    [simseries, varnames] = mobj.computeSimulationMoments(simseries,varnames,SDFRmat);
else
    [simseries, varnames] = mobj.computeSimulationMoments(simseries,varnames);
end
    
nvars = length(varnames);

% Create table object for easier access
simtable=array2table(simseries);
[~,ia,~]=unique(varnames);
simtable=simtable(:,ia);
simtable.Properties.VariableNames=varnames(ia);
dispnames=varnames(ia);

% make HashMap with mapping of names to indices
indexmap=java.util.HashMap;
for i=1:nvars
    indexmap.put(varnames{i},i);
end

% Check transition function errors
idx = sub2ind([NT_sim+1,mobj.Exogenv.exnpt],(1:NT_sim+1)',[statevec(2:end);1]);
idx=idx(1:end-1);

KBtrans=Wshtrans(:,1:mobj.Exogenv.exnpt);
KB_err=simseries(:,indexmap.get('KB')) - KBtrans(idx); % column index: 6
errmat = [errmat, [KB_err;0]];

LBtrans=Wshtrans(:,mobj.Exogenv.exnpt+1:end);
LB_err=simseries(:,indexmap.get('LB')) - LBtrans(idx); % column index: 6
errmat = [errmat, [LB_err;0]];

if ~CBS
    WItrans=Wshtrans(:,2*mobj.Exogenv.exnpt+1:end);
    WI_err=simseries(:,indexmap.get('WI')) - WItrans(idx);
    errmat = [errmat, [WI_err;0]];
    BGadd=1;
else
    BGadd=0;
end

BGtrans=Wshtrans(:,(2+BGadd)*mobj.Exogenv.exnpt+1:end);
BG_err=simseries(:,indexmap.get('BG')) - BGtrans(idx);
errmat = [errmat, [BG_err;0]];

if ~CBS
    state_range=6:9;
else
    state_range=6:8;
end
% SDFR_rlz = SDFRmat(idx);
% SDFS_rlz = SDFSmat(idx);

% Check value function convergence
if computeVFconvergence
	Dfct=grid.LinearInterpFunction(mobj.Vfct.SSGrid,distmat);
	simdist=Dfct.evaluateAt([statevec(2:end),simseries(:, state_range)])';
	errmat = [errmat, [zeros(1,size(simdist,2));simdist] ];
end

% Winsorize

if winsorize
    if isempty(winsfile)
        prct=prctile(abs(errmat(2:end,:)),cutoff);
        bad_idx = any(abs(errmat(2:end,:)) > repmat(prct,NT_sim-1,1) , 2);
    else
       load([respath,winsfile,'.mat'],'bad_idx');
    end
    simseries = simseries(~bad_idx,:);
    errmat = errmat([true;~bad_idx],:);
    NT_sim = size(simseries,1) + 1;
    save([respath,resfile,'.mat'],'bad_idx','-append');
end

%--------------------------------------------------------------------------
% calculate stats
%--------------------------------------------------------------------------
varst=zeros(length(startpt_vec)-1,1);
for i=1:length(startpt_vec)-1
    varst(i)=indexmap.get(mobj.En_names{i});
end
            
% state variable means in stationary distribution
stvstat=mean(simseries(:,varst));

% calculate business cycle stats
% first for all periods, then separately for low and high sigma_omega states
% 1. condition on exogenous states (smpsel_exog)
% 2. condition on endogenous states: GDP growth, default rate (smpsel_endog)
statsout_exog=cell(4,1);
statsout_endog=cell(4,1);

% first for all periods, then separately for low and high sigma_omega states
smpsel_exog={true(NT_sim-1,1), simseries(:,2)==mobj.Params.sig2_om(1) & simseries(:,1) >= 1+mobj.Params.zA, ...
                          simseries(:,2)==mobj.Params.sig2_om(1)  & simseries(:,1) < 1+mobj.Params.zA, ...
                          simseries(:,2)==mobj.Params.sig2_om(2)  & simseries(:,1) < 1+mobj.Params.zA, ...
                          simseries(:,2)==mobj.Params.sig2_om(1)};

% condition on endogenous variables
% thresholds
Y_gr_thr = 0;
Drate_thr = 0;
smpsel_endog={true(NT_sim-1,1), simseries(:,indexmap.get('Y_gr')) >= Y_gr_thr & simseries(:,indexmap.get('Drate')) < Drate_thr, ...
                          simseries(:,indexmap.get('Y_gr')) < Y_gr_thr  & simseries(:,indexmap.get('Drate')) < Drate_thr, ...
                          simseries(:,indexmap.get('Y_gr')) < Y_gr_thr  & simseries(:,indexmap.get('Drate')) >= Drate_thr};

gdp_idx = indexmap.get('Y');                     
gdpgr_idx = indexmap.get('Y_gr'); 
for j=1:numel(smpsel_exog)
    % subsample
    simtmp=real(simseries(smpsel_exog{j},:));
    statstmp=zeros(nvars,17);
    statstmp(:,1)=nanmean(simtmp)';
    statstmp(:,2)=nanstd(simtmp)';
    % contemp and first-order autocorrelations
    autocorrm=corrcoef([simtmp(2:end,:),simtmp(1:end-1,:)]);
    conm=autocorrm(1:nvars,1:nvars);
    lagm=autocorrm(nvars+1:end,1:nvars);
    % corr with shocks
    statstmp(:,3:4)=[conm(:,1),lagm(1,:)'];
    statstmp(:,5:6)=[conm(:,2),lagm(2,:)'];
    % corr with Y_gr
    statstmp(:,7:8)=[conm(:,gdp_idx),lagm(gdp_idx,:)'];
    statstmp(:,9:10)=[conm(:,gdpgr_idx),lagm(gdpgr_idx,:)'];
    % vector with fo autocorr
    statstmp(:,11)=diag(lagm);
	statstmp(:,12:22)=prctile(simtmp,[0,1,5,10,25,50,75,90,95,99,100])';
    statsout_exog{j}=statstmp;
    
    if ~CBS
        % percentiles of intermediary wealth depending on subsample
        WI_pct_exog(:,j) = prctile(simseries(:,indexmap.get('WI')),[5 10 20 25 30 35 40 45 50 90 95]);
    end
end


for j=1:numel(smpsel_endog)
    % subsample
    simtmp=simseries(smpsel_endog{j},:);
    statstmp=zeros(nvars,18);
    statstmp(:,1)=nanmean(simtmp)';
    statstmp(:,2)=nanstd(simtmp)';
    % contemp and first-order autocorrelations
    autocorrm=corrcoef([simtmp(2:end,:),simtmp(1:end-1,:)]);
    conm=autocorrm(1:nvars,1:nvars);
    lagm=autocorrm(nvars+1:end,1:nvars);
    % corr with shocks
    statstmp(:,3:4)=[conm(:,1),lagm(1,:)'];
    statstmp(:,5:6)=[conm(:,2),lagm(2,:)'];
    % corr with Y_gr
    statstmp(:,7:8)=[conm(:,gdp_idx),lagm(gdp_idx,:)'];
    statstmp(:,9:10)=[conm(:,gdpgr_idx),lagm(gdpgr_idx,:)'];
    % vector with fo autocorr
    statstmp(:,11)=diag(lagm);
	statstmp(:,12:22)=prctile(simtmp,[0,1,5,10,25,50,75,90,95,99,100])';
    statsout_endog{j}=statstmp;
    
    if ~CBS
        % percentiles of intermediary wealth depending on subsample
        WI_pct_endog(:,j) = prctile(simseries(:,indexmap.get('WI')),[5 10 20 25 30 35 40 45 50 90 95]);
    end
end

if ~CBS
    disp('Percentiles of intermediary wealth (exog subsamples): ');
    fprintf('\n');
    disp('unconditional: '); disp( WI_pct_exog(:,1));
    disp('expansion: '); disp( WI_pct_exog(:,2));
    disp('non-fin recession: '); disp( WI_pct_exog(:,3));
    disp('fin recession: '); disp( WI_pct_exog(:,4));
    
    disp('Percentiles of intermediary wealth (endog subsamples): ');
    fprintf('\n');
    disp('unconditional: '); disp( WI_pct_endog(:,1));
    disp('expansion: '); disp( WI_pct_endog(:,2));
    disp('non-fin recession: '); disp( WI_pct_endog(:,3));
    disp('fin recession: '); disp( WI_pct_endog(:,4));    
end

%--------------------------------------------------------------------------
% output
%--------------------------------------------------------------------------

% overview output for eyeball check against analytic st.st. values
% make one big structure with steady-state values
stvbig=model.HelperCollection.combineStructs({stv{1}.Sol,stv{1}.State,stv{1}.Add,stv{1}.statsout});

% output table
% make index vector
[displist,dispnames]=model.HelperCollection.makeListFromNames(indexmap,dispnames);
ndvars=length(displist);

disp(' ');
disp('Simulation steady state');

% overview output 
fprintf('Frequency (exog subsamples): ');
for j=1:numel(smpsel_exog)
    % select vars
    tabout_exog{j}=statsout_exog{j}(displist,:);
    fprintf('%f\t',sum(smpsel_exog{j}));
end
fprintf('\n');
fprintf('Frequency (endog subsamples): ');
for j=1:numel(smpsel_endog)
    % select vars
    tabout_endog{j}=statsout_endog{j}(displist,:);
    fprintf('%f\t',sum(smpsel_endog{j}));
end
fprintf('\n');
disp('-------------');

for s=1:ndvars
    if isfield(stvbig,dispnames{s})
        ststval=stvbig.(dispnames{s});
    else
        ststval=0;
    end
    if numel(dispnames{s}) > 7
        fprintf('%d\t%4s\t\t\t\t%f |',displist(s),dispnames{s},ststval);
    else
        fprintf('%d\t%4s\t\t\t\t\t%f |',displist(s),dispnames{s},ststval);
    end
%     disp('Exog subsamples')
    for j=1:numel(smpsel_exog)
        fprintf('\t%f, %f |',tabout_exog{j}(s,1),tabout_exog{j}(s,2));
    end
    
%     fprintf('\n');
%     disp('Endog subsamples')
%     for j=1:numel(smpsel_endog)
%         fprintf('\t%f, %f |',tabout_endog{j}(s,1),tabout_endog{j}(s,2));
%     end    
    fprintf('\n');    
end

if compEEErr
    avg_err=mean(abs(errmat))';
    med_err=median(abs(errmat))';
    p75_err=prctile(abs(errmat),75)';
    p95_err=prctile(abs(errmat),95)';
    p99_err=prctile(abs(errmat),99)';
    p995_err=prctile(abs(errmat),99.5)';
    max_err=max(abs(errmat))';
    errtab=table(avg_err,med_err,p75_err,p95_err,p99_err,p995_err,max_err);
    errarr=table2array(errtab);
    disp(' ');
    disp('-----------------------------------------------');
    disp('Average and maximum Euler equation error');
    fprintf('Equ.no.\t\tAvg.\t\tMed.\t\tp75\t\t\tp95\t\t\tp99\t\t\tp99.5\t\tMax.\n');
    for s=1:length(avg_err)
        fprintf('%d\t\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',s,errarr(s,1),errarr(s,2),errarr(s,3), ...
            errarr(s,4),errarr(s,5),errarr(s,6),errarr(s,7));
    end
    
    % plot EE error for these equations
	if ~batchmode
		plotEE_pol=[3,4];
		plotEE_state=[0,0];
		for i=1:length(plotEE_pol)
			points=simseries(:,[6,7]);
			errvals=abs(errmat(1:end-1,plotEE_pol(i)));
			if plotEE_state(i)>0
				itmp=(statvec==plotEE_state(i));
				points=points(itmp,:);
				errvals=errvals(itmp,:);
			end
			model.HelperCollection.scatterPoints2D(points,errvals);
		end
	end
end

% check grid bounds
min_vec=min(simseries(:,state_range));
max_vec=max(simseries(:,state_range));
disp('State bounds:');
disp(mobj.Pfct.SSGrid.StateBounds(:,2:end));
disp('Simulation mins:');
disp(min_vec);
disp('Simulation max:');
disp(max_vec);


% write to file
values=struct2cell(mobj.Params);
paramout=cell2table(values,'RowNames',fieldnames(mobj.Params));
colnames={'mean','std','corrG','corrG_1','corrOm','corrOm_1','corrY','corrY_1','corrY_gr','corrY_gr_1','AC', ...
	'min','p01','p05','p10','p25','p50','p75','p90','p95','p99','max'};
for j=1:4
    tableout_exog=array2table(tabout_exog{j},'RowNames',dispnames,'VariableNames',colnames);
    writetable(tableout_exog,outstats_exog,'WriteRowNames',1,'FileType','spreadsheet','Sheet',j);
    tableout_endog=array2table(tabout_endog{j},'RowNames',dispnames,'VariableNames',colnames);
    writetable(tableout_endog,outstats_endog,'WriteRowNames',1,'FileType','spreadsheet','Sheet',j);    
end
writetable(paramout,outstats_exog,'WriteRowNames',1,'FileType','spreadsheet','Sheet','params');
writetable(paramout,outstats_endog,'WriteRowNames',1,'FileType','spreadsheet','Sheet','params');
writetable(errtab,errstats,'FileType','spreadsheet');
if force_output
    params=mobj.Params;
    disp(['Saving simulation data to .mat file: ',['sim_',resfile,'.mat']]);
    save(['sim_',resfile,'.mat'],'simseries','displist','dispnames','errmat','tabout_exog','tabout_endog','outstats_exog','outstats_endog', ...
                'errstats','errtab','indexmap','NT_ini','NT_sim','smpsel_exog','smpsel_endog','statevec','statsout_exog','statsout_endog','varnames','winsorize','params');
end    

if expdata
    disp(' ');
    disp('Exporting simseries...');
    model.HelperCollection.tableExport(outdata,varnames,simseries);
end

% save model file with stationary state values
save([respath,resfile,'.mat'],'stvstat','-append');

%% Make graphs
if makeBRGraphs
    % file names for graphs (set to empty for no printing)
    % printfiles={'impulse_bankruptcy1_lowgf','impulse_bankruptcy2_lowgf','impulse_bankruptcy3_lowgf','impulse_bankruptcy4_lowgf'};        
    printfiles={[],[],[],[]};        
    % which variables
    brsel1=[indexmap.get('AB'),indexmap.get('Lspr'),indexmap.get('p'),...
           indexmap.get('ndebt_byY'),indexmap.get('bIstart_byY'),indexmap.get('KB')];  
    brsel2=[indexmap.get('cS_gr'),indexmap.get('wS'),...
            indexmap.get('cS_byY'),indexmap.get('C_gr')];  
    brsel3=[indexmap.get('G'),indexmap.get('Om'),...
            indexmap.get('Y_gr'),indexmap.get('XbyY')];
    brsel4=[indexmap.get('Imktlev'),indexmap.get('Bmktlev'),...
            indexmap.get('bind_lamB'),indexmap.get('Ibklev')];
    brsel_all=[brsel1,brsel2,brsel3,brsel4];
    nvar=length(brsel_all);   
    % find all indices with bankruptcy
    % brind=find(simseries(:,2)==mobj.Params.sig2B(2));
    % Expansion to Recession + Crisis
    brind=1+find(simseries(2:end,2)==mobj.Params.sig2_om(2) & simseries(1:end-1,2)==mobj.Params.sig2_om(1) & simseries(1:end-1,1) > 1 );
    % Expansion to Recession
    % brind=1+find(simseries(2:end,1) < 1 & simseries(1:end-1,1) > 1 );
    % Expansion to Recession + No Crisis
    %brind=1+find(simseries(2:end,1) < 1 & simseries(1:end-1,1) > 1 & simseries(2:end,2)==mobj.Params.sig2_om(1) & simseries(1:end-1,2)==mobj.Params.sig2_om(1));
    % Non-crisis to Crisis
    %brind=1+find(simseries(2:end,2)==mobj.Params.sig2_om(2) & simseries(1:end-1,2)==mobj.Params.sig2_om(1));
    % Recession Non-crisis to Recession Crisis
    %brind=1+find(simseries(2:end,1) < 1 & simseries(2:end,2)==mobj.Params.sig2_om(2) & simseries(1:end-1,1) < 1 & simseries(1:end-1,2)==mobj.Params.sig2_om(1));
    %brind=find(brupt(2:end)==1);
    %brind=find(idx_both_bind);
    % High leverage
    % brind=find(simseries(:,indexmap.get('LB'))>1);
    % 5 before and 10 after
    tvec=-2:7;
    stind=brind-2;
    eind=brind+7;
    valrows=logical((stind>0).*(eind<=NT_sim));
    numepi=sum(valrows);
    stind=stind(valrows);
    eind=eind(valrows);
    brseries=zeros(numepi,10,nvar);
    for i=1:numepi
        brseries(i,:,:)=simseries(stind(i):eind(i),brsel_all);
    end
    brseries_gr=zeros(10,nvar);
    for v=1:nvar
        brtmp=squeeze(brseries(:,:,v));
%        brseries_gr(:,v)=quantile(brtmp,0.5);
        brseries_gr(:,v)=mean(brtmp);
    end
    colors={'-b'};
    titles1={'AB','Spread','Tobin''s Q',...
            'Corp.Debt (book)','Deposits','KB'};
    titles2={'cS gr','wS',...
            'cS/Y','Agg C gr'};    
    titles3={'G','Om'...
             'Output Growth','Investment/Output'};
    titles4={'Mkt Lev','B Mkt Lev', ...
             'lamB binds', 'Book Lev'};
    if ~batchmode
        makeImpRes(brseries_gr(:,1:6),tvec,titles1,colors,[2,3],[],printfiles{1});   
        makeImpRes(brseries_gr(:,7:10),tvec,titles2,colors,[2,2],[],printfiles{2});   
        makeImpRes(brseries_gr(:,11:14),tvec,titles3,colors,[2,2],[],printfiles{3});  
        makeImpRes(brseries_gr(:,15:18),tvec,titles4,colors,[2,2],[],printfiles{4}); 
    end
    save([respath,resfile,'.mat'],'brseries_gr','-append');
end