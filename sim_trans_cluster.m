if usejava('desktop')
   clear; 
else
    ver;
end
close all;

%--------------------------------------------------------------------------
% simulation setup
%--------------------------------------------------------------------------

% excel export
excel_export_means = true;

% file with model
respath='./';
outpath='./Results/';
if ~exist('resfile','var')
    resfile='res_20200310_bench_s130';
end

load([respath,resfile,'.mat']);
params=mobj.Params;

% Initial Economy Config
varlist={'simseries','statevec','indexmap','varnames'};
load(['sim_',resfile],varlist{:});

% set starting point
start_ini=5;
start_shock=[0,3,4,6]; 
%start_shock=[0,6,4]; 
statevec=statevec(2:end);
startvals=mean(simseries(statevec==start_ini,:));
N_vars=length(startvals);
% states or variables computed at market prices?


% number of periods and burn-in
N_shock=length(start_shock);
N_runs=10000;
NT_sim=25;
NT_ini=0;
NT_sim=NT_sim+1;

% cluster endogenous states
maxclust=15;
enstatemat=simseries(statevec==start_ini, ...
	[ indexmap.get('KB'), ...
	  indexmap.get('LB'), ...
	  indexmap.get('WI'), ...
	  indexmap.get('BG') ]);
cindex=clusterdata(enstatemat,'criterion','distance','maxclust',maxclust,'linkage','weighted');
sttot=size(enstatemat,1);
startptmat=[];
for c=1:maxclust
    cfrac=sum(cindex==c)/sttot;
    thismean=mean(enstatemat(cindex==c,:),1);
    disp([num2str(c),': ',num2str(thismean)]);
    thisc=repmat(thismean,floor(N_runs*cfrac),1);
    startptmat=[startptmat; thisc];    
end
if size(startptmat,1)<N_runs
    thismean=mean(enstatemat);
    startptmat=[startptmat; repmat(thismean,N_runs-size(startptmat,1),1)];
end


% report levels or grwoth rates for output variables
reportLevels=1;

% compute Euler equation error?
compEEErr=0;
% make graphs grayscale
grayscale=0;

% Compute term premium (slow!)
term_premium=0;

% output table file
% outfile=['GR_',resfile];
% outstats = ['statsirf_',resfile];
outfile=['GR_',resfile];
outstats = ['statsirf_',resfile];

varnames_store = varnames;

simseries_median = cell(N_shock,1);
simseries_mean = cell(N_shock,1);
simseries_std = cell(N_shock,1);

if usejava('desktop')
    if ~exist('no_par_processes','var')
        no_par_processes=16;
    end
    %disp(['PPN: ',num2str(no_par_processes)]);

    cp=gcp('nocreate');
    if ~isempty(cp) 
        if  cp.NumWorkers~=no_par_processes
            delete(cp);
            if no_par_processes>0
                parpool(no_par_processes);
            end
        end
    else
        if no_par_processes>0
            parpool(no_par_processes);
        end
    end
else
    %disp(['PPN Request: ',getenv('NTHREADS')]);
    %data_location = getenv('DATA_LOCATION');
    scheduler = parcluster('local');
    %scheduler.JobStorageLocation = data_location;
    parpool(scheduler, str2double(getenv('NTHREADS')));
    %disp(['PPN Launched: ',num2str(scheduler.NumWorkers)]);
end

for s=1:N_shock

    disp(['Shock ',num2str(s),' of ',num2str(N_shock)]);
    tens_simseries = zeros(NT_sim,N_vars,N_runs);

    % compute entry of random number matrix that sets first state
    % deterministically to start_shock
    if start_shock(s)>0
        transprob=cumsum(mobj.Exogenv.mtrans(start_ini,:));
        shock_prob=transprob(start_shock(s));
        if start_shock(s)>1
            shock_prob_minus=transprob(start_shock(s)-1);
        else
            shock_prob_minus=0;
        end
        rvar_next=(shock_prob+shock_prob_minus)/2;
    end

    % Create shock matrix
    rng(1);
    %shmatfull = rand(NT_sim*N_runs,1);
    shmatfull = lhsdesign(N_runs,NT_sim);

    SDFmat=zeros(NT_sim,mobj.Exogenv.exnpt);

    fprintf([repmat('.',1,100) '\n\n']);

    parfor n=1:N_runs       
	%for n=1:N_runs  
        %--------------------------------------------------------------------------
        % start simulation
        %--------------------------------------------------------------------------
        %fprintf('Run %d - Start \n',n);
        % simulate
        shmat = shmatfull(n,:)';
        if start_shock(s)>0
            shmat(1)=rvar_next;
        end
        
        startpt=struct;
        startpt.KB=startptmat(n,1);
        startpt.LB=startptmat(n,2);
        if params.CBS
            startpt.BG=startptmat(n,3);
        else
            startpt.WI=startptmat(n,3);
            startpt.BG=startptmat(n,4);
        end
        startpt=orderfields(startpt,mobj.En_names);
        startpt_vec=model.DSGEModel.structToVec(startpt)';
        startpt_vec=[start_ini,startpt_vec];
        
        [simseries,varnames,~,~,~]=mobj.simulate(NT_sim,NT_ini,startpt_vec,compEEErr,shmat);
        simseries_orig=simseries;
        varnames_orig=varnames;
        statevec = simseries(:,1);
        %fprintf('Run %d - After simulation \n',n);
      
        [simseries, varnames] = mobj.computeSimulationMoments(simseries,varnames,[]);

        nvars = length(varnames);
        %fprintf('Run %d - After computation \n',n);
%         disp(size(startvals))
%         disp(size(simseries))
        tens_simseries(:,:,n) = [startvals; simseries];
        if mod(n,N_runs/100)==0
            %disp([num2str(n),'/',num2str(N_runs),': ',num2str(round(1000*n/N_runs)/10),'% complete']);
            fprintf('\b|\n');
        end

    end
    fprintf('\n');
    varnames = varnames_store;
    nvars = length(varnames);

    % make HashMap with mapping of names to indices
    indexmap=java.util.HashMap;
    for i=1:nvars
        indexmap.put(varnames{i},i);
    end
%     varst=zeros(length(startpt_vec)-1,1);
%     for i=1:length(startpt_vec)-1
%         varst(i)=indexmap.get(mobj.En_names{i});
%     end

    %save(outfile,'tens_simseries','indexmap');

    simseries_median{s} = median(tens_simseries,3);
    simseries_mean{s} = mean(tens_simseries,3);
    simseries_std{s} = std(tens_simseries,[],3);

end

save(outfile,'simseries_mean','simseries_median','simseries_std','indexmap','NT_sim','N_shock');

if excel_export_means
	[~,ia,~]=unique(varnames);
	colnames = arrayfun(@(i)sprintf('t%02d',i),0:NT_sim-1,'UniformOutput',false);
	
	for ii=1:N_shock
		tableout=array2table(simseries_mean{ii}');
		tableout = tableout(ia,:);
		tableout.Properties.RowNames = varnames(ia);
		tableout.Properties.VariableNames = colnames;
		writetable(tableout,[outpath,outstats],'WriteRowNames',1,'FileType','spreadsheet','Sheet',ii);
	end
end

% if usejava('desktop')
%    plot_trans; 
% end
