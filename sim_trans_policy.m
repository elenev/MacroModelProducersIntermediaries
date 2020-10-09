if usejava('desktop')
   clear; 
else
    ver;
end
close all;

%--------------------------------------------------------------------------
% simulation setup
%--------------------------------------------------------------------------

% file with model
respath='./';
outpath='./Results/';
fromfile = 'res_20200910_bench_s130';

if ~exist('resfile','var')
    resfile='res_20200910_xi85_s130';
	%resfile='res_20200910_xi95benchgrid_s130';
end

load([respath,resfile,'.mat']);

% Initial Economy Config
varlist={'simseries','statevec','indexmap','varnames'};
load(['sim_',fromfile],varlist{:});

% set starting point
start_ini=5;
start_shock=0; %[0,3,4];
statevec=statevec(2:end);
startvals=mean(simseries(statevec==start_ini,:));
N_vars=length(startvals);
% states or variables computed at market prices?


% number of periods and burn-in
N_shock=length(start_shock);
N_runs=50000;
NT_sim=75;
NT_ini=0;
NT_sim=NT_sim+1;

start_method = 'endog_cluster';
[~,endog_var_idx] = intersect(varnames,mobj.En_names,'stable');

if strcmp(start_method,'exog_cluster')
	startptmat = zeros(N_runs,mobj.NSTEN+1);
	counter_start=0;
	for ii=1:mobj.Exogenv.exnpt
		n_start_states = floor( sum(statevec==ii)/size(simseries,1) * N_runs );
		startptmat(counter_start+(1:n_start_states),:) = repmat( ...
			[ii, mean( simseries(statevec==ii, endog_var_idx) ) ], n_start_states, 1);
		counter_start = counter_start + n_start_states;
	end
	if size(startptmat,1)<N_runs
		thismean=[5,mean( simseries(statevec==5, endog_var_idx) )];
		startptmat=[startptmat; repmat(thismean,N_runs-size(startptmat,1),1)];
	end
elseif strcmp(start_method,'endog_cluster')
	% cluster endogenous states
	maxclust=1;
	startptmat=[];
	maxclust_per_state = 10;
	simseries_to_cluster = simseries; % replace with simseries(:,statevec==start_ini)
	% to only cluster points for a particular exog. state e.g. 5
	statevec_to_cluster = statevec; % replace as above for a particular state
	[clusters,probs,~,rsq] = clusterEachExogState( statevec_to_cluster, ...
		simseries_to_cluster( :, endog_var_idx), maxclust_per_state, false );
	maxclust = size(clusters,1);
	for c=1:maxclust
		cfrac=probs(c);
		thismean=clusters(c,:);
		disp([num2str(c),': ',num2str(thismean)]);
		thisc=repmat(thismean,floor(N_runs*cfrac),1);
		startptmat=[startptmat; thisc];    
	end
	if size(startptmat,1)<N_runs
		thismean=[5,mean(simseries(statevec==5,endog_var_idx))];
		startptmat=[startptmat; repmat(thismean,N_runs-size(startptmat,1),1)];
	end
	
	startvals = arrayfun(@(i)mean(simseries(statevec==i,:))',unique(clusters(:,1)), ...
		'UniformOutput', false);
	startvals = [startvals{:}]';
else
	error('No start method selected');
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
outfile=['PT_',resfile];

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

    

    SDFmat=zeros(NT_sim,mobj.Exogenv.exnpt);

    fprintf([repmat('.',1,100) '\n\n']);
	
	N_proc = no_par_processes;
	N_runs_per_set = N_runs / N_proc;
	if mod(N_runs,N_proc) > 0
		error('N_runs must be a multiple of N_proc');
	end
	
	tens_simseries_mean = zeros(NT_sim,N_vars,N_proc);
	
	rng(1);
	seeds = randi(1e6,N_proc,1);
    parfor nset=1:N_proc  
	%for nset=1:N_proc
		startptmat_set = startptmat(1:N_proc:end,:);
		% Create shock matrix
		%shmatfull = rand(NT_sim*N_runs,1);
		rng(seeds(nset));
		shmatfull = lhsdesign(N_runs_per_set,NT_sim+1);
		simseries_mean_set = [];
		for n=1:N_runs_per_set
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
			start_ini=startptmat_set(n,1);
			startpt.KB=startptmat_set(n,2);
			startpt.LB=startptmat_set(n,3);
			startpt.WI=startptmat_set(n,4);
			startpt.BG=startptmat_set(n,5);
			startpt=orderfields(startpt,mobj.En_names);
			startpt_vec=model.DSGEModel.structToVec(startpt)';
			startpt_vec=[start_ini,startpt_vec];
			
			[simseries,varnames,~,~,~]=mobj.simulate(NT_sim,NT_ini,startpt_vec,compEEErr,shmat(2:end));
			simseries_orig=simseries;
			varnames_orig=varnames;
			statevec = simseries(:,1);
			%fprintf('Run %d - After simulation \n',n);

			[simseries, varnames] = mobj.computeSimulationMoments(simseries,varnames,[]);

			nvars = length(varnames);
			%fprintf('Run %d - After computation \n',n);
	%         disp(size(startvals))
	%         disp(size(simseries))
			
			tmp = [startvals(start_ini,:); simseries];
			if n==1
				simseries_mean_set = tmp;
			else
				simseries_mean_set = ( (n-1)*simseries_mean_set + tmp ) / n;
			end
			if mod((nset-1)*N_runs_per_set+n,N_runs/100)==0
				%disp([num2str(n),'/',num2str(N_runs),': ',num2str(round(1000*n/N_runs)/10),'% complete']);
				fprintf('\b|\n');
			end
		end
		tens_simseries_mean(:,:,nset) = simseries_mean_set;
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

    %simseries_median{s} = median(tens_simseries,3);
    simseries_mean{s} = mean(tens_simseries_mean,3);
    %simseries_std{s} = std(tens_simseries,[],3);

end

save(outfile,'simseries_mean','indexmap','NT_sim','N_shock');

% if usejava('desktop')
%    plot_trans; 
% end
