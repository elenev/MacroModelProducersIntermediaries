if usejava('desktop')
   clear; 
end
close all;

% ===========================================
% program control
% ===========================================

% parallel processes? set this to zero for just one process
if ~exist('no_par_processes','var')
    no_par_processes=16; %10
end
disp(['PPN: ',num2str(no_par_processes)]);

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



% path to file with experiment definition
if ~exist('exper_path','var')
    exper_path='env_bench_ini0.mat';
	%exper_path='res_2019_9_19_11_9.mat';
end
if ~exist('maxit','var')
    maxit=2;
end
if ~exist('tol_avg','var')
    tol_avg=0.0001;
end
if ~exist('price_zns','var')
    price_zns=true;
end

% print mode
printmode=1;

% dampening and TF bounds
damp=[0.5,0.5];
enforceTFbounds=false;

% ===========================================
% load environment structure
% ===========================================

load(exper_path);

% ===========================================
% outer convergence loop
% ===========================================

% policy iteration
[mobj,failedPoints,dist,distT,distmat,distTmat]=mobj.polIter(maxit,1,printmode,damp,tol_avg,enforceTFbounds);

% price assets in zero net supply
if price_zns
	MAXITER_ZNS=200;
	mobj=mobj.priceZNS(stv{1},tol_avg,MAXITER_ZNS);
end

% make date string
curr_date=clock;
date_str='';
for i=1:5
    date_str=[date_str,'_',num2str(curr_date(i))];
end

if ~exist('outname','var')
    outname=['res',date_str,'.mat'];
end

save(outname,'mobj','stv','failedPoints','dist','distT','distmat','distTmat');
