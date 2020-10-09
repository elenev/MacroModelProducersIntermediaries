This file provides step by step instructions to replicate all results in the paper, as submitted for publication to Econometrica.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Part 1: Solve and Simulate the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This part can be done locally or on the cluster. There are separate instructions for each method below. If recomputing all 93 economies required for the paper, it is much faster to do it on the cluster.

------------------------------
Part 1a: Computing locally
------------------------------

Repeat the steps below for each economy in experdef_20200910.txt. The definitions of these economies are in experdef_20200910.m. Let the placeholder "econ" represent the name of the economy from this text file.

First, run the economy for up to 100 iterations on the coarse grid.
1.  Configure main_create_env.m as follows (leave other variables as is):
	- expername = 'econ_ini0';
	- guess_mode = 'no_guess';
2.  Run main_create_env.m and it will produce a file called env_econ_ini0.mat
3.  Configure main_run_exper.m as follows (leave other variables as is):
	- no_par_processes = XX; % where XX is the number of processors on the machine you're running it on
	- exper_path = 'env_econ_ini0.mat';
	- maxit = 100;
	- price_zns = false;
4.  Run main_run_exper.m. On a machine with 16 cores, this should take about 30 min. It will create a file named res_TIMESTAMP.mat, where TIMESTAMP is the time at which the calcuation is finished expressed as YYYY_MM_DD_hh_mm
5.  Rename the res_TIMESTAMP.mat file to res_20200910_econ_i100.mat.

Next, run the economy for up to 100 iterations on the fine grid
6.  Configure main_create_env.m as follows (leave other variables as is):
	- expername = 'econ';
	- guess_mode = 'guess';
	- guess_path = 'res_20200910_econ_i100';
7.  Run main_create_env.m and it will produce a file called env_econ.mat.
8.  Configure main_run_exper.m as follows (leave other variables as is):
	- no_par_processes = XX; % where XX is the number of processors on the machine you're running it on
	- exper_path = 'env_econ.mat';
	- maxit = 30;
	- price_zns = true;
9.  Run main_run_exper.m. On a machine with 16 cores, this should take about 40 min. It will create a file named res_TIMESTAMP.mat, where TIMESTAMP is the time at which the calcuation is finished expressed as YYYY_MM_DD_hh_mm
10. Rename the res_TIMESTAMP.mat file to res_20200910_econ_s130.mat.

Next, simulate. The model must be solved and res* file must exist.
11. Configure sim_stationary.m as follows (leave other variables as is):
	- resfile = 'res_20200910_econ_s130';
12. Run sim_stationary.m. If you are running sim_stationary having run it before during the current MATLAB session, run "clear" beforehand. This operation will create the following files:
	- sim_res_20200910_econ_s130.mat: all simulation results incl. full series, statistics, errors, and parameters
	- Results/statsexog_res_20200910_econ_s130.xls: statistics using exogenous subsampling to define crises (one sample per worksheet)
	- Results/statsendog_res_20200910_econ_s130.xls: statistics using endogenous subsampling to define crises (one sample per worksheet)
	- Results/errstats_res_20200910_econ_s130.xls: statistics of EE, VF, and TF errors

Next, compute impulse response functions. Model must be solved and simulated. Both res* and sim_res* files must exist.
13. Configure sim_trans_cluster.m as follows (leave other variables as is):
	- resfile = 'res_20200910_econ_s130';	
14. Run sim_trans_cluster.m. On a machine with 16 cores, this should take about 15 min. It will create the following files:
	- GR_res_20200910_econ_s130.mat: mean, median, and sd of IRF paths for each of 4 shocks (no shock, non-fin rec, fin rec, and pure uncertainty)
	- statsirf_res_20200910_econ_s130.mat: means of IRF paths (one sheet per shock)
	
------------------------------------
Part 1b: Computing on the cluster
------------------------------------

The following instructions work for a cluster that uses the SLURM job manager, has an environment variable $SCRATCH defined and pointed to a writeable folder. They are an ALTERNATIVE to Part 1a. 

1. Configure the cluster job file slurmcombined0312.sh processor, memory, and time requests. Current config of 2 hrs, 49 min is generous with 20 cores, but better to be safe. At least 2 gigs of memory per core.
2. Copy the folder ecma to the cluster into $SCRATCH/code/InterProd (e.g. using WinSCP)
3. If it doesn't exist, create $SCRATCH/InterProd folder.
3. Connect to the cluster's terminal (e.g. using PuTTY) and cd into that folder.
4. Run "sbatch -a 2-94 slurmcombined0910.sh." This will submit jobs to solve, simulate, and compute IRFs for each economy used in the paper. This will create for each economy in experdef_20200910.txt (let "econ" denote the economy name) 
	- res_20200910_econ_i100.mat: results file after just the coarse grid iterations
	- res_20200910_econ_s130.mat: results file
	- sim_res_20200910_econ_s130.mat: all simulation results incl. full series, statistics, errors, and parameters
	- GR_res_20200910_econ_s130.mat: mean, median, and sd of IRF paths for each of 4 shocks (no shock, non-fin rec, fin rec, and pure uncertainty)
	- Results/statsexog_res_20200910_econ_s130.xls: statistics using exogenous subsampling to define crises (one sample per worksheet)
	- Results/statsendog_res_20200910_econ_s130.xls: statistics using endogenous subsampling to define crises (one sample per worksheet)
	- Results/errstats_res_20200910_econ_s130.xls: statistics of EE, VF, and TF errors
	- statsirf_res_20200910_econ_s130.mat: means of IRF paths (one sheet per shock)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			Part 2: Compute and Save Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This part pre-supposes that the res*, sim_res*, GR_res* and Excel files for each economy created by Part 1 exist.

The steps below re-create all the results used in the paper.

1.  Compute CEV welfare by running "welfare(10)" and "welfare_appd5(10)" in MATLAB. "10" is referring to the number of data clusters to create for evaluating transitions from benchmark to one of the other economies. welfare.m computes CEV welfare relative to benchmark. welfare_appd5.m computes CEV welfare at xi=95 relative to xi=91 for several different parameter combinations, creating the welfare numbers in Table D.2.

	Note: this requires res* files for each economy, which are very large (~350MB). If they're stored on the cluster and not locally, log into the cluster and issue the following commands to run the welfare function on the cluster as a interactive job:
	
	cd $SCRATCH/code/InterProd/version20200320
	srun -c8 -t2:00:00 --mem=32000 --pty /bin/bash
	module load matlab/2018b
	matlab -nodesktop -nodisplay -r "welfare(10);"
	matlab -nodesktop -nodisplay -r "welfare_appd5(10);"
	
	This will create Results/welfare_20200910_bench_s130.xls, containing the CEV welfare changes from going to bench to altrnative economies for each computed economy in the folder.
	
2.  Run writeStats.m. This will create Results/simres.mat, which contains a variable called "mapDict," a hash-map of every simulation and IRF statistic and model parameter. They can be retrieved by running mapDict("key"), where "key" has the following format:
	[command] - [econ] - [subsmaple] - [variable] - [statistic]
	
	command:
		- sim: simulation statistic
		- comp: simulation statistic relative to bench (either level or percent change)
		- irf: IRF value at t=1
		
	econ: name of the economy (as in experdef_20200910.txt)
	
	subsample: if command is "sim" or "comp"
		- u: unconditional
		- e: expansions (not available if command = "sim" or "comp")
		- r: recessions
		- c: crises
		- l: low-uncertainty (not available if command = "sim" or "comp")
		
	variable: all variables reported by sim_stationary as well as "cvwelfare." Different variables are reported differently e.g. in percentage terms. See the code of writeStats.m for details.
	
	statistic: for command = "sim" and "comp", all statistics reported by sim_stationary e.g. "mean," "std," "p50". For a few ratios, there is "nmean" which is ratio of means (whereas "mean" is mean of ratios). 
	
		For "command" = irf, choices are "mean" (level at t=1) or "change" (deviation of t=1 from t=0 value, either in levels or percent)
		
	After simres.mat has been created, the build.m LaTeX parser uses these values to fill in placeholders in the *.tex file (see below). 
		
3.  Create benchmark economy IRFs (Figures 2 and 3). Configure plot_trans.m as follows:
	- outfile = ['GR_',resfile];
	- plot_shocks = 1:3;
	Run plot_trans.m. This will create figures in Results/GR_res_20200910_bench_s130_IRF#.(pdf|eps)with # = {1,2,3,4}.
	

4.  Create IRFs comparing financial crises in various economies (Figure 4). Configure plot_trans_finrec.m as follows:
	- econ{1} = 'FLbench';
	- econ{2} = 'FLwithdef';
	- econ{3} = 'bench';
	Comment out econ{4} if needed and run. This will create figures in Results/GR_res_20200910_finrec_FLbenchFLwithdefbench_s130_IRF1.(pdf|eps)
	
5.  Create IRFs comparing financial crises in benchmark vs counter-cyclical cap reqs economies (Figure 7). Configure plot_trans_finrec.m as follows:
	- econ{1} = 'bench';
	- econ{2} = 'xi9195';
	Comment out econ{3} and econ{4} if needed and run. This will create figures in Results/GR_res_20200910_finrec_benchxi9195_s130_IRF1.(pdf|eps)

6.  Create graphs comparing equilibrium quantities and welfare across macroprudential experiments. run makeMacropruPlots.m. This will create figures in Figures/macropru*.(pdf|eps).

7.  Compute transitions from benchmark to other capital requirements (Figure 8). For each of econ = xi85 and econ = xi95benchgrid,
		a. Configure sim_trans_policy.m as follows:
			- resfile = 'res_20200910_econ_s130';
		b. Run sim_trans_policy.m. This will create PT_res_20200910_econ_s130.mat.
		c. Configure plot_trans_policy.m as follows:
			- resfiles = {'res_20200910_econ_s130'};
			- labels = {'\xi=XX'}; % where XX is the value of xi in that economy
		   This will create figures in Results/PT_res_20200910_econ_s130.(pdf|eps).

8.  Create policy function plots (Figure B.1) by running plotPolicyFunctions.m. This will create figures in Results/polPlot#.(eps|pdf), where # = {,1,2}.

9.  Create state space histograms with grid point lines (Figure B.2) by running compareStateSpaces.m. This will create figures in Results/stateSpaceHistograms.(pdf|eps).

10. Create Gentskow-Shapiro sensitivity analysis bar. This will create figures in Figures/GS_1.(pdf|eps).

11. Create benchmark economy IRFs comparing financial recession to pure uncertainty (Figures D.1 and D.2). Configure plot_trans.m as follows:
	- outfile = ['GR_unc_',resfile];
	- plot_shocks = [1,4,3];
	Run plot_trans.m. This will create figures in Results/GR_unc_res_20200910_bench_s130_IRF#.(pdf|eps)with # = {1,2,3,4}.
	
12. Create a figure showing how the credit spread and expected excess return changes with intermediary wealth (Figure D.3) by running makeEERhist.m. This will create figures in Results/WI_EERhistres_20200910_bench_s130.(pdf|eps).

13. Create a figure showing the IRFs after an unanticipated shock to intermediary wealth ("mortgage" shock, Figure D.4).  For each of econ = bench and econ = benchriskysafesigI0,
		a. Configure sim_trans_mit.m as follows:
			- resfile = 'res_20200910_econ_s130';
		b. Run sim_trans_mit.m. This will create MIT_res_20200910_econ_s130.mat.
	Run plot_trans_mit.m. This will create figures in Results/GR_res_20200910_mit_benchbenchriskysafesigI0s130_IRF#.(pdf|eps), where # = {1,2}.
	
14. Create a table showing the robustness of macropru results to alternative parameters (Table D.2). 
	a. Copy every Excel file containing *xi91* and *xi95* from Results into Results_D5. Also copy welfareappd5_20200910_s130.xls.
	b. Open table_d2.xlsm, enabling macros.
	c. In the "Overview" sheet, click "Update from sim_stationary." The macro will run for a minute or so.
	d. When done, go to "Drivers of Macropru." Use the Excel2Latex marcro to convert the two shaded tables into Panel A (top) and Panel B (bottom)