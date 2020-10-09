#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2:49:00
#SBATCH --mem=62GB
#SBATCH --job-name=xpcombined
#SBATCH --output=xp-%A_%a.out


if [ -z "${EXPERDEF+xxx}" ];  then
    EXPERDEF=20200910
fi
EXPER_FILE=experdef_"$EXPERDEF".m
GUESS_STATUS=no_guess
if [ -z "${USEINI+xxx}" ];  then
    USEINI=1
fi
MAXITERINI=100

PPN=$SLURM_CPUS_PER_TASK

module purge
module load matlab/2017b
if [ "$SLURM_JOBTMP" == "" ]; then
    export SLURM_JOBTMP=/state/partition1/$USER/$$
    mkdir -p $SLURM_JOBTMP
fi
 
export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)

# Create run directory
RUNDIR=$SCRATCH/InterProd/run-${SLURM_ARRAY_JOB_ID/.*}-${SLURM_JOB_NAME/.*}-${SLURM_ARRAY_TASK_ID/.*}
mkdir $RUNDIR

# Copy files to $RUNDIR
DATADIR=$SCRATCH/code/InterProd/version20200310
rsync -az --exclude '*.mat' --exclude '*.xls' "$DATADIR"/ "$RUNDIR"
cd $RUNDIR

# Delete any existing output files
rm -f output.log
rm -f output_sim.log

matlab -nodisplay -nodesktop -r "run('$EXPER_FILE');exit;" >> output.log
ECONOMY="`sed -n "$SLURM_ARRAY_TASK_ID"p experdef_"$EXPERDEF".txt`"

if [ -z "$ECONOMY" ]; then
    echo "SLURM array index exceeds number of economies"
    exit 1
fi

# Define an environment variable for the number of available cores
export NTHREADS=$SLURM_CPUS_ON_NODE

echo "Requested economy $ECONOMY" >> output.log

# ini0
GUESS_NAME=none
if [ $USEINI -eq 1 ]; then
	EXPER_NAME="$ECONOMY"_ini0
	EXPER_NAME_i100="$ECONOMY"_i100
	RES_NAME=res_"$EXPERDEF"_"$EXPER_NAME_i100".mat
	MAXITER=$MAXITERINI
	echo "Starting ini0 grid run: $MAXITER iterations using guess: $GUESS_NAME" >> output.log
	matlab -nodisplay -nodesktop -r "experdef_file='$EXPER_FILE'; expername='$EXPER_NAME'; guess_mode='$GUESS_STATUS'; guess_path='$GUESS_NAME'; ststonly=0; run('main_create_env.m');exit;" >> output.log
    ENV_NAME=env_"$EXPER_NAME".mat
	matlab -nodisplay -nodesktop -r "maxit=$MAXITER;no_par_processes=$PPN;outname='$RES_NAME';exper_path='$ENV_NAME';price_zns=false;run('main_run_exper.m');exit;" >> output.log
	cp $RES_NAME $DATADIR/
	GUESS_NAME=$RES_NAME
	GUESS_STATUS=guess
fi

MAXITER=30

SUFFIX=s130
RES_NAME=res_"$EXPERDEF"_"$ECONOMY"_"$SUFFIX".mat
cp "$DATADIR"/"$GUESS_NAME" "$RUNDIR"/"$GUESS_NAME"

EXPER_NAME="$ECONOMY"
ENV_NAME=env_"$EXPER_NAME".mat

echo "Starting fine grid run: $MAXITER iterations using guess: $GUESS_NAME" >> output.log
matlab -nodisplay -nodesktop -r "experdef_file='$EXPER_FILE'; expername='$EXPER_NAME'; guess_mode='$GUESS_STATUS'; guess_path='$GUESS_NAME'; ststonly=0; run('main_create_env.m');maxit=$MAXITER;no_par_processes=$PPN;outname='$RES_NAME';exper_path='$ENV_NAME';run('main_run_exper.m');exit;" >> output.log
cp $RES_NAME $DATADIR/


# Simulate
RES_NAME=res_"$EXPERDEF"_"$ECONOMY"_"$SUFFIX"
matlab -nodisplay -nodesktop -r "batchmode=true; resfile='$RES_NAME'; run('sim_stationary.m');exit;" >> output_sim.log
cp sim_"$RES_NAME".mat $DATADIR/
cp Results/statsexog_"$RES_NAME".xls $DATADIR/Results/
cp Results/statsendog_"$RES_NAME".xls $DATADIR/Results/
cp Results/errstats_"$RES_NAME".xls $DATADIR/Results/

# IRFs
matlab -nodisplay -nodesktop -r "resfile='$RES_NAME';sim_trans_cluster;exit;" >> output_sim.log
cp GR_"$RES_NAME".mat $DATADIR/
cp  Results/statsirf_"$RES_NAME".xls $DATADIR/Results/

# Clean up
rm -rf $data_location
rm -rf $SLURM_JOBTMP/*


