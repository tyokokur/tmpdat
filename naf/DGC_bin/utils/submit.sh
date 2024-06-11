#!/bin/sh
# FILENAME: submit.sh

## INPUT ##
#SBATCH --time=00:01:00          # Maximum 96:00:00 for shared, wholenode
#SBATCH --cpus-per-task=1        # ntasks * cpus-per-task <= CPU cores / node (max 128)
##SBATCH --mem-per-cpu=5G         # Alternate method for determining CPU cores (1.85G each)
###########

#SBATCH -A chm240035             # TJY 2024
#SBATCH -p shared                # shared, wholenode
#SBATCH -o logs/slurm-%j.out     # Name of stdout output file
#SBATCH -e logs/slurm-%j.out     # Name of stderr output file
#SBATCH --open-mode=append
#SBATCH --nodes=1                # Must be 1 for OpenMP jobs
#SBATCH --ntasks-per-node=1      # Tasks per node
#SBATCH --mail-user=takashi.yokokura@berkeley.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND,REQUEUE,STAGE_OUT

# Clean (Same as make clean but don't destroy current output file)
echo "*** Clean ***"
shopt -s extglob
rm -vf "("* ph*.dat it*.dat el*.dat printout.dat
rm -vf ./logs/!(*$SLURM_JOBID*)
shopt -u extglob
echo ""

# Setup 
echo "*** Setup ***"
module purge     # Unload all loaded modules and reset everything to original state
module load gcc  # 
module load fftw # 
module list      # List currently loaded modules
hostname         # Print the hostname of the compute node on which this job is running
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK # Set thread count
echo "OMP_NUM_THREADS: ${SLURM_CPUS_PER_TASK}"
export STDOUT_REDIR=1 # Redirect stdout and refresh
module load monitor 
monitor cpu percent > logs/cpu-percent.log & CPU_PID=$!
echo "Monitor to ./logs/cpu-percent.log"
echo ""

## RUN ##
echo "*** Run ***"
make
echo "Running..."
srun ./tjygo
echo "Run complete"
echo ""
#########

# Monitor off
kill -s INT $CPU_PID # shut down the resource monitors
