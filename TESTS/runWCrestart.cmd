#!/bin/bash
#SBATCH --output=./slurms/WCrun_slurm-%j.out
#SBATCH --job-name="WCintf"
#SBATCH --time=72:00:00
#SBATCH --mail-type=all

set -e

module load matlab
echo "Running MATLAB"
nohup matlab -nodisplay -nosplash -nodesktop -r "run('restart_testing.m'); exit" > WCRestartOutput20intfine2.log # 2>&1 &
echo "Done"
echo ""