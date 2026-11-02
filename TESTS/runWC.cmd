#!/bin/bash
#SBATCH --output=./slurms/WCrun_slurm-%j.out
#SBATCH --job-name="WC20"
#SBATCH --time=72:00:00
#SBATCH --mail-type=all

set -e

module load matlab
echo "Running MATLAB"
nohup matlab -nodisplay -nosplash -nodesktop -r "run('resolution_testing20.m'); exit" > WCResOutput20dt.log # 2>&1 &
echo "Done"
echo ""