#!/bin/bash
#SBATCH --output=./slurms/WCrun_slurm-%j.out
#SBATCH --job-name="WC"
#SBATCH --nodelist=tethys-3g-21,tethys-3g-22
#SBATCH --time=48:00:00
#SBATCH --mail-type=all

set -e

module load matlab
echo "Running MATLAB"
nohup matlab -nodisplay -nosplash -nodesktop -r "run('resolution_testing.m'); exit" > WCResOutput.log # 2>&1 &
echo "Done"
echo ""