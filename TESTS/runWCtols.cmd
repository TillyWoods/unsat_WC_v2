#!/bin/bash
#SBATCH --output=./slurms/WCrun_slurm-%j.out
#SBATCH --job-name="WCtols5"
#SBATCH --time=72:00:00
#SBATCH --mail-type=all

set -e

#tolphi0 = [1e-3,1e-4,1e-5];
#tolS1 = $1

module load matlab
echo "Running MATLAB"
nohup matlab -nodisplay -nosplash -nodesktop -r "run('tolerance_testing(1e-5)'); exit" > WCtolsOutputtolS11e-5.log # 2>&1 &
echo "Done"
echo ""