#!/bin/bash
#SBATCH --job-name=test_CA80_postFrac0.1_job   # Job name
#SBATCH --output=test_CA80_postFrac0.1/output.log # Standard output log file
#SBATCH --error=test_CA80_postFrac0.1/error.log   # Standard error log file
#SBATCH --time=48:00:00           # Time limit hrs:min:sec
#SBATCH --partition=standard      # Partition name
#SBATCH --account=sc139-wetting
#SBATCH --qos=lowpriority
#SBATCH --ntasks=1                # Number of tasks (processes)
#SBATCH --cpus-per-task=32        # Number of CPU cores per task

# Change to the directory
cd "$SLURM_SUBMIT_DIR/test_CA80_postFrac0.1"

# Run the executable
srun ./run.exe
