#!/bin/bash
#SBATCH --job-name=test_CA82_postFrac0.3_job   # Job name
#SBATCH --output=test_CA82_postFrac0.3output.log # Standard output log file
#SBATCH --error=test_CA82_postFrac0.3error.log   # Standard error log file
#SBATCH --time=48:00:00           # Time limit hrs:min:sec
#SBATCH --partition=standard      # Partition name
#SBATCH --account=sc139-wetting
#SBATCH --qos=lowpriority
#SBATCH --ntasks=25                # Number of tasks (processes)
#SBATCH --cpus-per-task=32         # Number of CPU cores per task

# Change to the directory
cd test_CA82_postFrac0.3

export OMP_NUM_THREADS=
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_ENV=TRUE   # Optional: print OpenMP environment info at startup

# Run the executable
srun ./run.exe
