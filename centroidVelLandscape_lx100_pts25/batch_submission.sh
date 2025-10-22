#!/bin/bash

# Loop through all directories
for dir in test_CA*; do
  # Check if the directory exists and contains an executable file
  if [[ -d "$dir" && -x "$dir/run.exe" ]]; then
    # Create a SLURM batch script
    cat <<EOF > "${dir}job.slurm"
#!/bin/bash
#SBATCH --job-name=${dir%/}_job   # Job name
#SBATCH --output=${dir}output.log # Standard output log file
#SBATCH --error=${dir}error.log   # Standard error log file
#SBATCH --time=48:00:00           # Time limit hrs:min:sec
#SBATCH --partition=standard      # Partition name
#SBATCH --account=sc139-wetting
#SBATCH --qos=lowpriority
#SBATCH --ntasks=25                # Number of tasks (processes)
#SBATCH --cpus-per-task=32         # Number of CPU cores per task

# Change to the directory
cd $dir

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_ENV=TRUE   # Optional: print OpenMP environment info at startup

# Run the executable
srun ./run.exe
EOF

    # Submit the job script
    sbatch "${dir}job.slurm"
  else
    echo "Skipping $dir: Not a directory or no executable found."
  fi
done

echo "Job submission complete."
