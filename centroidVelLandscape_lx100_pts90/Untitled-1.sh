#!/bin/bash

# Loop through all directories
for CA in 94 96 98 100; do 
  for pf in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do
    dir="test_CA${CA}_postFrac${pf}"

    # Check if the directory exists and contains an executable file
    if [[ -d "$dir" && -x "$dir/run.exe" ]]; then
      # Create a SLURM batch script
      cat <<EOF > "${dir}/job.slurm"
#!/bin/bash
#SBATCH --job-name=${dir}_job   # Job name
#SBATCH --output=${dir}/output.log # Standard output log file
#SBATCH --error=${dir}/error.log   # Standard error log file
#SBATCH --time=48:00:00           # Time limit hrs:min:sec
#SBATCH --partition=standard      # Partition name
#SBATCH --account=sc139-wetting
#SBATCH --qos=lowpriority
#SBATCH --ntasks=1                # Number of tasks (processes)
#SBATCH --cpus-per-task=32        # Number of CPU cores per task

# Change to the directory
cd "\$SLURM_SUBMIT_DIR/$dir"

# Run the executable
srun ./run.exe
EOF

      # Submit the job script
      sbatch "${dir}/job.slurm"
    else
      echo "Skipping $dir: Not a directory or no executable found."
    fi
  done
done

echo "Job submission complete."
