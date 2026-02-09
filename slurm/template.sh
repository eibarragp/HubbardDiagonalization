#!/bin/bash
# Basic Info
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1

# Customizable Params
#SBATCH --job-name={{ job_name }}
#SBATCH --output={{ log_file }}
#SBATCH --mail-type={{ mail_type }}
#SBATCH --time={{ time_limit }}
#SBATCH --mem={{ memory_limit }}
#SBATCH --cpus-per-task={{ cpus_per_task }}

# Possibly null params
#{{ mail_user }}
#{{ dependencies }}
#{{ array_info }}

# Setup Environment
date;hostname;pwd
module load julia
echo "Running task $SLURM_ARRAY_TASK_ID on $SLURM_CPUS_PER_TASK cores"
# Use our own package repository bc the system one is write-protected
# But allow using the system one to enable bundled packages to work (Expanded by the empty entry after ':')
export JULIA_DEPOT_PATH="$NLCE_HOME/.julia:"
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK  # Use all cores!
cd $NLCE_HOME/HubbardDiagonalization

{{ command }}
