#!/bin/bash
# Basic Info
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1

# Customizable Params
#SBATCH --job-name={{ job_name }}
#SBATCH --output={{ log_file }}
#SBATCH --mail-type=END,FAIL
#SBATCH --time=1:0
#SBATCH --mem=1gb
#SBATCH --cpus-per-task=1

# Possibly null params
#{{ mail_user }}

date;hostname;pwd
CONFIG_FILE={{ project_root }}/SimulationConfig.toml
OUTPUT_DIR={{ nlce_home }}/output
U=$(grep 'U = ' $CONFIG_FILE | awk '{print $3}')

mkdir $OUTPUT_DIR/{{ output_group }}_U$U
mv $OUTPUT_DIR/{{ output_name_base }}* $OUTPUT_DIR/{{ output_group }}_U$U

if [[ "$U" == "{{ U_max }}" ]]
then
	exit 0
else
	sed -i "s/U = $U/U = $(((U + {{ U_step }})))/" $CONFIG_FILE
	{{ nlce_home }}/slurm/{{ script_name }}.sh 'no_reset'
fi
