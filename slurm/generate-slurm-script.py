#region Setup

import argparse
import json
import math
import os
import re

# Calculate the project root and NLCE home directory based on
# 1) Environment variables if set
# 2) The location of this script

script_dir = os.path.dirname(os.path.abspath(__file__))
on_cluster = script_dir.startswith('/hopper')

NLCE_HOME = os.getenv('NLCE_HOME',
					  os.path.abspath(os.path.join(script_dir, '..', '..')
					  if on_cluster else os.path.join(script_dir, '..', 'nlce_home')))
PROJECT_ROOT = os.getenv('PROJECT_ROOT', os.path.abspath(os.path.join(script_dir, '..')))

# Always run from the project root
julia_base_command = f'julia --project=. src/Main.jl'

cli_parser = argparse.ArgumentParser(
	prog="generate-slurm-script",
	description="Generate a SLURM script for performing an NLCE run"
)

cli_parser.add_argument('clusterfile', type=str, help="Path to the cluster info file")
cli_parser.add_argument('jobname', type=str, help="Name of the SLURM job")
cli_parser.add_argument('--test', metavar='order', type=int, help="If specified, create a test script that runs one cluster of the specified order.")
cli_parser.add_argument('--mail-user', type=str, action='append', help="Email address to send notifications to. May be specified multiple times for multiple addresses.")
cli_parser.add_argument('--no-mail', action='store_true', help="If specified, do not send email notifications.")

args = cli_parser.parse_args()

#endregion

#region Helpful Functions

def binomial(n, k):
	return math.factorial(n) / (math.factorial(k) * (math.factorial(n - k)))

def input_or_default(input_value, default_value):
	value = input(f'{input_value} [{default_value}]: ')
	return default_value if len(value.strip()) == 0 else value

def generate_script_from_template(filename, job_params):
	template_file = os.path.join(NLCE_HOME, f'{PROJECT_ROOT}/slurm/template.sh')
	with open(template_file, 'r') as f:
		template = f.read()
	for key, value in job_params.items():
		# Use a lambda to ensure that the replacement value is treated as a literal string
		template = re.sub(rf'{{{{\s*{key}\s*}}}}', lambda x: value, template)

	with open(filename, 'w') as f:
		f.write(template)

def input_time(prompt):
	value = input(f'{prompt}: ').strip().lower()

	if m := re.match(r'(?:(?:(?:(\d):)?(\d):)?(\d):)?(\d)$', value):
		def parse_pos(pos, mult):
			val = m.group(pos)
			return 0 if val is None else mult * int(val)

		time = parse_pos(4, 1)
		mult = 60
		time += parse_pos(3, mult)
		mult *= 60
		time += parse_pos(2, mult)
		mult *= 24
		time += parse_pos(1, mult)
		return time
	if m := re.match(r'(\d+)(dhms)$', value):
		value = int(m.group(1))
		dur = m.group(2)
		if dur == 'd':
			mult = 24 * 60 * 60
		elif dur == 'h':
			mult = 60 * 60
		elif dur == 'm':
			mult = 60
		else:
			mult = 1

		return value * mult

	raise ValueError('Invalid Time Specification!')

def format_time(time):
	if time < 0:
		raise ValueError('Cannot Format a Negative Time!')
	scale = 24 * 60 * 60
	days = time // scale
	time %= scale
	scale /= 24
	hours = time // scale
	time %= scale
	scale /= 60
	minutes = time // scale
	time %= scale
	scale /= 60
	seconds = time
	return f'{days}:{hours}:{minutes}:{seconds}'


# Re-implemented from the Julia code
def color_configurations(N_fermions, num_sites, num_colors):
	if N_fermions == 0:
		return [0] * num_colors
	if N_fermions > num_colors * num_sites:
		return []
	if num_colors == 1:
		return [[N_fermions]]
	configurations = []
	for N_fermions_color1 in range(min(N_fermions, num_sites) + 1):
		new_configs = color_configurations(
			N_fermions-N_fermions_color1,
			min(num_sites, N_fermions_color1),
			num_colors - 1
		)
		new_configs = [[N_fermions_color1, *rest] for rest in new_configs]
		configurations.extend(new_configs)
	return configurations

# Assume matrix diagonalization scales w/ L^2
def calculate_scaling_factor(order, num_colors):
	matrix_size = 0
	N_max_fermions = order * num_colors
	for N_fermions in range(N_max_fermions + 1):
		for color_config in color_configurations(N_fermions, order, num_colors):
			L = math.prod(binomial(order, n) for n in color_config)
			matrix_size += L ** 2
	return matrix_size

#endregion

#region Initial Parsing

with open(args.clusterfile, 'r') as f:
	cluster_data = json.load(f)

cluster_file_absolute_path = os.path.abspath(args.clusterfile)

# Sorting by multiple conditions https://stackoverflow.com/a/14299498
# Sort first by the order number and then by the cluster id to stay in-sync with the main code
# sorted_cluster_ids = sorted(cluster_data.keys(), key=lambda id: (len(cluster_data[id][0]), int(id)))

sorted_cluster_lengths = sorted([len(cluster[0]) for cluster in cluster_data.values()])

job_name = args.jobname
mail_user = args.mail_user
no_mail = args.no_mail

if re.match(r'^\w+$', job_name) is None:
	print("Error: Job name can only contain letters, numbers, and underscores!")
	exit(1)

if mail_user is None and not no_mail:
	"Please enter an email address to receive notifications, use 'none' to disable notifications or <Enter> to use the default for your account."
	mail_user = [input_or_default("Email Address", f'{os.getlogin()}@g.hmc.edu')]
	if mail_user[0].strip().lower() == 'none':
		no_mail = True

general_params = {
	'nlce_home': NLCE_HOME,
	'project_root': PROJECT_ROOT,
	'mail_user': f'SBATCH --mail-user={",".join(mail_user)}' if mail_user else '<no mail-user>',
	'mail_type': 'NONE' if no_mail else 'END,FAIL',
}

#endregion

#region Test Mode

if args.test is not None:
	test_order = args.test
	test_cluster_idx = sorted_cluster_lengths.index(test_order)

	test_cluster_idx += 1 # Cluster IDs are 1-indexed in Julia but 0-indexed in Python
	job_params = {
		'job_name': f'NLCE_{job_name}',
		'log_file': f'{NLCE_HOME}/logs/{job_name}_%j.out',
		'time_limit': input('Time Limit: '),
		'memory_limit': input('Memory Limit: '),
		'cpus_per_task': input('CPUs per Task: '),
		'array_info': '<no array>',
		'command': f'{julia_base_command} -o "{NLCE_HOME}/output/{job_name}" diagonalize {cluster_file_absolute_path} {test_cluster_idx}',
		**general_params,
	}

	generate_script_from_template(f'{NLCE_HOME}/slurm/{job_name}.sh', job_params)
	exit(0)

#endregion

cpus = int(input('Num Cores Per Task: '))
num_colors = int(input_or_default('Num Colors', '3'))
mem_scaling = float(input('Memory Scaling Factor (GiB / el)'))
time_scaling = input_time('Time Scaling Factor: ')
cpu_max = int(input('Max Total Num CPUs: '))
mem_max = int(input('Max Total Memory: '))
min_batched_order = int(input_or_default('Min Batched Order', '4'))

max_order = sorted_cluster_lengths[-1]
max_order = int(input_or_default('Max NLCE Order', str(max_order)))

if max_order < min_batched_order or max_order not in sorted_cluster_lengths:
	raise ValueError('Max Order is less than Min Batched Order or not a valid cluster order!')

order_scales = {
	order: calculate_scaling_factor(order, num_colors)
	for order in range(min_batched_order, max_order+1)
}

if (max_time := time_scaling * order_scales[max_order]) > 24 * 60 * 60:
	raise ValueError(f'Max time for highest order would be {format_time(max_time)} > 2 days!')

batches = {
	order: {
		'time_limit': format_time(time_scaling * scale),
		'memory_limit': f'{mem_scaling * scale}gb',
		'cpus_per_task': cpus,
		'max_concurrent_tasks': min(cpu_max // cpus, mem_max // (mem_scaling * scale))
	} for order, scale in order_scales.items()
}

next_cluster_id = 1
batch_idx = 0
for order in range(min_batched_order, max_order+1):
	batch_params = batches[order]

	if order == min_batched_order:
		if min_batched_order == max_order:
			num_clusters_in_batch = len(sorted_cluster_lengths)
		else:
			num_clusters_in_batch = sorted_cluster_lengths.index(order + 1)
	else:
		num_clusters_in_batch = sorted_cluster_lengths.count(order)

	if num_clusters_in_batch == 0:
		continue

	array_range_start = next_cluster_id
	array_range_end = array_range_start + num_clusters_in_batch - 1
	next_cluster_id += num_clusters_in_batch

	job_params = {
		'job_name': f'NLCE_{job_name}_batch_{batch_idx}',
		'log_file': f'{NLCE_HOME}/logs/{job_name}_batch_{batch_idx}_cluster_%a_%j.out',
		'array_info': f'SBATCH --array={array_range_start}-{array_range_end}%{batch_params['max_concurrent_tasks']}',
		'command': f'{julia_base_command} -o '
			f'"{NLCE_HOME}/output/{job_name}_cluster_%a" '
			f'diagonalize {cluster_file_absolute_path} %a',
		**batch_params,
		**general_params
	}

	generate_script_from_template(f'{NLCE_HOME}/slurm/{job_name}_batch_{batch_idx}.slurm', job_params)

	batch_idx += 1

output_dirs = [f'"{NLCE_HOME}/output/{job_name}_cluster_{batch_num}"' for batch_num in range(batch_idx)]

merge_job_params = {
	'time_limit': '20:00',
	'memory_limit': '4gb',
	'cpus_per_task': 1,  # Merging does not benefit from multithreading
	'log_file': f'{NLCE_HOME}/logs/{job_name}_merge_%j.out',
	'array_info': '<no array>',
	'command': f'{julia_base_command} -o '
		f'"{NLCE_HOME}/output/{job_name}_merged '
		f'diagonalize {cluster_file_absolute_path} {' '.join(output_dirs)}',
}

generate_script_from_template(f'{NLCE_HOME}/slurm/{job_name}_merge.slurm', merge_job_params)

job_run_script = ['#!/bin/bash\n\n']

for batch_num in range(batch_idx):
	dependency = '' if batch_num == 0 else f'--dependency=afterok:$batch_{batch_num-1}_jobid'
	job_run_script.append(f'batch_{batch_num}_jobid=$(sbatch {dependency} {NLCE_HOME}/slurm/{job_name}_batch_{batch_num}.slurm)\n')

job_run_script.append(f'sbatch --dependency=afterok:$batch_{batch_idx-1}_jobid {NLCE_HOME}/slurm/{job_name}_merge.slurm\n')

with open(f'{NLCE_HOME}/slurm/{job_name}.sh', 'w') as f:
	f.writelines(job_run_script)
