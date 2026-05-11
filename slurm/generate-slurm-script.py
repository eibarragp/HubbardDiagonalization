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
cli_parser.add_argument('resourcefile', type=str, help="Path to the resource info file")
cli_parser.add_argument('jobname', type=str, help="Name of the SLURM job")
cli_parser.add_argument('--test', metavar='order', type=int, help="If specified, create a test script that runs one cluster of the specified order.")
cli_parser.add_argument('--mail-user', type=str, action='append', help="Email address to send notifications to. May be specified multiple times for multiple addresses.")
cli_parser.add_argument('--no-mail', action='store_true', help="If specified, do not send email notifications.")

args = cli_parser.parse_args()

heap_size_default_scaling = 0.7  # Default to setting Julia's heap size hint to x% of the requested memory

#endregion

#region Helpful Functions

def binomial(n, k):
	return math.factorial(n) / (math.factorial(k) * (math.factorial(n - k)))

def input_or_default(input_value, default_value):
	value = input(f'{input_value} [{default_value}]: ')
	return default_value if len(value.strip()) == 0 else value

def generate_script_from_template(filename, job_params):
	template_file = os.path.join(NLCE_HOME, f'{PROJECT_ROOT}/slurm/template.slurm')
	with open(template_file, 'r') as f:
		template = f.read()
	for key, value in job_params.items():
		# Use a lambda to ensure that the replacement value is treated as a literal string
		template = re.sub(rf'{{{{\s*{key}\s*}}}}', lambda _: str(value), template)

	with open(filename, 'w') as f:
		f.write(template)

def parse_range_or_int(value):
	if '-' in value:
		start, end = map(int, value.split('-'))
		return range(start, end + 1)
	else:
		return [int(value)]

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

	generate_script_from_template(f'{NLCE_HOME}/slurm/{job_name}.slurm', job_params)
	exit(0)

#endregion

#region Batch Creation

# Load resource allocation info
with open(args.resourcefile, 'r') as f:
	resource_data = json.load(f)

max_num_cpus = resource_data['max_cpus']
max_memory_gib = resource_data['max_mem_gb']

batches = [
	{
		'orders': parse_range_or_int(orders),
		**batchinfo
	} for orders, batchinfo in resource_data.items() if orders[0].isdigit()
]

batches_for_order = {}
for batch in batches:
	for order in batch['orders']:
		if order in batches_for_order:
			print(f"Error: Order {order} is included in multiple batches!")
			exit(1)
		batches_for_order[order] = batch

for order in batches_for_order.keys():
	if order > 1 and order - 1 not in batches_for_order:
		print(f"Error: Resource file specifies a discontinuous range of orders!")
		exit(1)

#endregion

#region SLURM Script Generation

merged_cpus = 0
merged_memory = 0
mergeable_batches = []

next_cluster_id = 1
for batch_idx, batch in enumerate(sorted(batches, key=lambda b: min(b['orders']))):
	relevant_clusters = list(filter(lambda id: len(cluster_data[id][0]) in batch['orders'], cluster_data.keys()))
	num_clusters_in_batch = len(relevant_clusters)

	array_range_start = next_cluster_id
	array_range_end = array_range_start + num_clusters_in_batch - 1
	next_cluster_id += num_clusters_in_batch

	ncpus = batch['ncpus']
	mem_gb = batch['mem_gb']
	num_julia_threads = batch.get('julia_threads', ncpus)
	# Fall back on Julia's defaults if unspecified
	num_mark_threads = batch.get('mark_threads', num_julia_threads / 2)
	num_sweep_threads = batch.get('sweep_threads', 0)

	num_mark_threads = max(1, math.ceil(num_mark_threads))

	max_concurrent_tasks = min(max_num_cpus // ncpus, max_memory_gib // mem_gb, num_clusters_in_batch)
	max_concurrent_tasks = max(1, max_concurrent_tasks)

	if max_concurrent_tasks == num_clusters_in_batch and \
		(batch_idx == 0 or batch_idx-1 in mergeable_batches) and \
		merged_cpus + (ncpus * max_concurrent_tasks) <= max_num_cpus and \
		merged_memory + (mem_gb * max_concurrent_tasks) <= max_memory_gib:
		# This batch can be merged with the previous batch
		merged_cpus += ncpus * max_concurrent_tasks
		merged_memory += mem_gb * max_concurrent_tasks
		mergeable_batches.append(batch_idx)

	job_params = {
		'job_name': f'NLCE_{job_name}_batch_{batch_idx}',
		'log_file': f'{NLCE_HOME}/logs/{job_name}_batch_{batch_idx}_cluster_%a_%j.out',
		'array_info': f'SBATCH --array={array_range_start}-{array_range_end}%{max_concurrent_tasks}',
		'command': f'{julia_base_command} -o '
			f'"{NLCE_HOME}/output/{job_name}_cluster_$SLURM_ARRAY_TASK_ID" '
			f'diagonalize {cluster_file_absolute_path} $SLURM_ARRAY_TASK_ID',
		'cpus_per_task': ncpus,
		'memory_limit': mem_gb,
		'time_limit': batch['time'],
		'num_julia_threads': num_julia_threads,
		'num_mark_threads': num_mark_threads,
		'num_sweep_threads': num_sweep_threads,
		'heap_size_hint': batch.get('heap_size_hint', round(heap_size_default_scaling * mem_gb, 3)),
		**general_params
	}

	generate_script_from_template(f'{NLCE_HOME}/slurm/{job_name}_batch_{batch_idx}.slurm', job_params)

output_dirs = [f'"{NLCE_HOME}/output/{job_name}_cluster_{cluster_id}"' for cluster_id in range(1, next_cluster_id)]

if "merge" not in resource_data:
	resource_data["merge"] = {
		'time': '20:00',
		'ncpus': 1,
		'mem_gb': 4
	}

if resource_data["merge"]['ncpus'] > max_num_cpus or resource_data["merge"]['mem_gb'] > max_memory_gib:
	print("Warning: Merge job resource requirements exceed set maximums!")

num_merge_julia_threads = resource_data["merge"].get('julia_threads', resource_data["merge"]['ncpus'])
num_merge_mark_threads = resource_data["merge"].get('mark_threads', num_merge_julia_threads // 2)
num_merge_sweep_threads = resource_data["merge"].get('sweep_threads', 0)

num_merge_mark_threads = max(1, math.ceil(num_merge_mark_threads))

merge_job_params = {
	'job_name': f'NLCE_{job_name}_merge',
	'log_file': f'{NLCE_HOME}/logs/{job_name}_merge_%j.out',
	'time_limit': resource_data['merge']['time'],
	'memory_limit': resource_data['merge']['mem_gb'],
	'cpus_per_task': resource_data['merge']['ncpus'],
	'array_info': '<no array>',
	'command': f'{julia_base_command} -o '
		f'"{NLCE_HOME}/output/{job_name}_merged" '
		f'merge {cluster_file_absolute_path} {" ".join(output_dirs)}',
	'num_julia_threads': num_merge_julia_threads,
	'num_mark_threads': num_merge_mark_threads,
	'num_sweep_threads': num_merge_sweep_threads,
	'heap_size_hint': resource_data['merge'].get('heap_size_hint', round(heap_size_default_scaling * resource_data['merge']['mem_gb'], 3)),
	**general_params
}

generate_script_from_template(f'{NLCE_HOME}/slurm/{job_name}_merge.slurm', merge_job_params)

#endregion

#region Job Run Script Generation

job_run_script = ['#!/bin/bash\n\n']

for batch_num in range(batch_idx+1):
	if batch_num == 0 or batch_num in mergeable_batches:
		dependency = ''
	elif batch_num - 1 in mergeable_batches:
		dependency = ' '.join(
			f'--dependency=afterok:$batch_{merged_batch_num}_jobid' for merged_batch_num in mergeable_batches
		)
	else:
		dependency = f'--dependency=afterok:$batch_{batch_num-1}_jobid'

	job_run_script.append(f'batch_{batch_num}_jobid=$(sbatch --parsable {dependency} --kill-on-invalid-dep=yes {NLCE_HOME}/slurm/{job_name}_batch_{batch_num}.slurm)\n')

job_run_script.append(f'sbatch --dependency=afterok:$batch_{batch_idx}_jobid --kill-on-invalid-dep=yes {NLCE_HOME}/slurm/{job_name}_merge.slurm\n')

# Print the queued jobs
job_run_script.append('squeue -u $USER\n')

with open(f'{NLCE_HOME}/slurm/{job_name}.sh', 'w') as f:
	f.writelines(job_run_script)

os.chmod(f'{NLCE_HOME}/slurm/{job_name}.sh', 0o755)

#endregion
