#region Setup

import argparse
import json
import os
import re

NLCE_HOME="."

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

def input_or_default(input_value, default_value):
	value = input(f'{input_value} [{default_value}]: ')
	return default_value if len(value.strip()) == 0 else value

def generate_script_from_template(filename, job_params):
	template_file = os.path.join(NLCE_HOME, f'{NLCE_HOME}/HubbardDiagonalization/slurm/template.sh')
	with open(template_file, 'r') as f:
		template = f.read()
	for key, value in job_params.items():
		# Use a lambda to ensure that the replacement value is treated as a literal string
		template = re.sub(rf'{{{{\s*{key}\s*}}}}', lambda x: value, template)

	with open(filename, 'w') as f:
		f.write(template)

#endregion

#region Initial Parsing

with open(args.clusterfile, 'r') as f:
	cluster_data = json.load(f)

cluster_file_absolute_path = os.path.abspath(args.clusterfile)

# Sorting by multiple conditions https://stackoverflow.com/a/14299498
# Sort first by the order number and then by the cluster id to stay in-sync with the main code
sorted_cluster_ids = sorted(cluster_data.keys(), key=lambda id: (len(cluster_data[id][0]), int(id)))

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

mail_params = {
	'mail_user': f'SBATCH --mail-user={",".join(mail_user)}' if mail_user else '<no mail-user>',
	'mail_type': 'NONE' if no_mail else 'END,FAIL',
}

#endregion

#region Test Mode

if args.test is not None:
	test_order = args.test
	test_cluster_idx = None
	for (i, cluster_id) in enumerate(sorted_cluster_ids):
		if len(cluster_data[cluster_id][0]) == test_order:
			test_cluster_idx = i
			break

	if test_cluster_idx is None:
		print(f"Error: No cluster of order {test_order} found in the cluster file!")
		exit(1)

	test_cluster_idx += 1 # Cluster IDs are 1-indexed in Julia but 0-indexed in Python
	job_params = {
		'job_name': f'NLCE_{job_name}',
		'log_file': f'{NLCE_HOME}/logs/{job_name}_%j.out',
		'time_limit': input('Time Limit: '),
		'memory_limit': input('Memory Limit: '),
		'cpus_per_task': input('CPUs per Task: '),
		'dependencies': '<no dependencies>',
		'array_info': '<no array>',
		'command': f'{julia_base_command} -o "{NLCE_HOME}/output/{job_name}" diagonalize {cluster_file_absolute_path} {test_cluster_idx}',
		**mail_params,
	}

	generate_script_from_template(f'{job_name}.sh', job_params)
	exit(0)

#endregion

print("Not implemented yet!")
