using HubbardDiagonalization:
	ExactDiagonalization as ED

import TOML

using ArgParse
using Plots

function parse_cli()
	s = ArgParseSettings()

	@add_arg_table s begin
		# "--loglevel", "-l"
		# help = "Set the logging level (debug, info, warn, error, etc.)"
		# arg_type = Logging.LogLevel
		# default = Logging.Info
		"--outputdir", "-o"
		help = "Directory to save output files to"
		arg_type = String
		default = "output"
		"--validate"
		help = "How to validate input data before processing. (Options: yes, no, scan_only)"
		arg_type = String
		default = "yes"
		range_tester = x -> lowercase(x) in ("yes", "no", "scan_only")
		"datadirs"
        help = "Directories containing the output data to merge."
        arg_type = String
        nargs = '+'
	end

	parsed_args = parse_args(s)
	parsed_args["validate"] = lowercase(parsed_args["validate"])

	return parsed_args
end

function (@main)(args)
	cli_args = parse_cli()

	if cli_args["validate"] == "no"
		if length(cli_args["datadirs"]) > 1
			error("Cannot merge data from multiple directories without validation.")
		end
		@warn "Skipping validation of input data."

    	test_config = ED.TestConfiguration(num_colors = 0, t = 0, u_test = 0, U = 0)
		U = [0]
	else
		@info "Performing initial validation of input data..."
		test_config, U = validate_datasets(cli_args["datadirs"], cli_args["validate"] == "scan_only")
		@info "Validation complete."
	end
end

function load_config_from_output_dir(datadir::String)
	config_path = joinpath(datadir, "SimulationConfig.toml")
	if !isfile(config_path)
		error("SimulationConfig.toml not found in $datadir")
	end
	params = TOML.parsefile(config_path).get("parameters", nothing)
	if params === nothing
		error("No [parameters] section found in SimulationConfig.toml at $config_path")
	end
	config = ED.TestConfiguration(; Utils.convert_strings_to_symbols(params))
	config.u_test = 0  # We don't care about u_test
	return config
end

function validate_datasets(datadirs::Vector{String}, scan_only::Bool)
	test_config = load_config_from_output_dir(datadirs[1])
	Us = [test_config.U]

	for datadir in @view datadirs[2:end]
		config = load_config_from_output_dir(datadir)
		push!(Us, config.U)

		if !scan_only && config != test_config
			error("Test configuration mismatch between $datadir and $(datadirs[1]).")
		end
	end

	return test_config, Us
end
