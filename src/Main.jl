# Include Libraries
using HubbardDiagonalization:
	CSVUtil,
	DataHelpers,
	Graphs, Graphs.Graph,
	ExactDiagonalization as ED

import CSV
import JSON3 as JSON
import Logging
import TOML

using ArgParse
using Base.Threads

function parse_cli()
	s = ArgParseSettings()

	@add_arg_table s begin
		"--loglevel", "-l"
			help = "Set the logging level (debug, info, warn, error, etc.)"
			arg_type = Logging.LogLevel
			default = Logging.Info
		"--outputdir", "-o"
			help = "Directory to save output files to"
			arg_type = String
			default = "output"
		"diagonalize"
			help = "Diagonalize and compute observables for the specified cluster"
			action = :command
		"merge"
			help = "Merge the results of multiple diagonalization runs using NLCE"
			action = :command
		"simple"
			help = "Diagonalize and compute observables and graphs for the graph provided in SimulationConfig.toml"
			action = :command
	end

	@add_arg_table s["diagonalize"] begin
		"--generate-plots"
			help = "By default, plot generation is skipped for cluster runs. Use this flag to reenable it."
			action = :store_true
		"clusterfile"
			help = "Path to the cluster info file."
			arg_type = String
		"clusteridx"
			help = "Index of the cluster from the cluster file to diagonalize."
			arg_type = Int
	end

	@add_arg_table s["merge"] begin
		"clusterfile"
			help = "Path to the cluster info file."
			arg_type = String
		"datadirs"
			help = "Directories containing the output of diagonalization runs to merge."
			arg_type = String
			nargs = '+'
	end

	return parse_args(s)
end

"""
    main(args::Vector{String}) -> Int

The main entry point! Handles the high-level control flow of the program.
"""
function (@main)(args)
    parsed_args = parse_cli()

	@info "Args: $parsed_args"

    # Install our own logger for the duration of the program
    old_logger = Logging.global_logger(Logging.ConsoleLogger(stdout, parsed_args["loglevel"]))

    if nthreads() == 1
        @warn "Running in single-threaded mode. For better performance, consider setting the JULIA_NUM_THREADS environment variable to a higher value."
    else
        @info "Running with $(nthreads()) threads."
    end

    # Parse configuration file
    config = TOML.parsefile("SimulationConfig.toml")

    params = config["parameters"]
    plot_config = config["plot"]
    graph_config = config["graph"]

    test_config = ED.TestConfiguration(; DataHelpers.convert_strings_to_symbols(params)...)
    t_vals = plot_config["T_min"]:plot_config["T_step"]:plot_config["T_max"]
    u_vals = plot_config["u_min"]:plot_config["u_step"]:plot_config["u_max"]

	if parsed_args["%COMMAND%"] == "diagonalize"
		cluster_file = parsed_args["diagonalize"]["clusterfile"]
		cluster_idx = parsed_args["diagonalize"]["clusteridx"]

		cluster_data = JSON.read(read(cluster_file, String))
		if cluster_idx < 1 || cluster_idx > length(keys(cluster_data))
			error("Cluster index $cluster_idx out of bounds! File contains $(length(cluster_data["clusters"])) clusters.")
		end

		cluster = DataHelpers.sorted_clusters(cluster_data)[cluster_idx]
		graph = Graphs.from_cluster(cluster)

		# Unless told otherwise, don't generate plots
		if !parsed_args["diagonalize"]["generate-plots"]
			plot_config["T_fixed_plots"] = []
			plot_config["u_fixed_plots"] = []
		end
	elseif parsed_args["%COMMAND%"] == "simple"
    	graph = Graphs.linear_chain(graph_config["num_sites"])
	end

	if @isdefined(graph)
		observables, derived_observables, overlays = ED.default_observables(test_config, graph)
		@info "Defined observables: $(union(keys(observables), keys(derived_observables), keys(overlays)))"

		observable_data = ED.diagonalize_and_compute_observables(
			t_vals,
			u_vals,
			test_config,
			graph,
			observables,
			derived_observables,
			overlays,
		)

		# If we're not using NLCE,
		if parsed_args["%COMMAND%"] == "simple"
			# Normalize energy and entropy per site
			observable_data["Energy"] ./= Graphs.num_sites(graph)
			observable_data["Entropy"] ./= Graphs.num_sites(graph)
		end
		# Otherwise, we must report the extrinsic energy to the merger,
		# and NLCE will take care of the normalization for us.
	elseif parsed_args["%COMMAND%"] == "merge"
		cluster_file = parsed_args["merge"]["clusterfile"]
		data_dirs = parsed_args["merge"]["datadirs"]

		observable_names = String[]
		for ((observables, _)) in plot_config["observables"]
			append!(observable_names, observables)
		end

		observable_data = DataHelpers.merge_results_with_nlce(cluster_file, data_dirs, observable_names, plot_config["nlce_orders"])
	else
		error("Invalid State! 'graph' is not defined and command is not 'merge'. This should never happen!")
	end

    DataHelpers.export_observable_data(plot_config, t_vals, u_vals, observable_data, test_config,
		@isdefined(graph) ? string(Graphs.num_sites(graph)) : "NLCE",
		parsed_args["%COMMAND%"] == "merge",
		parsed_args["outputdir"])

    @info "Done."

    # Restore the old logger after we're done
    Logging.global_logger(old_logger)

    return 0
end

function ArgParse.parse_item(::Type{Logging.LogLevel}, s::AbstractString)
	level = lowercase(s)
	if level == "debug"
		return Logging.Debug
	elseif level == "info"
		return Logging.Info
	elseif level == "warn" ||
			level == "warning"
		return Logging.Warn
	elseif level == "error"
		return Logging.Error
	elseif level == "min" ||
			level == "all" ||
			level == "belowminlevel"
		return Logging.BelowMinLevel
	elseif level == "max" ||
			level == "none" ||
			level == "abovemaxlevel"
		return Logging.AboveMaxLevel
	else
		try
			return Logging.LogLevel(parse(Int, s))
		catch
			error("Invalid logging level: $s")
		end
	end
end
