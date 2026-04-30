using HubbardDiagonalization:
	CSVUtil, DataHelpers, ExactDiagonalization as ED

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
		"plot_config"
		help = "Path to a TOML file containing plotting configuration options."
		arg_type = String
		default = "PlotConfig.toml"
		"datadirs"
        help = "Directories containing the output data to merge."
        arg_type = String
        nargs = '+'
    end

    parsed_args = parse_args(s)
    parsed_args["validate"] = lowercase(parsed_args["validate"])

    return parsed_args
end
struct DataMatrix
	u_vals::Vector{Float64}
	T_vals::Vector{Float64}
	values::Matrix{Float64}
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
		test_config, U = DataHelpers.validate_datasets(cli_args["datadirs"], cli_args["validate"] == "scan_only")
		@info "Validation complete."
	end

	config = TOML.parsefile(cli_args["plot_config"])

	validate_config(config, U)

	variables = Set{String}(["u", "U", "T"])
	variable_names = Dict(
		"u" => "Chemical Potential (μ)",
		"U" => "Interaction Strength (U)",
		"T" => "Temperature (T)"
	)
	for (fig_name, fig_config) in config
		@info "Generating plot for $fig_name..."
		observable_ids = fig_config["observables"]["ids"]
		observable_names = fig_config["observables"]["names"]

		for observable_id in observable_ids
			observable_name = observable_names[observable_id]

			for axis in variables
				axis_range = get!(fig_config, "$(axis)_range", [])
				if isempty(axis_range)
					continue
				end

				is_overlay, secondary_var, secondary_range = if haskey(fig_config, "overlay")
					var = fig_config["overlay"][1]

					if var == axis
						continue
					end

					range = load_range(fig_config["overlay"][2:4])
					(true, var, range)
				else
					var = axis == "U" ? "u" : "U"
					(false, var, load_range(fig_config["fixed_$(var)"]))
				end

				fixed_var = setdiff(variables, [axis, var])[1]
				fixed_value_range = load_range(fig_config["fixed_$(fixed_var)"])

				new_figure() = plot(
					xlabel = variable_names[axis],
					ylabel = observable_name,
					title = fig_name,
					legend_position = true,
					xscale = is_log(axis_range) ? :log10 : :identity
				)

				for v1 in fixed_value_range
					if is_overlay
						figure = new_figure()
					end

					for v2 in secondary_range
						if !is_overlay
							figure = new_figure()
						end

						for prefix_name in keys(prefixes)
							create_plots_for_resummation_method!(
								figure,
								fig_config,
								cli_args["datadirs"],
								Us,
								observable_id,
								prefix_name,
								axis, axis_range,
								fixed_var, v1,
								secondary_var, v2,
								is_overlay
							)
						end

						if !is_overlay
							savefig(figure, filename)
						end
					end

					if is_overlay
						savefig(figure, filename)
					end
				end
			end
		end
	end
end

function validate_config(config::Dict{String, Any}, Us::Vector{Float64})
	if haskey(config, "defaults")
		defaults = config["defaults"]

		for figure in keys(config)
			if figure == "defaults"
				continue
			end
			config[figure] = merge(defaults, config[figure])
		end
	end

	for (fig_name, fig_config) in config
		if !haskey(fig_config, "observables") || !haskey(fig_config["observables"], "ids")
			error("Figure $fig_name is missing required 'observables.ids' field.")
		end

		get!(fig_config["observables"], "names", Dict())
		for observable_id in fig_config["observables"]["ids"]
			get!(observable_names, observable_id, observable_id)
		end

		get!(fig_config, "U_range", Us)
		get!(fig_config, "fixed_U", Us)

		get!(fig_config, "atol", 1e-3)
		get!(fig_config, "rtol", 0.05)

		if isempty(get!(fig_config, "prefixes", []))
			error("Figure $fig_name is missing required 'prefixes' field.")
		end

		for (prefix_name, prefix_data) in prefixes
			get!(prefix_data, "min_order", 1)
			if !haskey(prefix_data, "orders")
				error("Figure $fig_name is missing required 'orders' field in prefix $prefix_name.")
			end
		end
	end
end

function load_data_matrix(file::String)
	data_matrix = CSVUtil.load_csv_matrix(file)

	return DataMatrix(
		data_matrix[1, 2:end],
		data_matrix[2:end, 1],
		data_matrix[2:end, 2:end],
	)
end

function load_data_matrix(datadirs::Vector{String}, Us::Vector{Float64}, U::Float64, observable_id::String, prefix::String, order::Int)
	dir_idx = findfirst(==(U), Us)
	if dir_idx === nothing
		error("No data found for U = $U")
	end
	return load_data_matrix(joinpath(datadirs[dir_idx], "$(observable_id)_$(prefix)_$(order).csv"))
end

load_range(range) = Float64(range[1]):Float64(range[2]):Float64(range[3])
is_log(range) = length(range) == 4 && range[4] == "log"

function get_data_for_params(
	datadirs::Vector{String},
	Us::Vector{Float64},
	observable_id::String,
	prefix::String,
	order::Int,
	:U;
	T::Float64,
	u::Float64,
	U::Vector{Float64}
)
	valid_Us = Int[]
	data_series = Float64[]
	for U_val in U
		data_mat = load_data_matrix(datadirs, Us, U_val, observable_id, prefix, order)
		T_idx = findfirst(==(T), data_mat.T_vals)
		if T_idx === nothing
			@warn "No data found for U = $U_val and T = $T"
		end
		u_idx = findfirst(==(u), data_mat.u_vals)
		if u_idx === nothing
			@warn "No data found for U = $U_val and μ = $u"
		end
		data_series_val = data_mat.values[T_idx, u_idx]

		push!(valid_Us, U_val)
		push!(data_series, data_series_val)
	end

	return valid_Us, data_series
end

function get_data_for_params(
	datadirs::Vector{String},
	Us::Vector{Float64},
	observable_id::String,
	prefix::String,
	order::Int,
	:T;
	T::Vector{Float64},
	u::Float64,
	U::Float64
)
	data_mat = load_data_matrix(datadirs, Us, U, observable_id, prefix, order)
	idx = findfirst(==(u), data_mat.u_vals)
	if idx === nothing
		error("No data found for μ = $u")
	end
	axis_range = data_mat.T_vals
	data_series = data_mat.values[:, idx]

	min_idx = findfirst(>=(T[1]), axis_range)
	max_idx = findlast(<=(T[end]), axis_range)

	return axis_range[min_idx:max_idx], data_series[min_idx:max_idx]
end

function get_data_for_params(
	datadirs::Vector{String},
	Us::Vector{Float64},
	observable_id::String,
	prefix::String,
	order::Int,
	:u;
	T::Float64,
	u::Vector{Float64},
	U::Float64
)
	data_mat = load_data_matrix(datadirs, Us, U, observable_id, prefix, order)
	idx = findfirst(==(T), data_mat.T_vals)
	if idx === nothing
		error("No data found for T = $T")
	end
	axis_range = data_mat.u_vals
	data_series = data_mat.values[idx, :]

	min_idx = findfirst(>=(u[1]), axis_range)
	max_idx = findlast(<=(u[end]), axis_range)

	return axis_range[min_idx:max_idx], data_series[min_idx:max_idx]
end

function create_plots_for_resummation_method!(
	figure,
	fig_config::Dict{String, Any},
	datadirs::Vector{String},
	Us::Vector{Float64},
	observable_id::String,
	prefix::String,
	axis::String,
	axis_range::Vector{Float64},
	fixed_var::String,
	fixed_value::Float64,
	secondary_var::String,
	secondary_value::Float64,
	is_overlay::Bool
)
	atol = fig_config["atol"]
	rtol = fig_config["rtol"]

	datasets = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()

	# Also load data for the next order to enable convergence testing.
	orders_to_load = axis == "T" ? union(orders, orders .+ 1) : orders

	for order in orders_to_load
		try
			datasets[order] = get_data_for_params(
				datadirs,
				Us,
				observable_id,
				prefix,
				order,
				Symbol(axis);
				Symbol(axis) => axis_range,
				Symbol(fixed_var) => fixed_value,
				Symbol(secondary_var) => secondary_value,
			)
		catch
			# Don't raise an error if we fail to load an order that wasn't explicitly requested
			if order in orders
				rethrow()
			end
		end
	end

	cutoff_idxs = Dict{Int, Int}()
	if axis == "T"
		# Perform convergence tests
		for order in orders
			# Find an order to compare against
			# Check, in order:
			# 1) The next order
			# 2) The next largest order that we have data for
			# 3) The next smallest order that we have data for
			# 4) Give up and don't use a cutoff
			cmp_order = order + 1
			if !haskey(datasets, cmp_order)
				cmp_order = findmax(filter(<(order), keys(datasets)))
			end
			if !haskey(datasets, cmp_order)
				cmp_order = findmin(filter(>(order), keys(datasets)))
			end
			if !haskey(datasets, cmp_order)
				cutoff_idxs[order] = 1
				continue
			end

			values = datasets[order][2]
			cmp_values = datasets[cmp_order][2]

			# By high temperature, everything should be converged
			converged_value = values[end]

			# All orders (with the same params) should be generated from the same temperatures,
			# so the relevant indicies should be identical.
			cutoff = findlast(zip(values, cmp_values)) do (value, cmp_value)
				# Find the last non-converged point to avoid cases where the data
				# looks converged but re-diverges
				# Don't bother with adding 1 to get to the first converged point because
				# the step size *should* be small.
				!(
					# Converged if:
					# 1) Difference is within tolerance, or
					abs(value - cmp_value) <= atol ||
					# 2) Difference is small relative to distance from the converged value
					abs(value - cmp_value) <= rtol * abs(value - converged_value)
				)
			end

			if cutoff === nothing
				@warn "$observable_id with $prefix order $order failed to converge!"
				cutoff_idxs[order] = minimum(temps)
			else
				cutoff_idxs[order] = cutoff
			end
		end
	end

	for order in orders
		temps, values = datasets[order]
		cutoff = cutoffs[order]
		overlay_var_label = is_overlay ? "$secondary_var = $secondary_value, " : ""
		plot!(
			figure,
			temps[cutoff:end],
			values[cutoff:end],
			label = "$(overlay_var_label)$prefix Order $order"
		)
	end
end
