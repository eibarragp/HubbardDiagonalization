using HubbardDiagonalization: CSVUtil, DataHelpers, ExactDiagonalization as ED, Utils

import Measures
import TOML

using ArgParse
using Base: Fix1, Fix2
using ColorSchemes
using Plots

const variables = Set{String}(["u", "U", "T"])
const variable_names = Dict(
    "u" => "Chemical Potential (μ)",
    "U" => "Interaction Strength (U)",
    "T" => "Temperature (T)",
)

"""
	parse_cli() -> Dict{String, Any}

Parses the command line arguments and returns a dictionary of the parsed values.
"""
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
        "--figure", "-f"
        help =
            "Only generate plots for the specified figure(s) from the config file." *
            " If passed multiple times, will generate plots for all specified figures." *
            " By default, plots for all figures are generated."
        arg_type = String
        action = :append_arg
        "plot_config"
        help = "Path to a TOML file containing plotting configuration options."
        arg_type = String
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

"""
    @main(args)

Entrypoint for the plotting script!
"""
function (@main)(args)
    # Parse CLI arguments
    cli_args = parse_cli()
    @info "Args: $cli_args"

    # Validate input datasets and extract relevant parameters
    if cli_args["validate"] == "no"
        # We currently do not have any other way to determine the U values,
        # so in this case, we can't distinguish between different datasets
        if length(cli_args["datadirs"]) > 1
            error("Cannot merge data from multiple directories without validation.")
        end
        @warn "Skipping validation of input data."

        # We still don't know U, but at least there's only one dataset, so it doesn't really matter
        # Assume its 0
        test_config = ED.TestConfiguration(num_colors = 0, t = 1, u_test = 0, U = 0)
        Us = [0.0]
    else
        @info "Performing initial validation of input data..."
        test_config, Us = DataHelpers.validate_datasets(
            cli_args["datadirs"],
            cli_args["validate"] == "scan_only",
            false,
        )
        @info "Validation complete."
    end

    # Load the plotting configuration file
    config = TOML.parsefile(cli_args["plot_config"])
    validate_config(config, Us)

    # Create a directory for the output
    output_dir = cli_args["outputdir"]
    Base.Filesystem.mkpath(output_dir)

    # Now, generate the requested plots!

    # For every plot
    for (fig_name, fig_config) in config
        # Besides the reserved "defaults" config
        if fig_name == "defaults"
            continue
        end

        # If we were told to only generate specific figures, check if this is one of them
        figure_requested = isempty(cli_args["figure"]) || fig_name in cli_args["figure"]
        if !figure_requested
            continue
        end

        @info "Generating plots for $fig_name..."
        observable_ids = fig_config["observables"]["ids"]
        observable_names = fig_config["observables"]["names"]

        # Precompute the number of data-series we need for each configuration
        # of parameters to properly size the color palette.
        # This is just the total number of orders we're plotting
        num_series_per_param_set =
            sum(values(fig_config["prefixes"]) .|> Fix2(getindex, "orders") .|> length)

        # Put each observable on its own plot
        for observable_id in observable_ids
            observable_name = observable_names[observable_id]

            @info "Plotting observable $observable_id ($observable_name)..."

            # Check possibilities for the x-axis
            for axis in variables
                axis_range = get!(fig_config, "$(axis)_range", [])
                if isempty(axis_range)
                    # User did not tell us to use this axis. Skip it!
                    continue
                end

                # Select one of the other variables to use a a fixed axis
                is_overlay, secondary_var, secondary_range =
                    if haskey(fig_config, "overlay")
                        # The user wants us to overlay multiple values!
                        # Find out which ones!
                        var = fig_config["overlay"][1]

                        if var == axis
                            # Can't overlay multiple values if they're all already in the domain
                            continue
                        end

                        range = load_fixed_range(fig_config["overlay"][2:end])
                        (true, var, range)
                    else
                        # No explicit overlay specified. Just choose one that's not the x-axis
                        var = axis == "U" ? "u" : "U"
                        (false, var, load_fixed_range(fig_config["fixed_$(var)"]))
                    end

                # Now, select the remaining variable to be the last fixed value
                fixed_var = setdiff(variables, [axis, var]) |> only
                fixed_value_range = load_fixed_range(fig_config["fixed_$(fixed_var)"])

                is_logarithmic = is_log(axis_range)
                axis_range = axis_range[1:2] .|> Float64

                # For every value of the fixed variable, make a plot
                for v1 in fixed_value_range
                    # If we're overlaying multiple values for the secondary variable,
                    # we have to initialize the plot here so that they can all go into the same figure.
                    if is_overlay
                        figure, filename = init_figure(
                            fig_config,
                            # We'll be putting length(secondary_range) number of different
                            # plots on here. Each one will have num_series_per_param_set prefix-order pairs
                            length(secondary_range) * num_series_per_param_set,
                            test_config,
                            observable_id,
                            observable_name,
                            is_logarithmic,
                            axis,
                            axis_range,
                            fixed_var,
                            v1,
                            secondary_var,
                            nothing,
                        )
                    end

                    # Now, plot every possible value of the secondary variable
                    for v2 in secondary_range
                        # If we're overlaying, we already initialized the figure
                        # Otherwise, make one for this specific value of the secondary variable
                        if !is_overlay
                            figure, filename = init_figure(
                                fig_config,
                                # In ths case, we're only using one set of parameters,
                                # so each prefix-order pair only shows up once.
                                num_series_per_param_set,
                                test_config,
                                observable_id,
                                observable_name,
                                is_logarithmic,
                                axis,
                                axis_range,
                                fixed_var,
                                v1,
                                secondary_var,
                                v2,
                            )
                        end

                        # Now that we've decided on parameters, plot every requested order onto the graph.
                        for prefix_name in keys(fig_config["prefixes"])
                            create_plots_for_resummation_method!(
                                figure,
                                fig_config,
                                cli_args["datadirs"],
                                test_config,
                                Us,
                                observable_id,
                                prefix_name,
                                fig_config["prefixes"][prefix_name]["orders"],
                                axis,
                                axis_range,
                                fixed_var,
                                v1,
                                secondary_var,
                                v2,
                                is_overlay,
                            )
                        end

                        # If we're done with the figure, save it!
                        if !is_overlay
                            savefig(figure, joinpath(output_dir, filename))
                        end
                    end

                    # If we're overlaying, the figure is not finished until here,
                    # when all the different values of the secondary variable have been plotted.
                    if is_overlay
                        savefig(figure, joinpath(output_dir, filename))
                    end
                end
            end
        end
    end
end

"""
    validate_config(config::Dict{String,Any}, Us::Vector{Float64})

Validates the plotting configuration file and fills in any missing default values.
config: The plotting configuration dictionary loaded from the TOML file.
Us: The list of U values that we have data for. This is used as the default value for the "fixed_U" parameter.
"""
function validate_config(config::Dict{String,Any}, Us::Vector{Float64})
    # Defaults to be merged with the user-specified configs
    default_params = Dict(
        "observables" => Dict("ids" => [], "names" => Dict()),
        "fixed_U" => Us,
        "atol" => 1e-3,
        "rtol" => 0.05,
        "prefixes" => Dict(),
        "color_palette" => "coolwarm",
        "name_num_round_digits" => 3,
        "plot_params" => Dict("minorgrid" => true),
        "series_params" => Dict(),
    )

    # Merge our_defaults -> user_defaults -> figure_configs
    defaults = get(config, "defaults", Dict())
    defaults = Utils.recursive_merge(default_params, defaults)
    config["defaults"] = defaults

    for figure in keys(config)
        if figure == "defaults"
            continue
        end
        config[figure] = Utils.recursive_merge(defaults, config[figure])
    end

    # Plots.jl sometimes requires non-encodable types for certain parameters
    # define a standard way of making the relevant conversions
    type_parsers = Dict(Measures.Length => Fix1(Measures.Length, :mm), Symbol => Symbol)
    typed_plot_fields = Dict("margin" => Measures.Length)
    typed_series_fields = Dict("markershape" => Symbol)

    # Perform validation on each config
    for (fig_name, fig_config) in config
        # Set the display name for each observable to its id if none was specified
        for observable_id in fig_config["observables"]["ids"]
            get!(fig_config["observables"]["names"], observable_id, observable_id)
        end

        # Each prefix must have the `orders`` field defined
        for (prefix_name, prefix_data) in fig_config["prefixes"]
            if !haskey(prefix_data, "orders")
                error(
                    "Figure $fig_name is missing required 'orders' field in prefix $prefix_name.",
                )
            end
        end

        # Perform type conversions in plot_params and series_params
        plot_params = fig_config["plot_params"]
        for (field, field_type) in typed_plot_fields
            if haskey(plot_params, field) && !(plot_params[field] isa field_type)
                plot_params[field] = type_parsers[field_type](plot_params[field])
            end
        end
        series_params = fig_config["series_params"]
        for (field, field_type) in typed_series_fields
            if haskey(series_params, field) && !(series_params[field] isa field_type)
                series_params[field] = type_parsers[field_type](series_params[field])
            end
        end
    end
end

"""
    format_name(name, ...)

A wrapper around Utils.format_string that defines some standard variable types and performs some basic formatting/normalization
name: The name string to format
The remaining parameters are used to populate the variables that can be used in the name string.
They have the same meanings as in (@main)()
"""
function format_name(
    name::String,
    test_config::ED.TestConfiguration,
    observable_id::String,
    observable_name::String,
    is_logarithmic::Bool,
    axis::String,
    axis_range::Vector{Float64},
    fixed_var::String,
    fixed_value::Float64,
    secondary_var::String,
    secondary_value::Union{Float64,Nothing},
    name_num_round_digits::Int,
)
    # Round numeric values to avoid having unwieldy numbers in the formatted strings
    axis_range = round.(axis_range, digits = name_num_round_digits)
    fixed_value = round(fixed_value, digits = name_num_round_digits)
    secondary_value =
        isnothing(secondary_value) ? "overlay" :
        round(secondary_value, digits = name_num_round_digits)

    # Define the valid variable replacements
    available_variables = Dict(
        "oi" => observable_id,
        "on" => observable_name,
        "xn" => axis,
        "log" => is_logarithmic ? "logarithmic" : "linear",
        "xmin" => axis_range[1],
        "xmax" => axis_range[2],
        "f1n" => fixed_var,
        "f1v" => fixed_value,
        "f2n" => secondary_var,
        "f2v" => secondary_value,
        "t" => test_config.t,
        "ut" => test_config.u_test,
        "nc" => test_config.num_colors,
        # Ensures that the "T", "u", and "U" variables are mapped to their respective values regardless of which one is used where
        axis => axis_range,
        fixed_var => fixed_value,
        secondary_var => secondary_value,
    )
    return Utils.format_string(name, available_variables)
end

"""
    init_figure(fig_config, num_series, ...) -> (figure, filename)

Initializes a Plots.jl figure based on the provided configuration and parameters.
fig_config: The user-specified configuration for the figure
num_series: The number of data series that will be plotted on this figure. This is used to properly size the color palette.
The rest of the parameters have the same meanings as in (@main)()
"""
function init_figure(
    fig_config::Dict{String,Any},
    num_series::Int,
    test_config::ED.TestConfiguration,
    observable_id::String,
    observable_name::String,
    is_logarithmic::Bool,
    axis::String,
    axis_range::Vector{Float64},
    fixed_var::String,
    fixed_value::Float64,
    secondary_var::String,
    secondary_value::Union{Float64,Nothing},
)
    # Alias format_name with the relevant parameters pre-filled for convenience
    format(name) = format_name(
        name,
        test_config,
        observable_id,
        observable_name,
        is_logarithmic,
        axis,
        axis_range,
        fixed_var,
        fixed_value,
        secondary_var,
        secondary_value,
        fig_config["name_num_round_digits"],
    )

    # Initialize the figure with the appropriate defaults
    fig = plot(;
        xlabel = variable_names[axis] * "/t",  # Plots are all scaled by the hopping parameter
        ylabel = observable_name,
        title = format(fig_config["title"]),
        legend_position = num_series > 1 ? :best : :none,  # Only show legend if we have multiple series to distinguish between
        xscale = is_logarithmic ? :log10 : :identity,
        # We have to resample the colro palette to get the correct number of colors for the number of series we're plotting,
        # otherwise all the sampled colors are visually identical (because the palette is continuous)
        color_palette = resample(
            colorschemes[Symbol(fig_config["color_palette"])],
            num_series,
        ),
        # Apply user-specified parameters
        Utils.convert_strings_to_symbols(fig_config["plot_params"])...,
    )

    # Return the figure and the formatted filename
    return fig, format(fig_config["filename"])
end

"""
    load_data_matrix(file) -> DataMatrix

Loads the data from the specified CSV file and returns it as a DataMatrix struct.
The CSV file must be in the same format as created by DataHelpers
"""
function load_data_matrix(file::String)
    data_matrix = CSVUtil.load_csv_matrix(file)

    return DataMatrix(
        data_matrix[1, 2:end],
        # For some reason, this row only is parsed as strings rather than floats, so we have to manually convert it.
        parse.(Float64, data_matrix[2:end, 1]),
        data_matrix[2:end, 2:end],
    )
end

"""
    load_data_matrix(datadirs, Us, U, observable_id, prefix, order) -> DataMatrix

Computes the path to the relevant data file based on the provided parameters and loads it using load_data_matrix().
datadirs: The list of loaded data directories
Us: The corresponding U values for each data directory
U: The specific U value to load data for
observable_id: The id of the observable to load
prefix: The prefix of the resummation method to load
order: The order of the resummation method to load
"""
function load_data_matrix(
    datadirs::Vector{String},
    Us::Vector{Float64},
    U::Float64,
    observable_id::String,
    prefix::String,
    order::Int,
)
    # Find the index in datadirs that corresponds to the requested U value
    dir_idx = findfirst(==(U), Us)
    if dir_idx === nothing  # No loaded dataset corresponds to the requested U value!
        error("No data found for U = $U")
    end
    # Compute the filename, tag it onto the directory path, and load it
    return load_data_matrix(
        joinpath(datadirs[dir_idx], "$(observable_id)_$(prefix)_Order_$(order).csv"),
    )
end

"""
    load_fixed_range(range) -> Vector{Float64}

Parses one of the "fixed_" parameter specifications in the config file.
This can either be a literal range of values or a specification of the form ["range", start, stop, step].
"""
function load_fixed_range(range)
    if get(range, 1, nothing) == "range"
        return Float64(range[2]):Float64(range[3]):Float64(range[4])
    else
        return Float64.(range)
    end
end

"""
    is_log(range) -> Bool

Checks if the provided axis range specification indicates that the axis should be logarithmic.
"""
is_log(range) = length(range) == 3 && range[3] == "log"

"""
    get_data_for_params(datadirs, Us, observable_id, prefix, order, axis; T, u, U) -> (domain, range)

Loads the data for the specified parameters and returns the relevant domain and range of values for plotting.
The specific behavior of this function depends on the specified axis:
datadirs: The list of loaded data directories
Us: The corresponding U values for each data directory
observable_id: The id of the observable to load
prefix: The prefix of the resummation method to load
order: The order of the resummation method to load
axis: The variable to use as the x-axis. Must be one of :T, :u, or :U.

kwargs (all required):
    T: The value(s) of the temperature to load.
    u: The value(s) of the chemical potential to load.
    U: The value(s) of the interaction strength to load.
The type of the above values depend on the specified axis. The specified axis must be a range of the form [min value, max value],
indicating the bounds of the domain. The returned domain will lie entirely within these bounds, but might not span it.
The other values must be specific values, indicating which slice of the data to load.
"""

# See above docstring
function get_data_for_params(
    datadirs::Vector{String},
    Us::Vector{Float64},
    observable_id::String,
    prefix::String,
    order::Int,
    ::Val{:U};
    T::Float64,
    u::Float64,
    U::Vector{Float64},
)
    # Create output containers
    valid_Us = Int[]
    data_series = Float64[]

    # Search for valid interaction energies in the set that we have data for
    for U_val in Us
        # Clip to the specified domain
        if U_val < U[1] || U_val > U[2]
            continue
        end

        # Load the relevant data matrix
        data_mat = load_data_matrix(datadirs, Us, U_val, observable_id, prefix, order)

        # Pull out the requested T,u slice
        T_idx = findfirst(==(T), data_mat.T_vals)
        if T_idx === nothing
            @warn "No data found for U = $U_val and T = $T"
            continue
        end
        u_idx = findfirst(==(u), data_mat.u_vals)
        if u_idx === nothing
            @warn "No data found for U = $U_val and μ = $u"
            continue
        end
        data_series_val = data_mat.values[T_idx, u_idx]

        # Add it to the output
        push!(valid_Us, U_val)
        push!(data_series, data_series_val)
    end

    # Return the results
    return valid_Us, data_series
end

# See above docstring
function get_data_for_params(
    datadirs::Vector{String},
    Us::Vector{Float64},
    observable_id::String,
    prefix::String,
    order::Int,
    ::Val{:T};
    T::Vector{Float64},
    u::Float64,
    U::Float64,
)
    # Load the relevant data matrix
    data_mat = load_data_matrix(datadirs, Us, U, observable_id, prefix, order)

    # Pull out the requested u slice
    idx = findfirst(==(u), data_mat.u_vals)
    if idx === nothing
        error("No data found for μ = $u")
    end
    axis_range = data_mat.T_vals
    data_series = data_mat.values[:, idx]

    # Clip to the specified domain
    min_idx = findfirst(>=(T[1]), axis_range)
    max_idx = findlast(<=(T[2]), axis_range)

    # Return the results
    return axis_range[min_idx:max_idx], data_series[min_idx:max_idx]
end

# See above docstring
function get_data_for_params(
    datadirs::Vector{String},
    Us::Vector{Float64},
    observable_id::String,
    prefix::String,
    order::Int,
    ::Val{:u};
    T::Float64,
    u::Vector{Float64},
    U::Float64,
)
    # Load the relevant data matrix
    data_mat = load_data_matrix(datadirs, Us, U, observable_id, prefix, order)

    # Pull out the requested T slice
    idx = findfirst(==(T), data_mat.T_vals)
    if idx === nothing
        error("No data found for T = $T")
    end
    axis_range = data_mat.u_vals
    data_series = data_mat.values[idx, :]

    # Clip to the specified domain
    min_idx = findfirst(>=(u[1]), axis_range)
    max_idx = findlast(<=(u[2]), axis_range)

    # Return the results
    return axis_range[min_idx:max_idx], data_series[min_idx:max_idx]
end

"""
    create_plots_for_resummation_method!(figure, fig_config, ...)

Plots the data for a specific resummation method onto the provided figure.

figure: The Plots.jl figure object to plot onto
fig_config: The user-specified configuration for the figure
The rest of the parameters have the same meanings as in (@main)()
"""
function create_plots_for_resummation_method!(
    figure,
    fig_config::Dict{String,Any},
    datadirs::Vector{String},
    test_config::ED.TestConfiguration,
    Us::Vector{Float64},
    observable_id::String,
    prefix::String,
    orders::Vector{Int},
    axis::String,
    axis_range::Vector{Float64},
    fixed_var::String,
    fixed_value::Float64,
    secondary_var::String,
    secondary_value::Float64,
    is_overlay::Bool,
)
    # Preload and cache all the necessary datasets
    datasets = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()

    # Also load data for adjacent orders to enable convergence testing. (If necessary)
    orders_to_load = axis == "T" ? union(orders, orders .- 1, orders .+ 1) : orders

    for order in orders_to_load
        try
            datasets[order] = get_data_for_params(
                datadirs,
                Us,
                observable_id,
                prefix,
                order,
                Val(Symbol(axis));
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

    # Perform convergence tests if plotting against T.
    # This helps prevent the ranges from being dominated by non-converged data at low temperatures.
    cutoff_idxs = axis == "T" ?
        Dict(order => find_convergence_cutoff(datasets, order, fig_config["atol"], fig_config["rtol"]) for order in orders) :
        # Otherwise, just include the entire domain (no cutoff)
        Dict(order => 1 for order in orders)

    # Start creating labels for the plots
    # (For now, just the stuff that's common to all orders)

    # We only need to show the various bits of info if there's a possibility of ambiguity.
    show_prefix_label = length(fig_config["prefixes"]) > 1
    show_order_label = length(orders) > 1

    overlay_var_label = is_overlay ? "$secondary_var = $secondary_value" : ""
    prefix_label = show_prefix_label ? "$prefix" : ""

    # Plot the data for each order
    for order in orders
        # Load in the relevant data that we computed earlier
        x_vals, values = datasets[order]
        cutoff_idx = cutoff_idxs[order]
        # Also create a label for the order if necessary
        order_label = show_order_label ? "Order $order" : ""

        # Merge all the various label bits together
        label =
            [overlay_var_label, prefix_label, order_label] |>
            filter(!isempty) |>
            Fix2(join, ", ")

        # Plot it all onto the figure
        plot!(
            figure,
            x_vals[cutoff_idx:end] / test_config.t, # Normalize the x-axis by the hopping parameter
            values[cutoff_idx:end];
            label = label,
            # Apply user-specified parameters
            Utils.convert_strings_to_symbols(fig_config["series_params"])...,
        )
    end
end

"""
    find_convergence_cutoff(datasets, order, atol, rtol) -> Int

Finds the index at which the data for the specified order appears to converge (with respect to T).
This method first selects another order to compare against (see implementation for details) and then selects the lowest temperature
at which the two orders differ by no more than atol (absolutely) OR rtol (relatively) for all larger temperatures. If convergence fails,
this method prints a warning and returns 1, meaning that no cutoff is applied and all data points are plotted.
datasets: The preloaded datasets depending on T for all relevant orders.
order: The order to find the cutoff for
atol: The absolute tolerance for convergence
rtol: The relative tolerance for convergence
"""
function find_convergence_cutoff(datasets::Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}, order::Int, atol::Float64, rtol::Float64)
    # Find an order to compare against
    # Check, in order:
    # 1) The next order
    # 2) The next largest order that we have data for
    # 3) The previous order
    # 4) The next smallest order that we have data for
    # 5) Give up and don't use a cutoff
    cmp_order = order + 1
    if !haskey(datasets, cmp_order)
        larger_orders = filter(>(order), keys(datasets))
        if !isempty(larger_orders)
            cmp_order = minimum(larger_orders)
        end
    end
    if !haskey(datasets, cmp_order)
        cmp_order = order - 1
    end
    if !haskey(datasets, cmp_order)
        smaller_orders = filter(<(order), keys(datasets))
        if !isempty(smaller_orders)
            cmp_order = maximum(smaller_orders)
        end
    end
    if !haskey(datasets, cmp_order)
        @warn "Could not find another order for $observable_id with $prefix order $order to check convergence!"
        return 1
    end

    # Load the values for the order we're checking and the order we're comparing against
    values = datasets[order][2]
    cmp_values = datasets[cmp_order][2]

    # All orders (with the same params) should be generated from the same temperatures,
    # so the relevant indicies should be identical.
    cutoff = findlast(collect(zip(values, cmp_values))) do (value, cmp_value)
        # Find the last non-converged point to avoid cases where the data
        # looks converged but re-diverges
        # Don't bother with adding 1 to get to the first converged point because
        # the step size *should* be small.
        !(
            # Converged if:
            # 1) Difference is within tolerance, or
            abs(value - cmp_value) <= atol ||
            # 2) Difference is small relative to the current value
            abs((value - cmp_value) / value) <= rtol
        )
    end

    if cutoff === nothing  # No converged points found!
        @warn "$observable_id with $prefix order $order failed to converge!"
        return 1
    end

    return cutoff
end
