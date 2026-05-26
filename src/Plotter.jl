using HubbardDiagonalization: CSVUtil, DataHelpers, ExactDiagonalization as ED, Utils

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

function (@main)(args)
    cli_args = parse_cli()

    @info "Args: $cli_args"

    if cli_args["validate"] == "no"
        if length(cli_args["datadirs"]) > 1
            error("Cannot merge data from multiple directories without validation.")
        end
        @warn "Skipping validation of input data."

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

    config = TOML.parsefile(cli_args["plot_config"])

    validate_config(config, Us)

    output_dir = cli_args["outputdir"]
    Base.Filesystem.mkpath(output_dir)

    for (fig_name, fig_config) in config
        if fig_name == "defaults"
            continue
        end

        figure_requested = isempty(cli_args["figure"]) || fig_name in cli_args["figure"]
        if !figure_requested
            continue
        end

        @info "Generating plots for $fig_name..."
        observable_ids = fig_config["observables"]["ids"]
        observable_names = fig_config["observables"]["names"]

        num_series_per_param_set =
            sum(values(fig_config["prefixes"]) .|> Fix2(getindex, "orders") .|> length)

        for observable_id in observable_ids
            observable_name = observable_names[observable_id]

            @info "Plotting observable $observable_id ($observable_name)..."

            for axis in variables
                axis_range = get!(fig_config, "$(axis)_range", [])
                if isempty(axis_range)
                    continue
                end

                is_overlay, secondary_var, secondary_range =
                    if haskey(fig_config, "overlay")
                        var = fig_config["overlay"][1]

                        if var == axis
                            continue
                        end

                        range = load_range(fig_config["overlay"][2:end])
                        (true, var, range)
                    else
                        var = axis == "U" ? "u" : "U"
                        (false, var, load_range(fig_config["fixed_$(var)"]))
                    end

                fixed_var = setdiff(variables, [axis, var]) |> only
                fixed_value_range = load_range(fig_config["fixed_$(fixed_var)"])

                is_logarithmic = is_log(axis_range)
                axis_range = axis_range[1:2] .|> Float64

                for v1 in fixed_value_range
                    if is_overlay
                        figure, filename = init_figure(
                            fig_config,
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

                    for v2 in secondary_range
                        if !is_overlay
                            figure, filename = init_figure(
                                fig_config,
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

                        if !is_overlay
                            savefig(figure, joinpath(output_dir, filename))
                        end
                    end

                    if is_overlay
                        savefig(figure, joinpath(output_dir, filename))
                    end
                end
            end
        end
    end
end

function validate_config(config::Dict{String,Any}, Us::Vector{Float64})
    default_params = Dict(
        "observables" => Dict("ids" => [], "names" => Dict()),
        "fixed_U" => Us,
        "atol" => 1e-3,
        "rtol" => 0.05,
        "prefixes" => Dict(),
        "color_palette" => "coolwarm",
        "name_num_round_digits" => 3,
        "plot_params" => Dict(),
        "series_params" => Dict(),
    )

    defaults = get(config, "defaults", Dict())
    defaults = Utils.recursive_merge(default_params, defaults)
    config["defaults"] = defaults

    for figure in keys(config)
        if figure == "defaults"
            continue
        end
        config[figure] = Utils.recursive_merge(defaults, config[figure])
    end

    for (fig_name, fig_config) in config
        for observable_id in fig_config["observables"]["ids"]
            get!(fig_config["observables"]["names"], observable_id, observable_id)
        end

        for (prefix_name, prefix_data) in fig_config["prefixes"]
            if !haskey(prefix_data, "orders")
                error(
                    "Figure $fig_name is missing required 'orders' field in prefix $prefix_name.",
                )
            end
        end
    end
end

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
    axis_range = round.(axis_range, digits = name_num_round_digits)
    fixed_value = round(fixed_value, digits = name_num_round_digits)
    if secondary_value !== nothing
        secondary_value = round(secondary_value, digits = name_num_round_digits)
    end

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
        axis => axis_range,
        fixed_var => fixed_value,
        secondary_var => secondary_value,
    )
    return Utils.format_string(name, available_variables)
end

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

    fig = plot(;
        xlabel = variable_names[axis] * "/t",
        ylabel = observable_name,
        title = format(fig_config["title"]),
        legend_position = true,
        xscale = is_logarithmic ? :log10 : :identity,
        color_palette = resample(
            colorschemes[Symbol(fig_config["color_palette"])],
            num_series,
        ),
        Utils.convert_strings_to_symbols(fig_config["plot_params"])...,
    )
    return fig, format(fig_config["filename"])
end

function load_data_matrix(file::String)
    data_matrix = CSVUtil.load_csv_matrix(file)

    return DataMatrix(
        data_matrix[1, 2:end],
        # For some reason, this row only is parsed as strings rather than floats, so we have to manually convert it.
        parse.(Float64, data_matrix[2:end, 1]),
        data_matrix[2:end, 2:end],
    )
end

function load_data_matrix(
    datadirs::Vector{String},
    Us::Vector{Float64},
    U::Float64,
    observable_id::String,
    prefix::String,
    order::Int,
)
    dir_idx = findfirst(==(U), Us)
    if dir_idx === nothing
        error("No data found for U = $U")
    end
    return load_data_matrix(
        joinpath(datadirs[dir_idx], "$(observable_id)_$(prefix)_Order_$(order).csv"),
    )
end

function load_range(range)
    if get(range, 1, nothing) == "range"
        return Float64(range[2]):Float64(range[3]):Float64(range[4])
    else
        return Float64.(range)
    end
end

is_log(range) = length(range) == 3 && range[3] == "log"

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
    valid_Us = Int[]
    data_series = Float64[]
    for U_val in U
        data_mat = load_data_matrix(datadirs, Us, U_val, observable_id, prefix, order)
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
    ::Val{:T};
    T::Vector{Float64},
    u::Float64,
    U::Float64,
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
    ::Val{:u};
    T::Float64,
    u::Vector{Float64},
    U::Float64,
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
    atol = fig_config["atol"]
    rtol = fig_config["rtol"]

    datasets = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()

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

    cutoff_idxs = Dict{Int,Int}()
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
                larger_orders = filter(>(order), keys(datasets))
                if !isempty(larger_orders)
                    cmp_order = minimum(larger_orders)
                end
            end
            if !haskey(datasets, cmp_order)
                smaller_orders = filter(<(order), keys(datasets))
                if !isempty(smaller_orders)
                    cmp_order = maximum(smaller_orders)
                end
            end
            if !haskey(datasets, cmp_order)
                cutoff_idxs[order] = 1
                continue
            end

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
                    abs(value - cmp_value) / value <= rtol
                )
            end

            if cutoff === nothing
                @warn "$observable_id with $prefix order $order failed to converge!"
                cutoff_idxs[order] = 1
            else
                cutoff_idxs[order] = cutoff
            end
        end
    else
        cutoff_idxs = Dict(order => 1 for order in orders)
    end

    show_prefix_label = length(fig_config["prefixes"]) > 1
    show_order_label = length(orders) > 1

    overlay_var_label = is_overlay ? "$secondary_var = $secondary_value" : ""
    prefix_label = show_prefix_label ? "$prefix" : ""

    for order in orders
        x_vals, values = datasets[order]
        cutoff_idx = cutoff_idxs[order]
        order_label = show_order_label ? "Order $order" : ""

        label =
            [overlay_var_label, prefix_label, order_label] |>
            filter(!isempty) |>
            Fix2(join, ", ")
        label_kwarg = isempty(label) ? () : (label = label,)

        plot!(
            figure,
            x_vals[cutoff_idx:end] / test_config.t,
            values[cutoff_idx:end];
            label_kwarg...,
            Utils.convert_strings_to_symbols(fig_config["series_params"])...,
        )
    end
end
