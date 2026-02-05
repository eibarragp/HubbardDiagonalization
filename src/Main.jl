# Include Libraries
using HubbardDiagonalization:
	CSVUtil,
	Graphs, Graphs.Graph,
	ExactDiagonalization as ED

import CSV
import Logging
import TOML

using Base.Threads

# Set up plotting backend
using Plots
const use_unicode_plots = false
if use_unicode_plots
    import UnicodePlots
    unicodeplots()
end

"""
    convert_strings_to_symbols(dict::Dict{String,Any}) -> Dict{Symbol,Any}

Convert the keys of a dictionary from `String` to `Symbol`. This allows a dictionary loaded from a TOML file to be used as a keyword argument list.
"""
function convert_strings_to_symbols(dict::Dict{String,Any})
    new_dict = Dict{Symbol,Any}()
    for (key, value) in dict
        new_dict[Symbol(key)] = value
    end
    return new_dict
end

"""
    main(args::Vector{String}) -> Int

The main entry point! Handles the high-level control flow of the program.
"""
function (@main)(args)
    # Install our own logger for the duration of the program
    old_logger = Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Debug))
    if "--debug" in args
        @warn "Running in debug mode!"
        sleep(5)  # Give user time to see the warning
    else
        # In normal mode disable debug logging
        Logging.disable_logging(Logging.Debug)
    end

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

    test_config = ED.TestConfiguration(; convert_strings_to_symbols(params)...)
    t_vals = plot_config["T_min"]:plot_config["T_step"]:plot_config["T_max"]
    u_vals = plot_config["u_min"]:plot_config["u_step"]:plot_config["u_max"]
    graph = Graphs.linear_chain(graph_config["num_sites"])

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

    # Normalize energy and entropy per site
    observable_data["Energy"] ./= Graphs.num_sites(graph)
    observable_data["Entropy"] ./= Graphs.num_sites(graph)

    export_observable_data(plot_config, t_vals, u_vals, observable_data, test_config, graph)

    @info "Done."

    # Restore the old logger after we're done
    Logging.global_logger(old_logger)

    return 0
end

"""
    export_observable_data(
        plot_config::Dict{String,Any},
        T_vals::AbstractVector{Float64},
        u_vals::AbstractVector{Float64},
        observable_data::Dict{String,Matrix{Float64}},
        config::ED.TestConfiguration,
        graph::Graph,
    )

Exports the computed observable data to CSV files and generates plots based on the provided plot configuration.
"""
function export_observable_data(
    plot_config::Dict{String,Any},
    T_vals::AbstractVector{Float64},
    u_vals::AbstractVector{Float64},
    observable_data::Dict{String,Matrix{Float64}},
    config::ED.TestConfiguration,
    graph::Graph,
)
    @info "Exporting observable data..."

    num_sites = Graphs.num_sites(graph)
    plot_width = plot_config["width"]
    plot_height = plot_config["height"]

    Base.Filesystem.mkpath("output")
    for (observable_name, data_matrix) in observable_data
        labeled_matrix = hcat(["u/T", T_vals...], vcat(u_vals', data_matrix))
        CSV.write("output/$(observable_name).csv", CSV.Tables.table(labeled_matrix))
    end

    # Load csv (if requested)
    if haskey(plot_config, "overlay_data") &&
       lowercase(plot_config["overlay_data"]) != "none"
        csv_path = plot_config["overlay_data"]
        @info "Loading CSV overlay data from $csv_path..."
        csv_overlay_u_vals, csv_overlay_T_vals, csv_overlay_data =
            CSVUtil.load_overlay_data(csv_path)

        # Clip csv data to the range of our computed data
        u_indices =
            CSVUtil.get_clip_indices(csv_overlay_u_vals, minimum(u_vals), maximum(u_vals))
        T_indices =
            CSVUtil.get_clip_indices(csv_overlay_T_vals, minimum(T_vals), maximum(T_vals))

        csv_overlay_u_vals = csv_overlay_u_vals[u_indices]
        csv_overlay_T_vals = csv_overlay_T_vals[T_indices]
        for (key, value) in csv_overlay_data
            csv_overlay_data[key] = value[T_indices, u_indices]
        end
    end

    @info "Plotting observables..."

    """
        plot_data(fixed_value::Float64, fixed_value_type::Symbol)

    Creates plots for the given fixed value of either T or u.
    fixed_value: The value of either T or u that's fixed.
    fixed_value_type: Either :T or :u, indicating which variable is fixed.
    """
    function plot_data(fixed_value::Float64, fixed_value_type::Symbol)
        # Because we know which value is fixed, we can set up the axes accordingly
        if fixed_value_type == :T
            fixed_value_name = "T"  # Name of the fixed variable
            fixed_axis = T_vals     # Value array for the fixed variable (to get the index in the computed data)
            x_axis_name = "u"       # x-axis label
            x_axis = u_vals         # x-axis values for the generated data
            indexer = (index, data) -> data[index, :]  # Function to index into the data matrices
            if @isdefined(csv_overlay_data)
                fixed_overlay_axis = csv_overlay_T_vals  # Value array for the fixed variable (to get the index in the CSV data)
                x_overlay_axis = csv_overlay_u_vals  # x-axis values for the CSV data
            end
        elseif fixed_value_type == :u
            # Same as above, but altered for fixed u
            fixed_value_name = "u"
            fixed_axis = u_vals
            x_axis = T_vals
            x_axis_name = "T"
            indexer = (index, data) -> data[:, index]
            if @isdefined(csv_overlay_data)
                fixed_overlay_axis = csv_overlay_u_vals
                x_overlay_axis = csv_overlay_T_vals
            end
        else
            error("Invalid fixed_value_type: $fixed_value_type")
        end

        # Find the indices corresponding to the fixed value
        index = CSVUtil.find_index(fixed_axis, fixed_value, fixed_value_name)
        if index == -1
            @warn "Fixed value $fixed_value_name=$fixed_value not found in computed data; skipping plot."
            return
        end
        csv_index =
            @isdefined(csv_overlay_data) ?
            CSVUtil.find_index(fixed_overlay_axis, fixed_value, fixed_value_name) :
            # If no CSV data was loaded, just say we didn't find the axis
            -1

        # Array to store the subplots
        figures = []
        # For every subplot configuration
        for (observables, csv_overlays) in plot_config["observables"]
            # Create a plot
            figure = plot(
                xlabel = x_axis_name,
                ylabel = "Observable Value",
                title = join(observables, ", "),
                # Only show legend if it would be confusing without one
                legend_position = length(observables) > 1 || length(csv_overlays) > 1,
                size = (plot_width, plot_height),
            )
            # Add each observable
            for observable_name in observables
                if !haskey(observable_data, observable_name)
                    @warn "Observable $observable_name not found; skipping."
                    continue
                end
                figure = plot!(
                    figure,
                    x_axis,
                    indexer(index, observable_data[observable_name]),
                    labels = observable_name,
                )
            end
            # If valid data exists in the CSV, add those observables too
            if csv_index != -1
                for csv_overlay_name in csv_overlays
                    if !haskey(csv_overlay_data, csv_overlay_name)
                        @warn "Overlay $csv_overlay_name not found; skipping."
                        continue
                    end
                    data = indexer(csv_index, csv_overlay_data[csv_overlay_name])
                    filtered_indices = 1:length(x_overlay_axis)
                    figure = plot!(
                        figure,
                        x_overlay_axis[filtered_indices],
                        data[filtered_indices],
                        labels = "$(csv_overlay_name) (CSV Overlay)",
                        linestyle = :dash,
                    )
                end
            end

            # Add the subplot to the list
            push!(figures, figure)
        end

        # Merge all figures into a single plot
        combined_figure = plot(
            figures...;
            plot_title = "t=$(config.t), $fixed_value_name=$fixed_value, U=$(config.U), num_sites=$num_sites, num_colors=$(config.num_colors)",
            layout = length(figures),
        )

        # Save the plot
        savefig(
            combined_figure,
            "output/observable_data_$(fixed_value_name)_$fixed_value.png",
        )
    end

    # Create plots for each T value
    for T in plot_config["T_fixed_plots"]
        plot_data(T, :T)
    end

    # And create plots for each u value
    for u in plot_config["u_fixed_plots"]
        plot_data(u, :u)
    end

    return
end
