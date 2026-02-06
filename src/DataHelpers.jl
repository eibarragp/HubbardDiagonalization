module DataHelpers

import ..CSVUtil
import ..ExactDiagonalization as ED

import CSV
import JSON3 as JSON

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
    export_observable_data(
        plot_config::Dict{String,Any},
        T_vals::AbstractVector{Float64},
        u_vals::AbstractVector{Float64},
        observable_data::Dict{String,Matrix{Float64}},
        config::ED.TestConfiguration,
        num_sites_str::String,
        using_nlce::Bool,
        output_dir::String,
    )

Exports the computed observable data to CSV files and generates plots based on the provided plot configuration.
"""
function export_observable_data(
    plot_config::Dict{String,Any},
    T_vals::AbstractVector{Float64},
    u_vals::AbstractVector{Float64},
    observable_data::Dict{String,Matrix{Float64}},
    config::ED.TestConfiguration,
    num_sites_str::String,
	using_nlce::Bool,
	output_dir::String,
)
    @info "Exporting observable data..."

    plot_width = plot_config["width"]
    plot_height = plot_config["height"]

    Base.Filesystem.mkpath(output_dir)
    for (observable_name, data_matrix) in observable_data
        labeled_matrix = hcat(["u/T", T_vals...], vcat(u_vals', data_matrix))
        CSV.write("$(output_dir)/$(observable_name).csv", CSV.Tables.table(labeled_matrix), writeheader = false)
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
                legend_position = length(observables) > 1 || (using_nlce ? true : length(csv_overlays) > 0),
                size = (plot_width, plot_height),
            )
            # Add each observable
            for observable_name in observables
				if using_nlce
					for order in plot_config["nlce_orders"]
						nlce_observable_name = "$(observable_name)_NLCE_Order_$order"
                        if length(observables) > 1
                            label = "$observable_name (Order $order)"
                        else
                            label = "Order $order"
                        end
						plot!(
							figure,
							x_axis,
							indexer(index, observable_data[nlce_observable_name]),
							labels = label,
						)
					end
				else
                    if !haskey(observable_data, observable_name)
                        @warn "Observable $observable_name not found; skipping."
                        continue
                    end
					plot!(
						figure,
						x_axis,
						indexer(index, observable_data[observable_name]),
						labels = observable_name,
					)
				end
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
                    plot!(
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
            plot_title = "t=$(config.t), $fixed_value_name=$fixed_value, U=$(config.U), num_sites=$num_sites_str, num_colors=$(config.num_colors)",
            layout = length(figures),
        )

        # Save the plot
        savefig(
            combined_figure,
            "$(output_dir)/observable_data_$(fixed_value_name)_$fixed_value.png",
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

function merge_results_with_nlce(cluster_file::String, data_dirs::Vector{String}, observables::Vector{String}, nlce_orders::Vector{Int})
	cluster_data = JSON.read(read(cluster_file, String))

	coefficients = []
	for cluster in sorted_clusters(cluster_data)
		push!(coefficients, cluster[3])
	end

	# This will create a matrix of the form [NLCE Order, Cluster Index]
	coefficients = hcat(coefficients...)

	if length(data_dirs) > size(coefficients, 2)
		error("More clusters provided ($(length(data_dirs))) than clusters in the cluster file ($(size(coefficients, 2)))!")
	end

	@info "Coefficients loaded for $(length(data_dirs)) clusters"
	@debug "Coefficients: $coefficients"

	# Load the data for each cluster and store it in a dictionary of the form Dict{ObservableName => [ClusterIndex => DataMatrix]}
	cluster_observable_data = Dict{String,Vector{Matrix{Float64}}}()
	for observable in observables
		cluster_data = Vector{Matrix{Float64}}()
		for data_dir in data_dirs
			data_path = "$(data_dir)/$(observable).csv"
			data_matrix = CSVUtil.load_csv_matrix(data_path)[2:end, 2:end]  # Skip the first row and column which contain the u and T values
            push!(cluster_data, data_matrix)
		end

		cluster_observable_data[observable] = cluster_data
	end

	# Compute the NLCE results!
	observable_data = Dict{String,Matrix{Float64}}()
	for order in nlce_orders
		order_coefficients = coefficients[order, :]
		@info "Computing NLCE results for order $order..."
		for observable in observables
			data_matrices = cluster_observable_data[observable]

			# Compute the NLCE result for this observable and order by summing over the contributions from each cluster weighted by their coefficients
			nlce_result = zeros(size(data_matrices[1]))
			for (cluster_idx, cluster_data) in enumerate(data_matrices)
				nlce_result .+= order_coefficients[cluster_idx] .* cluster_data
			end

			observable_data["$(observable)_NLCE_Order_$order"] = nlce_result
		end
	end

	return observable_data
end

function sorted_clusters(cluster_data::AbstractDict)
	clusters = []

	for (id, cluster) in cluster_data
		push!(clusters, (cluster..., parse(Int128, String(id))))
	end
	return sort(clusters; lt = (a, b) -> begin
			# Sort primarily by size
			if length(a[1]) != length(b[1])
				return length(a[1]) < length(b[1])
			end
			# Then sort by id
			return a[4] < b[4]
		end
	)
end

end
