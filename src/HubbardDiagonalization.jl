module HubbardDiagonalization

export main

export TestConfiguration
export default_observables, diagonalize_and_compute_observables, export_observable_data

# Include our submodules
include("CSVUtil.jl")
include("Graphs.jl")
include("StateEnumeration.jl")
include("SymmetricMatrices.jl")

using .CSVUtil
using .Graphs
using .StateEnumeration
using .SymmetricMatrices

# Import libraries
const use_unicode_plots = false

import CSV
import LinearAlgebra
import Logging
import Plots
import TOML
import ZipFile

using Base.Threads

# Set up plotting backend
using Plots
if use_unicode_plots
    import UnicodePlots
    unicodeplots()
end

"""
    A structure to hold the configuration parameters for the simulation.
"""
@kwdef struct TestConfiguration
    num_colors::Int
    t::Float64
    u_test::Float64
    U::Float64
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

    test_config = TestConfiguration(; convert_strings_to_symbols(params)...)
    t_vals = plot_config["T_min"]:plot_config["T_step"]:plot_config["T_max"]
    u_vals = plot_config["u_min"]:plot_config["u_step"]:plot_config["u_max"]
    graph = linear_chain(graph_config["num_sites"])

    observables, derived_observables, overlays = default_observables(test_config, graph)
    @info "Defined observables: $(union(keys(observables), keys(derived_observables), keys(overlays)))"

    observable_data = diagonalize_and_compute_observables(
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
    count_double_occupancies(state::Vector{Int}, num_colors::Int) -> Int

Counts the number of double occupancies in a given basis state.
"""
function count_double_occupancies(state::Vector{Int}, num_colors::Int) :: Int
    total = 0
    # Consider two colors at a time to interact
    for color_pair in enumerate_states(num_colors, 2)
        # Get indicies of the colors
        first_color = trailing_zeros(color_pair)
        second_color = trailing_zeros(color_pair ⊻ (1 << first_color))  # leading_zeros would depend on the bit-width, so we can't use it
        # Correct for 1-based indexing
        first_color += 1
        second_color += 1
        # The number of double occupancies is the number of bits set on both colors
        total += count_ones(state[first_color] & state[second_color])
    end
    return total
end

"""
    default_observables(test_config::TestConfiguration, graph::Graph)

Defines the "standard" set of observables to compute for the Hubbard model.
Returns three dictionaries:
- observables: Maps observable names to functions that compute them from a basis state.
- derived_observables: Maps observable names to functions that compute them from the dictionary of computed results for the normal observables.
- overlays: Maps observable names to functions that compute them from (B, u) pairs.
"""
function default_observables(test_config::TestConfiguration, graph::Graph)
    num_colors = test_config.num_colors
    U = test_config.U
    num_sites = Graphs.num_sites(graph)

    observables = Dict{String,Function}()

    observables["Num_Particles"] = state -> sum(count_ones(color) for color in state)
    observables["Filled States"] = state -> count_ones(reduce(&, state))
    observables["Double Occupancies"] = state -> count_double_occupancies(state, num_colors)

    observables["Energy"] = _ -> 0.0  # Will be handled specially

    observables["P_a"] =
        state ->
            count_ones(state[1]) *
            prod((1 - count_ones(state[c]) for c in 2:num_colors), init = 1)
    if num_colors >= 2
        observables["P_ab"] =
            state ->
                prod(count_ones(state[c]) for c in 1:2) *
                prod((1 - count_ones(state[c]) for c in 3:num_colors), init = 1)
    end
    if num_colors >= 3
        observables["P_abc"] =
            state ->
                prod(count_ones(state[c]) for c in 1:3) *
                prod((1 - count_ones(state[c]) for c in 4:num_colors), init = 1)
    end

    # Observables that can be calculated from other observables
    derived_observables = Dict{String,Function}()
    derived_observables["Local Moment"] =
        observable_data ->
            @. observable_data["Num_Particles"] - 2 * observable_data["Filled States"]
    derived_observables["Density"] =
        observable_data -> observable_data["Num_Particles"] ./ num_sites
    derived_observables["Entropy"] =
        observable_data -> zeros(Float64, size(observable_data["Num_Particles"])...)  # Will be handled specially

    # Additional Plots that can be directly calculated
    overlays = Dict{String,Function}()

    if num_sites == 1
        # Single-site Hubbard model exact solutions
        e0(n, u) = U * binomial(n, 2) - (u + (U / 2) * (num_colors - 1)) * n
        weighted_sum(B, u, f) =
            sum(binomial(num_colors, n) * f(n) * exp(-B * e0(n, u)) for n in 0:num_colors)
        z0(B, u) = weighted_sum(B, u, n -> 1)
        rho(B, u) = (1/z0(B, u)) * weighted_sum(B, u, n -> n)
        energy(B, u) =
            (1/z0(B, u)) *
            weighted_sum(B, u, n -> e0(n, u) + (u + (U/2) * (num_colors - 1)) * n)
        overlays["Actual Energy"] = energy
        overlays["Actual Entropy"] =
            (B, u) ->
                log(z0(B, u)) +
                B * (energy(B, u) - (u + (U/2) * (num_colors - 1)) * rho(B, u))
    end


    # Calculating the spin-spin correlation function / number correlation:
    #   C = sum_{σ≠τ} [ <n(i,σ)n(j,σ)> - <n(i,σ)n(j,τ)> ]
    # Defining observable:
    observables["C_spin"] = state -> begin
        num_edges = 0
        total = 0.0

        # Iterating over the edge pairs (nearest neighbor pairs) in the graph
        for (i, j) in Graphs.edges(graph)
            # Keep track of the edges (nearest neighbors) to average over later
            num_edges += 1
            same_color = 0
            diff_color = 0

            # if they are the same color
            # n(i,σ)n(j,σ)
            for σ in 1:num_colors
                # finding the occupancy at the state (represented by a bitmask)
                niσ = (state[σ] >> (i - 1)) & 1
                njσ = (state[σ] >> (j - 1)) & 1

                #contributes (num_colors - 1) if there is a spin up up or a spin down down at i and j
                same_color += niσ * njσ * (num_colors - 1)
            end

            # if different colors
            # n(i,σ)n(j,τ), σ ≠ τ
            for σ in 1:num_colors
                # finding the occupancy at the state
                niσ = (state[σ] >> (i - 1)) & 1

                # looping over the possible other colors
                for τ in 1:num_colors
                    if σ == τ
                        ;
                        continue;
                    end
                    njτ = (state[τ] >> (j - 1)) & 1
                    # contributes non-zero only if the two sites have different spins
                    diff_color += niσ * njτ
                end
            end

            # finding the total contribution
            # positive if both neighbors have the same spins (up up, down down)
            # negative if both neighbors have different spins (up down, down up)
            total += (same_color - diff_color)
        end
        # return average
        return num_edges == 0 ? 0.0 : total / (num_edges)
    end

    return observables, derived_observables, overlays
end

"""
    diagonalize_and_compute_observables(
        T_vals::AbstractVector{Float64},
        u_vals::AbstractVector{Float64},
        config::TestConfiguration,
        graph::Graph,
        observables::Dict{String,Function},
        derived_observables::Dict{String,Function},
        overlays::Dict{String,Function},
    ) -> Dict{String,Matrix{Float64}}

Performs the exact diagonalization of the Hubbard Hamiltonian and computes the specified observables (see default_observables() for specification) over the given ranges of temperature and chemical potential.
"""
function diagonalize_and_compute_observables(
    T_vals::AbstractVector{Float64},
    u_vals::AbstractVector{Float64},
    config::TestConfiguration,
    graph::Graph,
    observables::Dict{String,Function},
    derived_observables::Dict{String,Function},
    overlays::Dict{String,Function},
)
    # Load parameters into the local scope
    num_colors = config.num_colors
    t = config.t
    u_test = config.u_test
    U = config.U

    # Compute some useful quantities
    num_temps = length(T_vals)
    num_us = length(u_vals)
    num_sites = Graphs.num_sites(graph)
    B = 1 ./ T_vals
    N_max_fermions = num_colors * num_sites

    @info "Initialized with t=$t, U=$(U), T_vals=$T_vals, u_vals=$u_vals"
    @debug begin
        "  Graph edges: $(Graphs.edges(graph))"
    end
    @debug "N_max_fermions=$N_max_fermions"

    """
    	create_observable_data_map(include_derived::Bool, include_overlays::Bool, size::Int...)

    include_derived: Whether to include derived observables in the map.
    include_overlays: Whether to include overlay functions in the map.
    size: The size of the vectors to create for each observable.

    A convenience function to create an empty map from observable names to data vectors.
    """
    # We're going to need a few of these. Might as well make it a function.
    function create_observable_data_map(
        include_derived::Bool,
        include_overlays::Bool,
        size::Int...,
    )
        storage_type = typeof(zeros(Float64, size))
        map = Dict{String,storage_type}()
        for observable_name in keys(observables)
            map[observable_name] = zeros(Float64, size)
        end
        if include_derived
            for derived_name in keys(derived_observables)
                map[derived_name] = zeros(Float64, size)
            end
        end
        if include_overlays
            for overlay_name in keys(overlays)
                map[overlay_name] = zeros(Float64, size)
            end
        end
        return map
    end

    # Precalculate configurations and block sizes so we can preallocate memory and do multithreading
    system_configurations = [
        (N_fermions, color_configuration) for N_fermions in 0:N_max_fermions for
        color_configuration in color_configurations(N_fermions, num_sites, num_colors)
    ]
    block_sizes = [
        prod(binomial(num_sites, n) for n in color_configuration) for
        (_, color_configuration) in system_configurations
    ]
    # Offset into the state arrays for each block
    size_offset = cumsum(block_sizes) .- block_sizes

    # Initialize some containers to store the generated data
    num_computed_states = sum(block_sizes)
    weights = zeros(Float64, num_temps, num_computed_states)  # Weights for each state
    n_fermion_data = zeros(Int, num_computed_states)  # Number of fermions for each state (used for re-weighting)
    observable_data = create_observable_data_map(false, false, num_computed_states)

    @info "Computing Hamiltonian blocks and observables..."
    # The number of fermions and the color configuration are conserved over tunneling,
    # so we can break the Hamiltonian into blocks labeled by these two quantities
    @threads :greedy for (config_idx, (N_fermions, color_configuration)) in
                         collect(enumerate(system_configurations))
        # Size of the Hamiltonian block
        L = block_sizes[config_idx]
        H = SymmetricMatrix(L)  # Use custom "SymmetricMatrix" type to save memory at the cost of speed
        observables_basis = create_observable_data_map(false, false, L)  # Compute the observables for each state as we build the matrix

        # Compute Hamiltonian matrix elements between all pairs of states
        # enumerate_multistate returns elements in a consistent order, so
        # as long as we're consistent, the matrix elements will be in the right place
        # state_i and state_j are arrays of integers, where each integer is a bitmask
        # representing the occupation of each site for a given color
        for (i, state_i) in enumerate(enumerate_multistate(num_sites, color_configuration))
            # Note: We're going to cut this inner loop off early since the matrix is symmetric
            for (j, state_j) in
                enumerate(enumerate_multistate(num_sites, color_configuration))
                @debug begin
                    "Computing H[$i,$j] between states:\n  state_i=$(digits.(state_i, base=2, pad=num_sites))\n  state_j=$(digits.(state_j, base=2, pad=num_sites))"
                end
                if i == j  # Because enumerate_multistate is consistent, if the indices are equal, the states are equal
                    # Diagonal element
                    H[i, i] =
                        -u_test * N_fermions +
                        U * count_double_occupancies(state_i, num_colors)

                    break  # No need to compute upper-triangular elements
                else
                    # Off-diagonal element
                    H[j, i] = 0.0   # Use j,i to be efficient with column-major storage

                    # Hopping term
                    # First, compute the difference between the two states
                    # This is a bitmask where bits are 1 if a fermion appeared or disappeared
                    # Because the total number of fermions is fixed, if only two
                    # bits are set, then one fermion hopped from one site to another
                    diff = state_i .⊻ state_j
                    # Figure out which color hopped
                    # If more than one color hopped, then each individual hopping term is zero
                    # so their sum is also zero
                    hopped_color = 0
                    for color in 1:num_colors
                        if diff[color] != 0
                            if count_ones(diff[color]) != 2 || hopped_color != 0
                                # More than one color or more than one fermion hopped,
                                # so this matrix element is zero
                                hopped_color = -1
                                break
                            end
                            hopped_color = color
                        end
                    end
                    @assert hopped_color != 0  # States must be different
                    if hopped_color == -1
                        continue
                    end

                    # Get the sites involved in the hop
                    site_1 = trailing_zeros(diff[hopped_color])
                    site_2 = trailing_zeros(diff[hopped_color] ⊻ (1 << site_1))  # leading_zeros would depend on the bit-width, so we can't use it
                    # Correct for 1-based indexing
                    site_1 += 1
                    site_2 += 1
                    @assert site_1 != site_2

                    @debug begin
                        "Considering hop of color $hopped_color from site $site_1 to site $site_2"
                    end

                    # Verify that the hop is allowed by the graph
                    if has_edge(graph, site_1, site_2)
                        # If so, the sign will be flipped if an odd number of
                        # spots are occupied *between* the two sites.
                        # Note that this refers to the representation of the
                        # spots not how they're related by the graph

                        occupied_sites = state_i[hopped_color] & state_j[hopped_color]
                        # Create a mask with 1s in all bits between site_1 and site_2.
                        # Note that julia's one-based indexing actually works out here
                        # If site=2 (so it's referring to the 2nd bit in the bitmask),
                        # Then, (1 << site) = 0b100, and (1 << site) - 1 = 0b011
                        # (Keep in mind that digits is interpreting this in little-endian,
                        # so the second-least-significant bit is the "2nd" bit)
                        # With this in mind, we can calculate the mask by taking the
                        # mask for all bits at or below site_2 (which is larger than site_1)
                        # and ANDing it with all of the bits above site_1 (which is just the
                        # same algorithm negated)
                        # Note that because the hop occurs between site_1 and site_2,
                        # both site_1 and site_2 are 0 in occupied_sites, so it doesn't
                        # matter if we include them in the mask or not.
                        bitween_mask = ((1 << site_2) - 1) & ~((1 << site_1) - 1)
                        sign = iseven(count_ones(occupied_sites & bitween_mask)) ? 1 : -1
                        @debug begin
                            "  Hop is allowed by graph! sign=$sign"
                        end
                        H[j, i] = sign * (-t)
                    end
                end
            end

            # Now that we've constructed the row for state_i, compute the observables
            for (observable_name, observable_function) in observables
                # Pre-compute the observable for this basis state
                observables_basis[observable_name][i] = observable_function(state_i)
            end
        end

        num_permutations = num_configuration_permutations(color_configuration)

        @debug begin
            "N_fermions=$N_fermions, color_configuration=$(color_configuration), L=$L, num_configuration_permutations=$(num_permutations), H=$H"
        end

        # Diagonalize the Hamiltonian block
        H_symmetric_view = LinearAlgebra.Symmetric(H, :U)
        # Annoyingly, eigen() forces us to store all the eigenvectors
        # in memory at once, but I can't find a good way around this
        # Even the builtin `eigvals`/`eigvecs` functions are just
        # wrappers around this.
        eigen_data = LinearAlgebra.eigen(H_symmetric_view)

        @debug begin
            msg = "observable_basis_data:\n"
            for (name, data) in observables_basis
                msg *= "  $name: $data\n"
            end
            msg
        end

        # Compute and store observables for each eigen-state
        offset = size_offset[config_idx]
        for (i, (eigen_val, eigen_vec)) in
            enumerate(zip(eigen_data.values, eachcol(eigen_data.vectors)))
            @debug begin
                "  eigen_val=$eigen_val, eigen_vec=$eigen_vec"
            end

            idx = offset + i

            # eigen() returns normalized eigenvectors, so we don't need to do any normalization here

            # Weight data of this state in the partition function
            weight = num_permutations * exp.(-B * eigen_val)
            weights[:, idx] = weight
            n_fermion_data[idx] = N_fermions

            # Compute each observable for this state
            for (observable_name, observable_basis_data) in observables_basis
                if observable_name == "Energy"
                    # The energy is just the eigenvalue
                    observable_data[observable_name][idx] = eigen_val
                else
                    # Because we already computed the observables for each basis state,
                    # we can just do a weighted sum over those based on the eigenvector components
                    observable_value = sum(@. observable_basis_data * eigen_vec * eigen_vec)
                    observable_data[observable_name][idx] = observable_value
                end
            end
        end
    end

    @info "Computed data for $num_computed_states states."

    @info "Computing derived observables..."

    for (observable_name, observable_function) in derived_observables
        observable_data[observable_name] = observable_function(observable_data)
    end

    @debug begin
        msg = "Computed data:\n"
        msg *= "  weights: $weights\n"
        msg *= "  n_fermion_data: $n_fermion_data\n"
        for (name, data) in observable_data
            msg *= "  $name: $data\n"
        end
        msg
    end

    @info "Computing observables over range of u..."

    # For computing the observables, its more efficient to have the temperatures in the rows
    weights = weights'

    u_shift = (U/2) * (num_colors - 1)  # Shift observables so that density=N/2 at u=0
    # Create a new container to store the observable values at each u
    computed_observable_values = create_observable_data_map(true, true, num_temps, num_us)
    @threads for (i, u) in collect(enumerate(u_vals))
        # The value that has to be added to u_test to shift to the desired u
        u_datapoint_shift = u - u_test + u_shift

        # Re-weight the data according to the new u value
        # Because of the way broadcasting works, we have to construct the exponents first
        # (or Julia gets confused). The below syntax makes a matrix of all combinations of
        # elements from B and n_fermion_data
        correction_exponents = permutedims(-B * -u_datapoint_shift) .* n_fermion_data
        # This makes corrected_weights a matrix with indexing [state, temp]
        corrected_weights = @. exp(correction_exponents) * weights

        # Sum over all the values in each column to get a vector of partition functions for each temperature
        Z = sum(corrected_weights, dims = 1)

        @debug begin
            "u=$u, corrected_weights=$corrected_weights, Z=$Z"
        end

        # Compute each observable
        for (observable_name, observable_values) in observable_data
            if observable_name == "Energy"
                # Now, update the energy. The free energy is H + u * N, so it works out to
                observable_values = (-u_test * n_fermion_data) .+ observable_values
                @debug begin
                    "  Updated Energy data: $observable_values"
                end

                # While we're here, calculate the entropy
                # The expectation value of the Hamiltonian depends on u, so we have to shift it here
                internal_energy_values =
                    (-u_datapoint_shift * n_fermion_data) .+ observable_values
                internal_energy_expectations =
                    sum(corrected_weights .* internal_energy_values; dims = 1)
                normalized_internal_energy_expectations = internal_energy_expectations ./ Z
                entropy_expectation =
                    @. normalized_internal_energy_expectations * B' + log(Z)
                computed_observable_values["Entropy"][:, i] = entropy_expectation
                @debug begin
                    "  Entropy: $entropy_expectation (Internal Energy Values: $internal_energy_values Internal Energy Expectation: $normalized_internal_energy_expectations)"
                end

                # Fallthrough/continue to store the energy value based on the corrected data
            elseif observable_name == "Entropy"
                # Calculated above
                continue
            end

            # Compute the expectation value of each observable
            expectation_values = sum(corrected_weights .* observable_values; dims = 1)
            normalized_expectation_values = expectation_values ./ Z

            @debug begin
                "  $observable_name: $normalized_expectation_values"
            end

            computed_observable_values[observable_name][:, i] .=
                normalized_expectation_values'
        end

        for (overlay_name, overlay_function) in overlays
            computed_observable_values[overlay_name][:, i] = overlay_function.(B, u)
        end
    end

    return computed_observable_values
end

"""
    export_observable_data(
        plot_config::Dict{String,Any},
        T_vals::AbstractVector{Float64},
        u_vals::AbstractVector{Float64},
        observable_data::Dict{String,Matrix{Float64}},
        config::TestConfiguration,
        graph::Graph,
    )

Exports the computed observable data to CSV files and generates plots based on the provided plot configuration.
"""
function export_observable_data(
    plot_config::Dict{String,Any},
    T_vals::AbstractVector{Float64},
    u_vals::AbstractVector{Float64},
    observable_data::Dict{String,Matrix{Float64}},
    config::TestConfiguration,
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

end
