module ExactDiagonalization

export TestConfiguration
export apply_normalizations!, default_observables, diagonalize_and_compute_observables

using ..Graphs
using ..StateEnumeration
using ..SymmetricMatrices
import ..Utils

import CSV
import Dates
import LinearAlgebra

using Base.Threads

"""
    Defines when and how an observable should be computed and plotted.
    - Direct: An observable that is computed purely from the basis states.
        - func should be a function that takes a basis state. args is unused.
    - Derived: An observable that is computed from the values of the direct observables. Must not be energy-dependent.
        - func should accept the observable values listed in args (in order) as arbitrary-length arrays
            and perform element-wise operations to compute the derived observable.
    - Energy_Dependent: An observable that is computed directly from the basis states and eigenvalues.
        - func should be a function that takes a value of mu; and an extrinsic energy, intrinsic energy, and the number of particles,
            all as arrays containing the value at each temperature. It should then perform element-wise operations to compute the derived observable.
            args is unused.
    - Expectation_Derived: An observable that is computed from the expectation values of other observables
        - func should accept (in-order) the value of u, a vector of beta values, the value of the partition function,
            and the observable values listed in args (in order) as arrays containing the expectation value at each temperature.
    - Overlay: An observable that is computed directly from (B, u) pairs and is meant to be plotted as an overlay on top of the other observables.
        - func should accept (in-order) a value of B and a value of u and return the value of the observable at that point. args is unused.
"""
@enum ObservableType begin
    ObservableType_Direct
    ObservableType_Derived
    ObservableType_Energy_Dependent
    ObservableType_Expectation_Derived
    ObservableType_Overlay
end

"""
    A structure describing an observable to be computed and plotted.
"""
struct Observable
    """
        The type of the observable. (See ObservableType for options)
    """
    type::ObservableType
    """
        The other observables that this observable depends on. (Used for Derived and Expectation_Derived types)
    """
    args::AbstractVector{String}
    """
        The function that computes the observable. See ObservableType for specifications.
    """
    func::Function
end
# Identity operator to make some typing magic work later
Observable(o::Observable) = o

"""
    A structure to hold the configuration parameters for the simulation.
"""
@kwdef struct TestConfiguration
    num_colors::Int
    t::Float64
    u_test::Float64
    U::Float64
end
Base.:(==)(a::TestConfiguration, b::TestConfiguration) =
    a.num_colors == b.num_colors && a.t == b.t && a.u_test == b.u_test && a.U == b.U

"""
    count_double_occupancies(state::Vector{Int}, num_colors::Int) -> Int

Counts the number of double occupancies in a given basis state.
"""
function count_double_occupancies(state::Vector{Int}, num_colors::Int)::Int
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

    observable_tmp_type = Function
    observable_with_args_tmp_type = Tuple{AbstractVector{String},Function}

    direct_observables = Dict{String,observable_tmp_type}(
        "Density" => state -> sum(count_ones(color) for color in state),
        "n^2" => state -> sum(count_ones(color) for color in state)^2,
        "Filled States" => state -> count_ones(reduce(&, state)),
        "Double Occupancies" => state -> count_double_occupancies(state, num_colors),
        "na_nb" => state -> count_double_occupancies(state, num_colors),
        # Correlation Functions
        "niσ_njσ" =>
            state -> begin
                total = 0.0

                # For every edge
                for (i, j) in Graphs.edges(graph)
                    # Make a mask with the bits for sites i and j set so we can easily extract them
                    mask = (1 << (i - 1)) | (1 << (j - 1))

                    # For every color, check if both sites are occupied
                    total += count(color -> color & mask == mask, state)
                end

                return total
            end,
        "niσ_njτ" => state -> begin
            total = 0.0

            # For every edge
            for (i, j) in Graphs.edges(graph)
                # Create bitmasks for site i and site j
                imask = 1 << (i - 1)
                jmask = 1 << (j - 1)

                # For every pair of colors, check if site i is occupied in color σ and site j is occupied in color τ
                for σ in 1:num_colors
                    niσ = state[σ] & imask == imask

                    for τ in 1:num_colors
                        if σ != τ
                            njτ = state[τ] & jmask == jmask

                            total += niσ && njτ  # Julia will implicitly convert Bool to Int here
                        end
                    end
                end
            end

            return total
        end,
    )

    # Create a function to generate P_a (1-P_b) + P_b (1-P_a) type observables for any combination of colors
    p_function(n_check_colors) =
        state -> begin
            # Create array of P_n
            color_count_state = count_ones.(state)
            # Create array of (1-P_n)
            color_count_state_inverse = 1 .- color_count_state
            total = 0.0
            # For every combination of n_check_colors colors,
            for color_combination in enumerate_states(num_colors, n_check_colors)
                active_colors = digits(color_combination, base = 2, pad = num_colors)
                # Take the product of P_n for the active colors and (1-P_n) for the inactive colors
                total +=
                    prod(active_colors .* color_count_state) *
                    prod((1 .- active_colors) .* color_count_state_inverse)
            end
            # Normalize by num permutations
            return total / binomial(num_colors, n_check_colors)
        end

    direct_observables["P_a"] = p_function(1)
    if num_colors >= 2
        direct_observables["P_ab"] = p_function(2)
    end
    if num_colors >= 3
        direct_observables["P_abc"] = p_function(3)
    end

    # Observables that can be calculated from other observables
    derived_observables = Dict{String,observable_with_args_tmp_type}(
        "C_spin" => (
            ["niσ_njσ", "niσ_njτ"],
            (same_neighbors, different_neighbors) ->
                (num_colors - 1) * same_neighbors - different_neighbors,
        ),
    )

    # Observables that depend on the energy
    energy_dependent_observables = Dict{String,observable_tmp_type}(
        "Energy" => (_, energy, _, _) -> energy,
        "H" => (_, _, H, _) -> H,
        "H^2" => (_, _, H, _) -> H .^ 2,
        "Hn" => (_, _, H, n_particles) -> H .* n_particles,
    )

    # Observables that depend on the expectation values of other observables
    expectation_derived_observables = Dict{String,observable_with_args_tmp_type}(
        "Entropy" => (
            ["H"],
            (_, B, Z, internal_energy_expectation) ->
                @. log(Z)' + B * internal_energy_expectation
        ),
        "ΔH^2" => (["H", "H^2"], (_, _, _, H, H2) -> @. H2 - (H^2)),
        "Δn^2" => (["Density", "n^2"], (_, _, _, n, n2) -> @. n2 - (n^2)),
        "Specific Heat" => (
            ["Density", "n^2", "H", "H^2", "Hn"],
            (_, B, _, n, n2, H, H2, Hn) ->
                @. B^2 * ((H2 - H^2) - (Hn - n * H)^2 / (n2 - n^2))
        ),
        "Compressibility" =>
            (["Density", "n^2"], (_, B, _, n, n2) -> @. B * (n2 - n^2)),
    )

    # Additional Plots that can be directly calculated
    overlays = Dict{String,observable_tmp_type}()

    if num_sites == 1
        # Single-site Hubbard model exact solutions
        e0(n, u) = U * binomial(n, 2) - u * n
        weighted_sum(B, u, f) =
            sum(binomial(num_colors, n) * f(n) * exp(-B * e0(n, u)) for n in 0:num_colors)
        z0(B, u) = weighted_sum(B, u, n -> 1)
        rho(B, u) = (1/z0(B, u)) * weighted_sum(B, u, n -> n)
        energy(B, u) = (1/z0(B, u)) * weighted_sum(B, u, n -> e0(n, u) + u * n)
        overlays["Actual Energy"] = energy
        overlays["Actual Entropy"] =
            (B, u) -> log(z0(B, u)) + B * (energy(B, u) - u * rho(B, u))
    end

    # Define some temporary functions to convert the above functions into Observables in what I think is a more-readable way.
    to_observable(type, with_args) =
        with_args ? ((args, func),) -> Observable(type, args, func) :
        func -> Observable(type, String[], func)
    map(f, d) = Utils.map_dict_values(Observable, f, d)
    return merge(
        map(to_observable(ObservableType_Direct, false), direct_observables),
        map(to_observable(ObservableType_Derived, true), derived_observables),
        map(
            to_observable(ObservableType_Energy_Dependent, false),
            energy_dependent_observables,
        ),
        map(
            to_observable(ObservableType_Expectation_Derived, true),
            expectation_derived_observables,
        ),
        map(to_observable(ObservableType_Overlay, false), overlays),
    )
end

"""
    apply_normalizations!(observable_data::Dict{String,Matrix{Float64}}, graph::Graph)

Applies standard normalization factors to computed observables. This step should be omitted when
using NLCE, since NLCE will take care of the normalization for us.
"""
function apply_normalizations!(observable_data::Dict{String,Matrix{Float64}}, graph::Graph)
    num_sites = Graphs.num_sites(graph)
    num_edges = length(Graphs.edges(graph))

    if haskey(observable_data, "Density")
        observable_data["Density"] ./= num_sites
    end
    if haskey(observable_data, "Energy")
        observable_data["Energy"] ./= num_sites
    end
    if haskey(observable_data, "Entropy")
        observable_data["Entropy"] ./= num_sites
    end
    if haskey(observable_data, "C_spin") && num_edges > 0
        observable_data["C_spin"] ./= num_edges
    end

    return observable_data
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
    observables::Dict{String,Observable},
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
    	create_observable_data_map(include_list::Union{Nothing, AbstractVector{ObservableType}}, size::Int...)

    include_list: list of observable types to include in the map or Nothing to include all types.
    size: The size of the vectors to create for each observable.

    A convenience function to create an empty map from observable names to data vectors.
    """
    # We're going to need a few of these. Might as well make it a function.
    function create_observable_data_map(
        include_list::Union{Nothing,AbstractVector{ObservableType}},
        size::Int...,
    )
        number_type = Float64
        storage_type = typeof(zeros(number_type, size))
        map = Dict{String,storage_type}()
        for (observable_name, observable) in observables
            if include_list === nothing || observable.type in include_list
                map[observable_name] = zeros(number_type, size)
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
    energy_data = zeros(Float64, num_computed_states)  # Energy for each state (used for energy-dependent observables)
    observable_data =
        create_observable_data_map([ObservableType_Direct], num_computed_states)

    @info "Computing Hamiltonian blocks and observables..."
    # The number of fermions and the color configuration are conserved over tunneling,
    # so we can break the Hamiltonian into blocks labeled by these two quantities
    for (config_idx, (N_fermions, color_configuration)) in enumerate(system_configurations)
        # Size of the Hamiltonian block
        L = block_sizes[config_idx]
        H = SymmetricMatrix(L)  # Use custom "SymmetricMatrix" type to save memory at the cost of speed
        observables_basis = create_observable_data_map([ObservableType_Direct], L)  # Compute the observables for each state as we build the matrix

        # Compute Hamiltonian matrix elements between all pairs of states
        # enumerate_multistate returns elements in a consistent order, so
        # as long as we're consistent, the matrix elements will be in the right place
        # state_i and state_j are arrays of integers, where each integer is a bitmask
        # representing the occupation of each site for a given color
        @threads :greedy for (i, state_i) in enumerate(enumerate_multistate(num_sites, color_configuration))
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
            for observable_name in keys(observables_basis)
                # Pre-compute the observable for this basis state
                observables_basis[observable_name][i] =
                    observables[observable_name].func(state_i)
            end
        end

        num_permutations = num_configuration_permutations(color_configuration)

        @debug begin
            "N_fermions=$N_fermions, color_configuration=$(color_configuration), L=$L, num_configuration_permutations=$(num_permutations), H=$H"
        end

        # Hand off thread control to BLAS for diagonalization
        # to make use of the substantially faster multi-threaded diagonalization routines
        LinearAlgebra.BLAS.set_num_threads(nthreads())

        # Diagonalize the Hamiltonian block
        H_symmetric_view = LinearAlgebra.Symmetric(H, :U)
        # Annoyingly, eigen() forces us to store all the eigenvectors
        # in memory at once, but I can't find a good way around this
        # Even the builtin `eigvals`/`eigvecs` functions are just
        # wrappers around this.
        eigen_data = LinearAlgebra.eigen(H_symmetric_view)

        # Reset BLAS to single-threaded mode
        LinearAlgebra.BLAS.set_num_threads(1)

        @debug begin
            msg = "observable_basis_data:\n"
            for (name, data) in observables_basis
                msg *= "  $name: $data\n"
            end
            msg
        end

        # Compute and store observables for each eigen-state
        offset = size_offset[config_idx]
        @threads for (i, (eigen_val, eigen_vec)) in
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
            energy_data[idx] = eigen_val

            # Compute each observable for this state
            for (observable_name, observable_basis_data) in observables_basis
                # Because we already computed the observables for each basis state,
                # we can just do a weighted sum over those based on the eigenvector components
                observable_value = sum(@. observable_basis_data * eigen_vec * eigen_vec)
                observable_data[observable_name][idx] = observable_value
            end
        end

        @info "Computed data for block $config_idx/$(length(system_configurations)) at $(Dates.now())"
    end

    @info "Computed data for $num_computed_states states."

    @info "Computing derived observables..."

    # First, update the energy so derived observables can use it.
    # The free energy is H + u * N. At each observed energy state, N is const,
    # we can just add u_test * N to each eigenvalue to get the corrected energy.
    energy_data .+= u_test * n_fermion_data
    @debug begin
        "  Updated Energy data: $energy_data"
    end

    for (observable_name, observable) in observables
        if observable.type == ObservableType_Derived
            observable_data[observable_name] =
                observable.func(getindex.((observable_data,), observable.args)...)
        end
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

    # The value of mu such that the density is N/2
    u_shift = (U/2) * (num_colors - 1)

    # Calculate the ground-state energy
    E0 = minimum(@. energy_data - (u_shift * n_fermion_data))

    @info "E0: $E0"

    # Create a new container to store the observable values at each u
    computed_observable_values = create_observable_data_map(nothing, num_temps, num_us)
    @threads for (i, u) in collect(enumerate(u_vals))
        # Shift observables so that density=N/2 at u=0
        u_shifted = u + u_shift

        # Re-weight the data according to the new u value
        # To explain the negatives: The partition function terms are of the form e^(-B * E),
        # The energy contains a term of the form -u * N. We seek to remove the contribution of
        # u_test and add the contribution of u_shifted, so we need E -> E + u_test * N - u_shifted * N
        # (again, because the u * N term is negative). We can rewrite this as E -> E - (u_shifted - u_test) * N.
        # Thus, we need to multiply by e^(-B * -(u_shifted - u_test) * N) = e^(B * (u_shifted - u_test) * N).
        # Because of the way broadcasting works, we have to construct the exponents first
        # (or Julia gets confused). The below syntax makes a matrix of all combinations of
        # elements from B and n_fermion_data
        correction_exponents = permutedims(B * (u_shifted - u_test)) .* n_fermion_data
        # This makes corrected_weights a matrix with indexing [state, temp]
        corrected_weights = @. exp(correction_exponents) * weights

        # Sum over all the values in each column to get a vector of partition functions for each temperature
        Z = sum(corrected_weights, dims = 1)

        @debug begin
            "u=$u, corrected_weights=$corrected_weights, Z=$Z"
        end

        # The intrinsic energy is also used a lot of these, so just precompute it here.
        H = energy_data - (u_shifted * n_fermion_data)

        # Compute any observables that depend on the energy
        energy_dependent_observable_values = create_observable_data_map(
            [ObservableType_Energy_Dependent],
            num_computed_states,
        )
        for observable_name in keys(energy_dependent_observable_values)
            observable = observables[observable_name]
            # Get the observable value for each state
            observable_values = observable.func(u_shifted, energy_data, H, n_fermion_data)
            energy_dependent_observable_values[observable_name] = observable_values
        end

        # Compute the expectation value of each observable
        for (observable_name, observable_values) in
            merge(observable_data, energy_dependent_observable_values)
            expectation_values = sum(corrected_weights .* observable_values; dims = 1)
            normalized_expectation_values = expectation_values ./ Z

            @debug begin
                "  $observable_name: $normalized_expectation_values"
            end

            computed_observable_values[observable_name][:, i] .=
                normalized_expectation_values'
        end

        # Finally, compute the remaining types observables (those that depend on the expectation values or are overlays)
        for (observable_name, observable) in observables
            if observable.type == ObservableType_Expectation_Derived
                # Get the arguments for this observable from the computed_observable_values
                args = getindex.((computed_observable_values,), observable.args)
                # Compute the observable value at this u for each temperature
                computed_observable_values[observable_name][:, i] .=
                    observable.func(u_shifted, B, Z, view.(args, :, i)...)
            elseif observable.type == ObservableType_Overlay
                computed_observable_values[observable_name][:, i] =
                    observable.func.(B, u_shifted)
            end
        end
    end

    return computed_observable_values
end

end
