module IsingMetropolis

using CircularArrays
using SciPy:signal as signal
using ProgressMeter

# Includes
include("MetropolisDataStructures.jl")
include("MetropolisUtilities.jl")
include("MetropolisAlgorithm.jl")

# Exports
export generate_circular_lattice
export MetropolisLattice
export energy
export dE_at_site
export metropolis!
export fill_data

# Methods on the main AbstractMetropolisLattice type go here.
# Other methods go in MetropolisUtilities.jl

"""
    Fill magnetization history of the lattice based on the tracked spins that
    were flipped over the course of execution of the metropolis algorithm.
"""
function _fill_M_history!(lattice::AbstractMetropolisLattice)
    l = copy(lattice.initial)

    # Clear magnetization history if there was any
    lattice.magnetization_hist = Float64[sum(l)]
    for s_k in lattice.spinsflipped
        if s_k == (-1, -1)
            # If no spin was flipped, magnetization has not changed from last state
            push!(lattice.magnetization_hist, lattice.magnetization_hist[end])
        else
            # Update lattice for the next calculation
            l[s_k...] *= -1

            # Calculation new magnetization and append
            push!(lattice.magnetization_hist, lattice.magnetization_hist[end] + 2*l[s_k...])
        end
    end
    return l
end

"""
    Fill internal energy history of the lattice based on the tracked spins that
    were flipped over the course of execution of the metropolis algorithm.
"""
function _fill_U_history!(lattice::AbstractMetropolisLattice)
    l = copy(lattice.initial)
    lattice.internalenergy_hist = Float64[sum(energy(l))]
    for s_k in lattice.spinsflipped
        if s_k == (-1, -1)
            # If no spin was flipped then energy has not changed
            push!(lattice.internalenergy_hist, lattice.internalenergy_hist[end])
        else
            # Calculate dE at the site and append it
            push!(lattice.internalenergy_hist, lattice.internalenergy_hist[end] + dE_at_site(l, s_k))

            # Update lattice for next calculation
            l[s_k...] *= -1
        end
    end
    return l
end

"""
    Fill the internal energy or magnetization history of a lattice evolved using the
    `IsingMetropolis.metropolis!` function.
"""
function fill_data(lattice::AbstractMetropolisLattice, data::Symbol)
    if data == :M
        return _fill_M_history!(lattice)
    elseif data == :U
        return _fill_U_history!(lattice)
    else
        throw(ArgumentError("Data parameter $(string(data)) not implimented!"))
    end
end

end # module IsingMetropolis