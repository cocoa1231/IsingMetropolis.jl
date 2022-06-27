module IsingMetropolis

using CircularArrays
using SciPy:signal as signal
using ProgressMeter

# Includes
include("MetropolisDataStructures.jl")
include("MetropolisAlgorithm.jl")

export generate_circular_lattice
export MetropolisLattice
export energy
export dE_at_site
export metropolis!

"""
    Generates a 2D square lattice of side N with circular boundary
    conditions. Boundary conditions are imposed by using the type
    `CircularArray` from `CircularArrays.jl`
"""
function generate_circular_lattice(N::Integer, T₀::Symbol)
    if T₀ == :zero
        return CircularArray(ones(Int, (N, N)))
    elseif T₀ == :infty
        return CircularArray( rand([-1, 1], (N, N)) )
    end
end
end # module IsingMetropolis