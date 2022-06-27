module IsingMetropolis

using CircularArrays
using SciPy:signal as signal
using ProgressMeter

# Includes
include("MetropolisDataStructures.jl")
include("MetropolisUtilities.jl")
include("MetropolisAlgorithm.jl")

export generate_circular_lattice
export MetropolisLattice
export energy
export dE_at_site
export metropolis!

end # module IsingMetropolis