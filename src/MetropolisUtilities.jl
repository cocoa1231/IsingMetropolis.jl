"""
    Return the energy associated with a 2D lattice represented by an array.
    Does not assume boundary conditions imposed on lattice structure.
    Uses `scipy.signal.convolve2d` to perform energy calculation with boundary
    condition "wrap" for circular boundary. mode = "same" ensures the output
    array is the same size as the input lattice.

    TODO: Add support for additional B Field interaction
"""
function energy(lattice::AbstractMatrix{Int64})::Int64
    kern = [
        0 1 0;
        1 0 1;
        0 1 0
    ]

    atom_energy = signal.convolve2d(lattice, kern, mode = "same", boundary = "wrap") .* lattice
    return -Int(sum(atom_energy) / 2)
end

"""
    Return the energy required to flip a spin on the lattice specified by a 2-tuple.
    Requires circular boundary conditions implimented by CircularArray.

    TODO: Add support for additional B Field interaction
"""
function dE_at_site(lattice::CircularArray, site)
    x, y = site
    spin_site = lattice[x, y]
    nn_sum = sum([lattice[x+1,y], lattice[x-1, y], lattice[x, y+1], lattice[x, y-1]])
    return 2*spin_site*nn_sum
end


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