"""
    Return the energy associated with a 2D lattice represented by an array.
    Does not assume boundary conditions imposed on lattice structure.
    Uses `scipy.signal.convolve2d` to perform energy calculation with boundary
    condition "wrap".
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
"""
function dE_at_site(lattice::CircularArray, site)
    x, y = site
    spin_site = lattice[x, y]
    nn_sum = sum([lattice[x+1,y], lattice[x-1, y], lattice[x, y+1], lattice[x, y-1]])
    return 2*spin_site*nn_sum
end

# Main Metropolis Algorithm
"""
    Evolve the given lattice `steps` Monte Carlo steps at inverse temperature β.
"""
function metropolis!(lattice::AbstractMetropolisLattice, steps::Integer, β::Float64; progressbar = true)

    # Create a list of all the exponentials beforehand
    possible_dE = [-8:2:8;]
    exponentials = Dict(possible_dE .=> exp.(-possible_dE .* β))
    N = size(lattice.initial)[1]

    p = Progress(steps)
    for _ in 1:steps

        # Pick Random Point on lattice
        x, y = rand(1:N, 2)

        # Calculate dE
        dE = dE_at_site(lattice.final, (x, y))

        # If new energy is lower accept
        if dE < 0
            lattice.final[x, y] *= -1
            push!(lattice.spinsflipped, (x, y))
        elseif dE >= 0
            # Otherwise probabilistically accept w/ probability exp(-β*dE)
            u = rand()
            if u < exponentials[dE]
                lattice.final[x, y] *= -1
                push!(lattice.spinsflipped, (x, y))
            else
                push!(lattice.spinsflipped, (-1, -1))
            end
        end

        if progressbar == true
            next!(p)
        end
    end

    return lattice
end