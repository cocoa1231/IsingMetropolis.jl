# Main Metropolis Algorithm
"""
    Evolve the given lattice `steps` Monte Carlo steps at inverse temperature β.
    TODO: Add support for additional B field interaction
"""
function metropolis!(lattice::AbstractMetropolisLattice, steps::Integer, β::Float64; progressbar = false)

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