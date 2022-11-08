# Main Metropolis Algorithm
"""
    Evolve the given lattice `steps` Monte Carlo steps at inverse temperature β.
"""
function metropolis!(lattice::AbstractMetropolisLattice, steps::Integer, β::Float64; progressbar = false)

    # Create a list of all the exponentials beforehand
    possible_dE = [-8:2:8;]
    exponentials = Dict(possible_dE .=> exp.(-possible_dE .* β))
    N = size(lattice.initial)[1]

    # The field can add -2B or 2B units of energy for arbitrary B
    # We cache the two exponentials that we may have to calculate
    # Along with exp(-β*0) in case field is absent.
    possible_field_dE = [-1, 0, 1]
    field_exp = Dict(possible_field_dE .=> exp.(-2β .* possible_field_dE)) 

    p = Progress(steps)
    for _ in 1:steps

        # Pick Random Point on lattice
        x, y = rand(1:N, 2)

        # Calculate dE
        dE, field_dE = dE_at_site(lattice.final, (x, y); Bfield_strength = lattice.external_field)

        # If new energy is lower accept
        if dE + field_dE < 0
            lattice.final[x, y] *= -1
            push!(lattice.spinsflipped, (x, y))
        elseif dE >= 0
            # Otherwise probabilistically accept w/ probability exp(-β*dE)
            u = rand()
            s = lattice.final[x, y]
            if u < exponentials[dE]*field_exp[Int(s*sign(lattice.external_field))]
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