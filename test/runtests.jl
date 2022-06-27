using IsingMetropolis
using Test
using CircularArrays
using Statistics

@testset "IsingMetropolis.jl" begin
    # Test if a lattice is generated correctly
    zero_temp_lattice = CircularArray(ones(Int, (1000, 1000)))
    inf_temp_lattice = CircularArray(rand([-1, 1], (1000, 1000)))
    @test generate_circular_lattice(1000, :zero) == zero_temp_lattice
    @test isapprox( mean(generate_circular_lattice(1000, :infty)), 
                    mean(inf_temp_lattice), atol = 0.01 )

    # Test energy calculations
    N = prod(size(zero_temp_lattice))
    @test energy(zero_temp_lattice) == -2*N

    # Test if ΔE is calculated properly
    true_ΔE = 4
    new_l = copy(zero_temp_lattice)
    new_l[10, 10] *= -1
    @test energy(new_l) - energy(zero_temp_lattice) == dE_at_site(zero_temp_lattice, (10, 10))

    # Test Metropolis Algorithm with simple case
    N = 32
    zero_temp_lattice = generate_circular_lattice(N, :zero)
    l = MetropolisLattice(zero_temp_lattice)
    β = 100.
    steps = N*N*10000
    metropolis!(l, steps, β)
    @test isapprox(energy(l.final), -2*N*N, atol = N)
end