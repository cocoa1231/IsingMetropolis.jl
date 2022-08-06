using IsingMetropolis
using Test
using CircularArrays
using Statistics

@testset "MetropolisUtilities.jl" begin
    # Test if a lattice is generated correctly
    N = 1000
    zero_temp_lattice = CircularArray(ones(Int, (N, N)))
    inf_temp_lattice = CircularArray(rand([-1, 1], (N, N)))
    @test generate_circular_lattice(1000, :zero) == zero_temp_lattice
    @test isapprox( mean(generate_circular_lattice(N, :infty)), 
                    mean(inf_temp_lattice), atol = 0.01 )

    # Test energy calculations with and without field
    Bfield = -20
    zero_temp_field_lattice = MetropolisLattice(generate_circular_lattice(N, :zero), Bfield)
    N = prod(size(zero_temp_lattice))
    @test energy(zero_temp_lattice) == (-2*N, 0)
    @test energy(zero_temp_field_lattice.initial, Bfield_strength = Bfield) == 
        (-2*N, -Bfield * sum(zero_temp_field_lattice.initial))

    # Test if ΔE is calculated properly without field
    new_l = copy(zero_temp_lattice)
    new_l[10, 10] *= -1
    @test energy(new_l) .- energy(zero_temp_lattice) == dE_at_site(zero_temp_lattice, (10, 10))

    # Test if ΔE is calculated properly with field
    new_l = copy(zero_temp_field_lattice.initial)
    new_l[10, 10] *= -1
    energy_diff = energy(new_l; Bfield_strength = Bfield) .- energy(zero_temp_field_lattice.initial, Bfield_strength = Bfield)
    @test energy_diff == dE_at_site(zero_temp_field_lattice.initial, (10, 10), Bfield_strength = Bfield)
end

@testset "MetropolisAlgorithm.jl" begin
    # Test Metropolis Algorithm with simple case without field
    N = 32
    zero_temp_lattice = generate_circular_lattice(N, :zero)
    l = MetropolisLattice(zero_temp_lattice)
    β = 100.
    steps = N*N*10000
    @info "Please hold on, simulating a 32x32 lattice for testing..."
    metropolis!(l, steps, β, progressbar = false)
    @test isapprox(energy(l.final)[1], -2*N*N, atol = N)
    @test isapprox(energy(l.final)[2], 0, atol = N)

    # Test Metropolis Algorithm with simple case with strong field
    N = 32
    Bfield = -20.
    zero_temp_lattice = generate_circular_lattice(N, :zero)
    l = MetropolisLattice(zero_temp_lattice, Bfield)
    β = 100.
    steps = N*N*10000
    @info "Please hold on, simulating a 32x32 lattice for testing..."
    metropolis!(l, steps, β, progressbar = false)
    @test isapprox(energy(l.final, Bfield_strength = Bfield)[1], -2*N*N, atol = N)

    # The final state we know should be all spins aligned with the field, so we flip all the spins
    # manually then calculate the field interaction energy and test if the energy of the final state
    # is the energy we expect
    @test isapprox(energy(l.final, Bfield_strength = Bfield)[2], -Bfield * sum(-1 .* l.initial), atol = N)
end

@testset "IsingMetropolis.jl" begin
    # Test if we can fill magnetization and internal energy history accurately
    N = 32
    zero_temp_lattice = generate_circular_lattice(N, :zero)
    l = MetropolisLattice(zero_temp_lattice)
    β = 100.
    steps = N*N*10000
    metropolis!(l, steps, β)

    # Internal energy test
    final_state = fill_data!(l, :U)
    @test l.final == final_state
    @test isapprox(energy(final_state)[1], -2*N*N, atol = N)
    @test isapprox(energy(final_state)[2], 0, atol = N)

    # Magnetization test
    final_state = fill_data!(l, :M)
    @test l.final == final_state
    @test isapprox(abs(sum(l.final)), N*N, atol = Int(N/4)) 
end