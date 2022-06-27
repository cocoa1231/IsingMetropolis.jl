abstract type AbstractMetropolisLattice end

mutable struct MetropolisLattice <: AbstractMetropolisLattice
    initial::AbstractMatrix{Int64}
    final::AbstractMatrix{Int64}
    spinsflipped::Vector{Tuple{Int64, Int64}}
    internalenergy_hist::Vector{Float64}
    magnetization_hist::Vector{Float64}
end

function MetropolisLattice(lattice::AbstractArray{Int64, 2})
    MetropolisLattice(
        copy(lattice),
        copy(lattice),
        Vector{Tuple{Int64, Int64}}(),
        Float64[],
        Float64[]
    )
end