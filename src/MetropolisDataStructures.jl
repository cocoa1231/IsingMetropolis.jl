import Base:size

abstract type AbstractMetropolisLattice end

mutable struct MetropolisLattice <: AbstractMetropolisLattice
    initial::AbstractMatrix{Int64}
    final::AbstractMatrix{Int64}
    external_field::Float64
    spinsflipped::Vector{Tuple{Int64, Int64}}
    internalenergy_hist::Vector{Float64}
    magnetization_hist::Vector{Float64}
end

# External constructors below for with or without external field

function MetropolisLattice(lattice::AbstractArray{Int64, 2})
    MetropolisLattice(
        copy(lattice),
        copy(lattice),
        0,
        Vector{Tuple{Int64, Int64}}(),
        Float64[],
        Float64[]
    )
end

function MetropolisLattice(lattice::AbstractArray{Int64, 2}, external_field)
    MetropolisLattice(
        copy(lattice),
        copy(lattice),
        external_field,
        Vector{Tuple{Int64, Int64}}(),
        Float64[],
        Float64[]
    )
end

size(l::MetropolisLattice) = size(l.initial)