export BindingSite
export Interaction
export Harmonic
export Bond
export LennardJones

mutable struct BindingSite
    position::MVector{3, Float64}
    force::MVector{3, Float64}

    energy::Float64
end

function BindingSite(position::V) where V <: AbstractVector
    return BindingSite(MVector{3}(position), MVector(0.0, 0.0, 0.0), 0.0)
end

abstract type Interaction end
struct Harmonic <: Interaction
    k::Float64
end
struct LennardJones <: Interaction
    ϵ::Float64
    σ::Float64
    cutoff::Float64
end

mutable struct Bond{I <: Interaction}
    pair::Tuple{BindingSite, BindingSite}

    interaction::I
end

function compute_forces!(bond::Bond{Harmonic})
    site1, site2 = bond.pair
    pos1, pos2 = site1.position, site2.position
    Δx, Δy, Δz = pos1[1] - pos2[1], pos1[2] - pos2[2], pos1[3] - pos2[3]

    f1, f2 = site1.force, site2.force
    k = bond.interaction.k
    f2[1] = -(f1[1] = -k * Δx)
    f2[2] = -(f1[2] = -k * Δy)
    f2[3] = -(f1[3] = -k * Δz)

    site1.energy = site2.energy = 0.5 * k * (Δx^2 + Δy^2 + Δz^2)
end

function compute_forces!(bond::Bond{LennardJones})
    site1, site2 = bond.pair
    pos1, pos2 = site1.position, site2.position
    Δx, Δy, Δz = pos1[1] - pos2[1], pos1[2] - pos2[2], pos1[3] - pos2[3]
    Δr² = Δx^2 + Δy^2 + Δz^2

    if Δr² < bond.interaction.cutoff^2
        f1, f2 = site1.force, site2.force

        invΔr² = bond.interaction.σ^2 / Δr^2
        invΔr⁶ = invΔr²^3
        coef = bond.interaction.ϵ * (48.0 * invΔr⁶ - 24.0) * invΔr⁶ / Δr²

        f2[1] = -(f1[1] = coef * Δx)
        f2[2] = -(f1[2] = coef * Δy)
        f2[3] = -(f1[3] = coef * Δz)
    end
end

function compute_forces!(bonds::Vector{Bond})
    Threads.@threads for bond in bonds
        compute_forces!(bond)
    end
end