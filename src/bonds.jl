export BindingSite
export Interaction
export Harmonic
export Bond
export compute_forces!

mutable struct BindingSite
    position::SVector{3, Float64}

    energy::Float64
end

abstract type Interaction end
mutable struct Harmonic <: Interaction
    k::Float64
end
mutable struct Cosine <: Interaction

end

mutable struct Bond{I <: Interaction}
    pair::Tuple{BindingSite, BindingSite}

    interaction::I
end

function compute_forces!(bonds::Vector{Bond{Harmonic}})
    Threads.@threads for bond in bonds
        site1, site2 = bond.pair
        position1, position2 = site1.position, site2.position
        Δx, Δy, Δz = position1[1] - position2[1], position1[2] - position2[2], position1[3] - position2[3]

        f1, f2 = site1.force, site2.force
        k = bond.interaction.k
        f2[1] = -(f1[1] = -k * Δx)
        f2[2] = -(f1[2] = -k * Δy)
        f2[3] = -(f1[3] = -k * Δz)

        site1.energy = site2.energy = 0.5 * k * (Δx^2 + Δy^2 + Δz^2)
    end
end

function compute_forces!(bonds::Vector{Bond{Cosine}})

end