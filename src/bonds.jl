export BindingSite
export Interaction
export Harmonic
export Bond
export compute_forces!

mutable struct BindingSite
    position::SVector{3, Float64}
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
        s1, s2 = bond.pair
        r1, r2 = s1.position, s2.position
        Δx, Δy, Δz = r1[1] - r2[1], r1[2] - r2[2], r1[3] - r2[3]

        f1, f2 = s1.force, s2.force
        k = bond.interaction.k
        f2[1] = -(f1[1] = -k * Δx)
        f2[2] = -(f1[2] = -k * Δy)
        f2[3] = -(f1[3] = -k * Δz)

        s1.energy = s2.energy = 0.5 * k * (Δx^2 + Δy^2 + Δz^2)
    end
end

function compute_forces!(bonds::Vector{Bond{Cosine}})

end