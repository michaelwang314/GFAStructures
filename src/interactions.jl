abstract type Interaction end

struct HarmonicBond <: Interaction
    k::Float64
    r0::Float64
    neighbor_list::FixedPairList
end

function compute_forces!(hb::HarmonicBond)
    Threads.@threads for (i, site) in collect(enumerate(hb.neighbor_list.interaction_sites))
        for neighbor in hb.neighbor_list.neighbors[i]
            Δx, Δy, Δz = site.position[1] - neighbor.position[1], site.position[2] - neighbor.position[2], site.position[3] - neighbor.position[3]
            coef = (hb.r0 == 0 ? -hb.k : -hb.k * (1.0 - hb.r0 / sqrt(Δx^2 + Δy^2 + Δz^2)))

            site.force[1] += coef * Δx
            site.force[2] += coef * Δy
            site.force[3] += coef * Δz
        end
    end
end

struct LennardJones{NL <: NeighborList} <: Interaction
    ϵ::Float64
    σ::Float64
    cutoff::Float64
    neighbor_list::NL
end

function compute_forces!(lj::LennardJones{FixedPairList})
    Threads.@threads for (i, site) in collect(enumerate(lj.neighbor_list.interaction_sites))
        for neighbor in lj.neighbor_list.neighbors[i]
            Δx, Δy, Δz = site.position[1] - neighbor.position[1], site.position[2] - neighbor.position[2], site.position[3] - neighbor.position[3]
            Δr² = Δx^2 + Δy^2 + Δz^2

            if Δr² < lj.cutoff^2
                invΔr⁶ = (lj.σ^2 / Δr²)^3
                coef = lj.ϵ * (48.0 * invΔr⁶ - 24.0) * invΔr⁶ / Δr²

                site.force[1] += coef * Δx
                site.force[2] += coef * Δy
                site.force[3] += coef * Δz
            end
        end
    end
end

function compute_forces!(lj::LennardJones{CellList})

end

struct HertzianSphere{NL <: NeighborList} <: Interaction
    ϵ::Float64
    σ::Float64
    neighbor_list::NL
end

function compute_forces!(hs::HertzianSphere{FixedPairList})
    Threads.@threads for (i, site) in collect(enumerate(hs.neighbor_list.interaction_sites))
        for neighbor in hs.neighbor_list.neighbors[i]
            Δx, Δy, Δz = site.position[1] - neighbor.position[1], site.position[2] - neighbor.position[2], site.position[3] - neighbor.position[3]
            Δr² = Δx^2 + Δy^2 + Δz^2

            if Δr² < hs.σ^2
                coef = hs.ϵ * sqrt(1 - sqrt(Δr²) / hs.σ)^5

                site.force[1] += coef * Δx
                site.force[2] += coef * Δy
                site.force[3] += coef * Δz
            end
        end
    end
end

function compute_forces!(hs::HertzianSphere{CellList})

end