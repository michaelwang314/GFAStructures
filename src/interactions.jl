abstract type Interaction end

struct HarmonicBond <: Interaction
    k::Float64
    r0::Float64
    neighbor_list::FixedPairList
end

function compute_forces!(hb::HarmonicBond)
    Threads.@threads for (site, neighbors) in hb.neighbor_map
        for neighbor in neighbors
            Δx, Δy, Δz = site.position[1] - neighbor.position[1], site.position[2] - neighbor.position[2], site.position[3] - neighbor.position[3]
            if hb.r0 == 0.0
                coef = -hb.k
                if !site.exclude
                    site.energy += 0.5 * hb.k * (Δx^2 + Δy^2 + Δz^2)
                end
            else
                coef = -hb.k * (1.0 - hb.r0 / sqrt(Δx^2 + Δy^2 + Δz^2))
                if !site.exclude
                    site.energy += 0.5 * hb.k * (sqrt(Δx^2 + Δy^2 + Δz^2) - hb.r0)^2
                end
            end

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
    Threads.@threads for (site, neighbors) in lj.neighbor_map
        for neighbor in neighbors
            Δx, Δy, Δz = site.position[1] - neighbor.position[1], site.position[2] - neighbor.position[2], site.position[3] - neighbor.position[3]
            Δr² = Δx^2 + Δy^2 + Δz^2

            if Δr² < lj.cutoff^2
                invΔr⁶ = (lj.σ^2 / Δr²)^3
                coef = lj.ϵ * (48.0 * invΔr⁶ - 24.0) * invΔr⁶ / Δr²

                site.force[1] += coef * Δx
                site.force[2] += coef * Δy
                site.force[3] += coef * Δz

                if !site.exclude
                    #add lj energy.  currently not needed
                end
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

end

function compute_forces!(hs::HertzianSphere{CellList})

end