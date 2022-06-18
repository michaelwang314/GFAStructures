abstract type Interaction end

struct HarmonicBond <: Interaction
    k::Float64
    r0::Float64
    neighbor_list::FixedPairList
end

function compute_forces!(hb::HarmonicBond)
    
end

struct LennardJones{NL <: NeighborList} <: Interaction
    ϵ::Float64
    σ::Float64
    cutoff::Float64
    neighbor_list::NL
end

function compute_forces!(lj::LennardJones{FixedPairList})

end

function compute_forces!(lj::LennardJones{CellList})

end