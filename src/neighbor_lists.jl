abstract type NeighborList end

struct CellList <: NeighborList
    neighbor_map::Vector{Tuple{InteractionSite, Vector{InteractionSite}}}
    start_id::Array{Int64, 2}
    next_id::Vector{Int64}
    interaction_matrix::Matrix{Bool}
end

function update_neighbor_list!(cell_list::CellList)

end

struct FixedPairList <: NeighborList
    neighbor_map::Vector{Tuple{InteractionSite, Vector{InteractionSite}}}
end