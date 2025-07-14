abstract type NeighborList end

struct CellList <: NeighborList
    interaction_sites::Vector{InteractionSite}
    neighbors::Vector{InteractionSite}
    start_id::Array{Int64, 2}
    next_id::Vector{Int64}
    interaction_matrix::Matrix{Bool}
end

struct FixedPairList <: NeighborList
    neighbor_map::Vector{Tuple{InteractionSite, Vector{InteractionSite}}}
end

function FixedPairList(subunits::Vector{RigidSubunit}, cutoff::Float64, interaction_matrix::Matrix{Bool}, site_ids::Vector{Int64}, neighbor_ids::Vector{Int64})
    neighbor_map = find_neighbors(subunits, cutoff, interaction_matrix)
    sorted_neighbor_map = sort_by_id(neighbor_map, site_ids, neighbor_ids)
    return FixedPairList(sorted_neighbor_map)
end

function create_interaction_matrix(num_ids::Int64, pairs::Vector{Tuple{Int64, Int64}})
    interaction_matrix = Matrix{Bool}(undef, num_ids, num_ids)
    for pair in pairs
        interaction_matrix[pair[1], pair[2]] = true
    end
    return interaction_matrix
end

function find_neighbors(subunits::Vector{RigidSubunit}, cutoff::Float64, interaction_matrix::Matrix{Bool})
    neighbor_map = Vector{Tuple{InteractionSite, Vector{InteractionSite}}}()

    N = length(subunits)
    for i = 1 : N, site_i in subunits[i].interaction_sites
        temp_neighbor_list = Vector{InteractionSite}()
        for j = 1 : N, site_j in subunits[j].interaction_sites
            if i != j && interaction_matrix[site_i.id, site_j.id] && sum((site_i.position .- site_j.position).^2) <= cutoff^2
                push!(temp_neighbor_list, site_j)
            end
        end

        if !isempty(temp_neighbor_list)
            push!(neighbor_map, (site_i, temp_neighbor_list))
        end
    end

    return neighbor_map
end

function sort_by_id(neighbor_map::Vector{Tuple{InteractionSite, Vector{InteractionSite}}}, site_ids::Vector{Int64}, neighbor_ids::Vector{Int64})
    sorted_neighbor_map = Vector{Tuple{InteractionSite, Vector{InteractionSite}}}()
    for (site, neighbors) in neighbor_map
        temp_neighbor_list = Vector{InteractionSite}()
        if site.id in site_ids
            for neighbor in neighbors
                if neighbor.id in neighbor_ids
                    push!(temp_neighbor_list, neighbor)
                end
            end
        end
        if !isempty(temp_neighbor_list)
            push!(sorted_neighbor_map, (site, temp_neighbor_list))
        end
    end
    return sorted_neighbor_map
end

function get_neighbor_statistics(neighbor_map::Vector{Tuple{InteractionSite, Vector{InteractionSite}}})
    counts = Dict{Int64, Int64}()
    for (_, neighbors) in neighbor_map
        size = length(neighbors)
        if haskey(counts, size)
            counts[size] += 1
        else
            counts[size] = 1
        end
    end

    println("Getting neighbor statistics...")
    println("There are $(length(neighbor_map)) interaction sites with neighbors:")
    for (size, count) in counts
        println("    $count interaction site(s) with $size neighbor(s)")
    end
end