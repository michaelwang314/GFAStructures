struct System
    subunits::Vector{RigidSubunit}
    interactions::Vector{<:Interaction}
    integrator::Integrator
end

function initialize_lattice(unit_cell::Vector{RigidSubunit}, lattice_vectors::NTuple{3, Vector{Float64}}, dims::NTuple{3, Int64})
    subunits = Vector{RigidSubunit}()
    a1, a2, a3 = lattice_vectors
    for i = 0 : dims[1] - 1, j = 0 : dims[2] - 1, k = 0 : dims[3] - 1
        for subunit in deepcopy(unit_cell)
            translate!(subunit, i .* a1 .+ j .* a2 .+ k .* a3)
            push!(subunits, subunit)
        end
    end
    return subunits
end
function initialize_lattice(unit_cell::Vector{RigidSubunit}, lattice_vectors::NTuple{2, Vector{Float64}}, dims::NTuple{2, Int64})
    return initialize_lattice(unit_cell, (lattice_vectors[1], lattice_vectors[2], [0.0, 0.0, 1.0]), (dims[1], dims[2], 1))
end

function find_neighbors(subunits::Vector{RigidSubunit}, subunit_cutoff::Float64, neighbor_cutoff::Float64, interaction_matrix::Matrix{Bool})
    neighbor_map = Vector{Tuple{InteractionSite, Vector{InteractionSite}}}()

    N = length(subunits)
    for i = 1 : N, site_i in subunits[i].interaction_sites
        temp_neighbor_list = Vector{InteractionSite}()
        for j = 1 : N, site_j in subunits[j].interaction_sites
            if i != j && interaction_matrix[site_i.id, site_j.id] && sum((subunits[i].position .- subunits[j].position).^2) <= subunit_cutoff^2 && sum((site_i.position .- site_j.position).^2) <= neighbor_cutoff^2
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
        for neighbor in neighbors
            if site.id in site_ids && neighbor.id in neighbor_ids
                push!(temp_neighbor_list, neighbor)
            end
        end
        if !isempty(temp_neighbor_list)
            push!(sorted_neighbor_map, (interaction_site, temp_neighbor_list))
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

function get_energy(subunits::Vector{RigidSubunit})
    energy = 0.0
    for subunit in subunits
        energy += subunit.energy
    end
    return energy / 2
end

function hr_min_sec(time::Float64)
    hours = trunc(Int64, time / 3600.0)
    minutes = trunc(Int64, mod(time, 3600.0) / 60.0)
    seconds = trunc(Int64, mod(time, 60.0))

    return string(hours < 10 ? "0" : "", hours, 
                  minutes < 10 ? ":0" : ":", minutes, 
                  seconds < 10 ? ":0" : ":", seconds)
end

function run_simulation!(system::System; num_steps::Int64 = 1, message_interval::Float64 = 10.0)
    println("Simulation started.............................................")
    println("Number of subunits: ", length(system.subunits))
    
    prev_step = 0
    time_elapsed = 0.0
    interval_start = time()
    for step = 1 : num_steps
        for interaction in system.interactions
            compute_forces!(interaction)
        end
        update_subunits!(system.integrator)

        interval_time = time() - interval_start
        if interval_time > message_interval || step == num_steps
            time_elapsed += interval_time
            rate = (step - prev_step) / interval_time
            println(hr_min_sec(time_elapsed), " | ",
                    step, "/", num_steps, " (", round(step / num_steps * 100, digits = 1), "%) | ",
                    round(rate, digits = 1), " steps/s | ",
                    hr_min_sec((num_steps - step) / rate), " | ", 
                    "energy = ", get_energy(system.subunits))
            prev_step = step
            interval_start = time()
        end
    end
    println("Average steps/s: ", round(num_steps / time_elapsed, digits = 1))
    println("Simulation done................................................")
end

function format_for_mathematica(system::System, file::String; params = [])
    if !isdir(dirname(file))
        mkpath(dirname(file))
    end

    param_str = ""
    for p in params
        param_str *= "$p,"
    end
    param_str = "{" * chop(param_str) * "}"

    subunit_data = "{"
    for subunit in system.subunits
        x, y, z = subunit.position
        b1, b2 = subunit.body_axes
        energy = subunit.energy
        subunit_data *= "{{$x, $y, $z}, {$(b1[1]), $(b1[2]), $(b1[3])}, {$(b2[1]), $(b2[2]), $(b2[3])}, $energy, $param_str},"
    end
    subunit_data = replace(chop(subunit_data) * "}", "e" => "*^")

    open(file, "w") do io
        write(io, subunit_data)
    end
end

function save!(system::System, file::String)
    if !isdir(dirname(file))
        mkpath(dirname(file))
    end

    open(file, "w") do io
        serialize(io, system)
    end
end

function load(file::String)
    return begin
        open(file, "r") do io
            deserialize(io)
        end
    end
end