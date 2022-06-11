export System
export initialize_lattice
export link
export run_simulation!
export format_for_mathematica!

mutable struct System
    subunits::Vector{Subunit}
    bonds::Vector{Bond}

    integrator::TemperatureQuench
end

function initialize_lattice(unit_cell::Vector{Subunit}, lattice_vectors::Vector{Vector{Float64}}, lattice_dimensions::Vector{Int64})
    subunits = Vector{Subunit}()

    a1, a2, a3 = lattice_vectors
    Δr = MVector(0.0, 0.0, 0.0)
    for i = 0 : lattice_dimensions[1] - 1, j = 0 : lattice_dimensions[2] - 1, k = 0 : lattice_dimensions[3] - 1
        for subunit in deepcopy(unit_cell)
            @. Δr = i * a1 + j * a2 + k * a3
            translate!(subunit, Δr)
            push!(subunits, subunit)
        end
    end

    return subunits
end

function link(subunits::Vector{Subunit}, interactions::Vector{Tuple{Int64, Int64, I}}, neighbor_cutoff::Float64, bond_cutoff::Float64) where I <: Interaction
    bonds = Vector{Bond}()
    
    N_sub = length(subunits)
    for i = 1 : N_sub - 1, j = i + 1 : N_sub
        sub1, sub2 = subunits[i], subunits[j]
        sub1_pos, sub2_pos = sub_i.position, sub_j.position
        if (sub1_pos[1] - sub2_pos[1])^2 + (sub1_pos[2] - sub2_pos[2])^2 + (sub1_pos[3] - sub2_pos[3])^2 <= neighbor_cutoff^2
            for (id1, id2, interaction) in interactions
                site1, site2 = sub1.binding_sites[id1], sub2.binding_sites[id2]
                site1_pos, site2_pos = site1.position, site2.position
                if (site1_pos[1] - site2_pos[1])^2 + (site1_pos[2] - site2_pos[2])^2 + (site1_pos[3] - site2_pos[3])^2 <= bond_cutoff^2
                    push!(bonds, Bond((site1, site2), interaction))
                end

                if id1 != id2
                    site1, site2 = sub1.binding_sites[id2], sub2.binding_sites[id1]
                    site1_pos, site2_pos = site1.position, site2.position
                    if (site1_pos[1] - site2_pos[1])^2 + (site1_pos[2] - site2_pos[2])^2 + (site1_pos[3] - site2_pos[3])^2 <= bond_cutoff^2
                        push!(bonds, Bond((site1, site2), interaction))
                    end
                end
            end
        end
    end

    return bonds
end

function get_energy(subunit::Subunit)
    energy = 0.0
    for site in subunit.binding_sites
        energy += site.energy
    end
    return energy
end

function get_energy(subunits::Vector{Subunit})
    energy = 0.0
    for subunit in subunits
        energy += get_energy(subunit)
    end
    return energy
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

    prev_step = 0
    time_elapsed = 0.0
    interval_start = time()
    for step = 1 : num_steps
        compute_forces!(system.bonds)
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
end

function format_for_mathematica!(system::System; params = [], file::String = "TEMP.txt")
    if !isdir(dirname(file))
        mkpath(dirname(file))
    end


end