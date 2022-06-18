struct System
    subunits::Vector{Subunit}
    interactions::Vector{<:Interaction}
    integrator
end

function get_energy(subunit::Subunit)
    energy = 0.0
    for site in subunit.interaction_sites
        if !site.exclude
            energy += site.energy
        end
    end
    return energy
end
function get_energy(subunits::Vector{Subunit})
    energy = 0.0
    for subunit in subunits
        energy += get_energy(subunit)
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