export System

mutable struct System
    subunits::Vector{Subunit}
    bonds::Vector{Bond}

    integrator::TemperatureQuench
end

function run_simulation!(system::System; num_steps::Int64 = 1, message_interval::Float64 = 10.0)

end