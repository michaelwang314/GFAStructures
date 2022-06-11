export TemperatureQuench

mutable struct TemperatureQuench
    subunits::Vector{Subunit}

    step_size::Float64

    T_start::Float64
    ΔT::Float64
    quench_steps::Int64
end

function update_subunits!(integrator::TemperatureQuench)
    
end