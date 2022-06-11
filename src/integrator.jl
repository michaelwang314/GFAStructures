export TemperatureQuench

mutable struct TemperatureQuench
    subunits::Vector{Subunit}

    step_size::Float64

    T::Float64
    ΔT::Float64
    quench_duration::Int64
end

function TemperatureQuench(subunits::Vector{Subunit}, step_size::Float64, T_start::Float64, quench_duration::Int64)
    ΔT = T_start / quench_duration

    return TemperatureQuench(subunits, step_size, T_start, ΔT, quench_duration)
end

function update_subunits!(integrator::TemperatureQuench)
    δ = integrator.step_size
    Threads.@threads for subunit in integrator.subunits
        sub_pos = subunit.position
        
        fx, fy, fz = 0.0, 0.0, 0.0
        τx, τy, τz = 0.0, 0.0, 0.0
        for site in subunit.binding_sites
            site_pos = site.position
            sf = site.force
            
            fx += sf[1]
            fy += sf[2]
            fz += sf[3]
            
            rx, ry, rz = site_pos[1] - sub_pos[1], site_pos[2] - sub_pos[2], site_pos[3] - sub_pos[3]
            τx += ry * sf[3] - sf[2] * rz 
            τy += -(rx * sf[3] - sf[1] * rz) 
            τz += rx * sf[2] - sf[1] * ry
        end

        translate!(subunit, δ * fx, δ * fy, δ * fz)
        τ = sqrt(τx^2 + τy^2 + τz^2)
        if τ > 0.0
            rotate!(subunit, τx / τ, τy / τ, τz / τ, δ * τ)
        end
    end
end