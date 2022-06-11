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
    Threads.@threads for subunit in integrator.subunits
        sub_pos = subunit.position
        sub_force = fill!(subunit.force, 0.0)
        sub_torque = fill!(subunit.torque, 0.0)
        
        for site in subunit.binding_sites
            site_pos = site.position
            site_force = site.force
            
            sub_force .+= site_force
            
            rx, ry, rz = site_pos[1] - sub_pos[1], site_pos[2] - sub_pos[2], site_pos[3] - sub_pos[3]
            sub_torque[1] += ry * site_force[3] - site_force[2] * rz 
            sub_torque[2] += -(rx * site_force[3] - site_force[1] * rz) 
            sub_torque[3] += rx * site_force[2] - site_force[1] * ry
        end

        δ = integrator.step_size
        translate!(subunit, δ * sub_force[1], δ * sub_force[2], δ * sub_force[3])

        τ = sqrt(sub_torque[1]^2 + sub_torque[2]^2 + sub_torque[3]^2)
        if τ > 0.0
            rotate!(subunit, sub_torque[1] / τ, sub_torque[2] / τ, sub_torque[3] / τ, δ * τ)
        end
    end
end