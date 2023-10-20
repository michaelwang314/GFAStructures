abstract type Integrator end

mutable struct GradientDescent <: Integrator
    subunits::Vector{RigidSubunit}

    step_size::Float64

    ampl::Float64
    Δampl::Float64
    quench_duration::Int64
end

function GradientDescent(subunits::Vector{RigidSubunit}, step_size::Float64, ampl_start::Float64, quench_duration::Int64)
    Δampl = ampl_start / quench_duration

    return GradientDescent(subunits, step_size, ampl_start, Δampl, quench_duration)
end
function GradientDescent(subunits::Vector{RigidSubunit}, step_size::Float64)
    return GradientDescent(subunits, step_size, 0.0, 0.0, 0)
end

function update_subunits!(integrator::GradientDescent)
    Threads.@threads for subunit in integrator.subunits        
        fx, fy, fz = 0.0, 0.0, 0.0
        τx, τy, τz = 0.0, 0.0, 0.0
        for site in subunit.interaction_sites            
            fx += site.force[1] + (integrator.ampl > 0.0 ? integrator.ampl * randn() : 0.0)
            fy += site.force[2] + (integrator.ampl > 0.0 ? integrator.ampl * randn() : 0.0)
            fz += site.force[3] + (integrator.ampl > 0.0 ? integrator.ampl * randn() : 0.0)
            
            rx, ry, rz = site.position[1] - subunit.position[1], site.position[2] - subunit.position[2], site.position[3] - subunit.position[3]
            τx += ry * site.force[3] - site.force[2] * rz 
            τy += -(rx * site.force[3] - site.force[1] * rz) 
            τz += rx * site.force[2] - site.force[1] * ry

            fill!(site.force, 0.0)
            site.energy = 0.0
        end

        translate!(subunit, integrator.step_size * fx, integrator.step_size * fy, integrator.step_size * fz)
        τ = sqrt(τx^2 + τy^2 + τz^2)
        if τ > 0.0
            rotate!(subunit, τx / τ, τy / τ, τz / τ, integrator.step_size * τ)
        end
    end
    integrator.ampl = (integrator.ampl > 0.0 ? integrator.ampl - integrator.Δampl : 0.0)
end