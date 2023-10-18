abstract type ExternalForce end

struct AnisotropicSpring <: ExternalForce
    sites::Vector{InteractionSite}

    kx::Float64
    ky::Float64
    kz::Float64

    x0::Float64
    y0::Float64
    z0::Float64
end

function compute_forces!(as::AnisotropicSpring)
    for site in as.sites
        site.force[1] += -as.kx * (site.position[1] - as.x0)
        site.force[2] += -as.ky * (site.position[2] - as.y0)
        site.force[3] += -as.kz * (site.position[3] - as.z0)
    end
end

struct ConstantForce <: ExternalForce
    sites::Vector{InteractionSite}

    force::MVector{3, Float64}
end

function compute_forces!(cf::ConstantForce)
    for site in cf.sites
        site.force .+= cf.force
    end
end