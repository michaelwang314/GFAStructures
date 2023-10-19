abstract type ExternalForce end

struct AnisotropicSpring <: ExternalForce
    sites::Vector{InteractionSite}

    k::MVector{3, Float64}
    r0::MVector{3, Float64}
end

function compute_forces!(as::AnisotropicSpring)
    for site in as.sites
        site.force[1] += -as.k[1] * (site.position[1] - as.r0[1])
        site.force[2] += -as.k[2] * (site.position[2] - as.r0[2])
        site.force[3] += -as.k[3] * (site.position[3] - as.r0[3])
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