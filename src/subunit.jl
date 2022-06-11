export Subunit

mutable struct Subunit
    position::SVector{3, Float64}
    binding_sites::Vector{BindingSite}

    force::SVector{3, Float64}
    torque::SVector{3, Float64}
end

function rotate!(subunit::Subunit, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)
    
end
rotate!(subunit::Subunit, axis::Vector{Float64}, θ::Float64) = rotate!(subunit, axis[1], axis[2], axis[3], θ)

function translate!(subunit::Subunit, Δx::Float64, Δy::Float64, Δz::Float64)
    position = subunit.position
    position[1] += Δx
    position[2] += Δy
    position[3] += Δz

    for site in subunit.binding_sites
        position = site.position
        position[1] += Δx
        position[2] += Δy
        position[3] += Δz
    end
end

function translate!(subunit::Subunit, Δr::Vector{Float64})
    subunit.position .+= Δr

    for site in subunit.binding_sites
        site.position .+= Δr
    end
end