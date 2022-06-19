mutable struct InteractionSite
    position::MVector{3, Float64}
    force::MVector{3, Float64}

    exclude::Bool
    energy::Float64
    id::Int64
end

struct RigidSubunit
    position::MVector{3, Float64}
    interaction_sites::Vector{InteractionSite}
    body_axes::MVector{2, MVector{3, Float64}}
end

function RigidSubunit(interaction_sites::Vector{InteractionSite})
    position = MVector{3}(0.0, 0.0, 0.0)
    for site in interaction_sites
        position .+= site.position
    end
    position ./= length(interaction_sites)

    return RigidSubunit(position, interaction_sites, ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]))
end

function rotate!(subunit::RigidSubunit, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)
    sin, cos = sincos(θ)
    Rxx, Rxy, Rxz = cos + axis_x^2 * (1 - cos), axis_x * axis_y * (1 - cos) - axis_z * sin, axis_x * axis_z * (1 - cos) + axis_y * sin
    Ryx, Ryy, Ryz = axis_y * axis_x * (1 - cos) + axis_z * sin, cos + axis_y^2 * (1 - cos), axis_y * axis_z * (1 - cos) - axis_x * sin
    Rzx, Rzy, Rzz = axis_z * axis_x * (1 - cos) - axis_y * sin, axis_z * axis_y * (1 - cos) + axis_x * sin, cos + axis_z^2 * (1 - cos)

    for site in subunit.interaction_sites
        rx, ry, rz = site.position[1] - subunit.position[1], site.position[2] - subunit.position[2], site.position[3] - subunit.position[3]

        site.position[1], site.position[2], site.position[3] = begin
            subunit.position[1] + Rxx * rx + Rxy * ry + Rxz * rz,
            subunit.position[2] + Ryx * rx + Ryy * ry + Ryz * rz,
            subunit.position[3] + Rzx * rx + Rzy * ry + Rzz * rz
        end 
    end

    for b in subunit.body_axes
        b[1], b[2], b[3] = Rxx * b[1] + Rxy * b[2] + Rxz * b[3], Ryx * b[1] + Ryy * b[2] + Ryz * b[3], Rzx * b[1] + Rzy * b[2] + Rzz * b[3]
    end
end
rotate!(subunit::RigidSubunit, axis::V, θ::Float64) where V <: AbstractVector = rotate!(subunit, axis[1], axis[2], axis[3], θ)

function translate!(subunit::RigidSubunit, Δx::Float64, Δy::Float64, Δz::Float64)
    subunit.position[1] += Δx
    subunit.position[2] += Δy
    subunit.position[3] += Δz

    for site in subunit.interaction_sites
        site.position[1] += Δx
        site.position[2] += Δy
        site.position[3] += Δz
    end
end
translate!(subunit::RigidSubunit, Δr::V) where V <: AbstractVector = translate!(subunit, Δr[1], Δr[2], Δr[3])

function get_energy(subunit::Subunit)
    energy = 0.0
    for site in subunit.interaction_sites
        if !site.exclude
            energy += site.energy
        end
    end
    return energy
end