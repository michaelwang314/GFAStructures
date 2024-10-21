mutable struct InteractionSite
    position::MVector{3, Float64}
    force::MVector{3, Float64}

    exclude::Bool
    energy::Float64
    id::Int64
end

function InteractionSite(position::V, id::Int64; exclude::Bool = false) where V <: AbstractVector
    return InteractionSite(position, MVector{3}(0.0, 0.0, 0.0), exclude, 0.0, id)
end

mutable struct RigidSubunit
    position::MVector{3, Float64}
    interaction_sites::Vector{InteractionSite}
    body_axes::MVector{2, MVector{3, Float64}}

    energy::Float64
end

function RigidSubunit(interaction_sites::Vector{InteractionSite})
    position = MVector{3}(0.0, 0.0, 0.0)
    for site in interaction_sites
        position .+= site.position
    end
    position ./= length(interaction_sites)

    return RigidSubunit(position, interaction_sites, ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]), 0.0)
end

function rotate!(subunit::RigidSubunit, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)
    sin, cos = sincos(θ)
    Rxx, Rxy, Rxz = cos + axis_x^2 * (1 - cos), axis_x * axis_y * (1 - cos) - axis_z * sin, axis_x * axis_z * (1 - cos) + axis_y * sin
    Ryx, Ryy, Ryz = axis_y * axis_x * (1 - cos) + axis_z * sin, cos + axis_y^2 * (1 - cos), axis_y * axis_z * (1 - cos) - axis_x * sin
    Rzx, Rzy, Rzz = axis_z * axis_x * (1 - cos) - axis_y * sin, axis_z * axis_y * (1 - cos) + axis_x * sin, cos + axis_z^2 * (1 - cos)

    for site in subunit.interaction_sites
        rx, ry, rz = site.position[1] - subunit.position[1], site.position[2] - subunit.position[2], site.position[3] - subunit.position[3]

        site.position[1] = subunit.position[1] + Rxx * rx + Rxy * ry + Rxz * rz
        site.position[2] = subunit.position[2] + Ryx * rx + Ryy * ry + Ryz * rz
        site.position[3] = subunit.position[3] + Rzx * rx + Rzy * ry + Rzz * rz
    end

    for b in subunit.body_axes
        b[1], b[2], b[3] = Rxx * b[1] + Rxy * b[2] + Rxz * b[3], Ryx * b[1] + Ryy * b[2] + Ryz * b[3], Rzx * b[1] + Rzy * b[2] + Rzz * b[3]
    end
end
rotate!(subunit::RigidSubunit, axis::V, θ::Float64) where V <: AbstractVector = rotate!(subunit, axis[1], axis[2], axis[3], θ)

function rotate!(subunit::RigidSubunit, origin_x::Float64, origin_y::Float64, origin_z::Float64, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)
    sin, cos = sincos(θ)
    Rxx, Rxy, Rxz = cos + axis_x^2 * (1 - cos), axis_x * axis_y * (1 - cos) - axis_z * sin, axis_x * axis_z * (1 - cos) + axis_y * sin
    Ryx, Ryy, Ryz = axis_y * axis_x * (1 - cos) + axis_z * sin, cos + axis_y^2 * (1 - cos), axis_y * axis_z * (1 - cos) - axis_x * sin
    Rzx, Rzy, Rzz = axis_z * axis_x * (1 - cos) - axis_y * sin, axis_z * axis_y * (1 - cos) + axis_x * sin, cos + axis_z^2 * (1 - cos)

    rx, ry, rz = subunit.position[1] - origin_x, subunit.position[2] - origin_y, subunit.position[3] - origin_z

    subunit.position[1] = origin_x + Rxx * rx + Rxy * ry + Rxz * rz
    subunit.position[2] = origin_y + Ryx * rx + Ryy * ry + Ryz * rz
    subunit.position[3] = origin_z + Rzx * rx + Rzy * ry + Rzz * rz

    for site in subunit.interaction_sites
        rx, ry, rz = site.position[1] - origin_x, site.position[2] - origin_y, site.position[3] - origin_z

        site.position[1] = origin_x + Rxx * rx + Rxy * ry + Rxz * rz
        site.position[2] = origin_y + Ryx * rx + Ryy * ry + Ryz * rz
        site.position[3] = origin_z + Rzx * rx + Rzy * ry + Rzz * rz
    end

    for b in subunit.body_axes
        b[1], b[2], b[3] = Rxx * b[1] + Rxy * b[2] + Rxz * b[3], Ryx * b[1] + Ryy * b[2] + Ryz * b[3], Rzx * b[1] + Rzy * b[2] + Rzz * b[3]
    end
end
rotate!(subunit::RigidSubunit, origin::Vo, axis::Va, θ::Float64) where {Vo <: AbstractVector, Va <: AbstractVector} = rotate!(subunit, origin[1], origin[2], origin[3], axis[1], axis[2], axis[3], θ)
rotate!(subunit::RigidSubunit, origin::V, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64) where V <: AbstractVector = rotate!(subunit, origin[1], origin[2], origin[3], axis_x, axis_y, axis_z, θ)
rotate!(subunit::RigidSubunit, origin_x::Float64, origin_y::Float64, origin_z::Float64, axis::V, θ::Float64) where V <: AbstractVector = rotate!(subunit, origin_x, origin_y, origin_z, axis[1], axis[2], axis[3], θ)

function rotate!(subunits::Vector{RigidSubunit}, origin_x::Float64, origin_y::Float64, origin_z::Float64, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)
    for subunit in subunits
        rotate!(subunit, origin_x, origin_y, origin_z, axis_x, axis_y, axis_z, θ)
    end
end
rotate!(subunits::Vector{RigidSubunit}, origin::Vo, axis::Va, θ::Float64) where {Vo <: AbstractVector, Va <: AbstractVector} = rotate!(subunits, origin[1], origin[2], origin[3], axis[1], axis[2], axis[3], θ)
rotate!(subunits::Vector{RigidSubunit}, origin::V, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64) where V <: AbstractVector = rotate!(subunits, origin[1], origin[2], origin[3], axis_x, axis_y, axis_z, θ)
rotate!(subunits::Vector{RigidSubunit}, origin_x::Float64, origin_y::Float64, origin_z::Float64, axis::V, θ::Float64) where V <: AbstractVector = rotate!(subunits, origin_x, origin_y, origin_z, axis[1], axis[2], axis[3], θ)

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

function translate!(subunits::Vector{RigidSubunit}, Δx::Float64, Δy::Float64, Δz::Float64)
    for subunit in subunits
        translate!(subunit, Δx, Δy, Δz)
    end
end
translate!(subunits::Vector{RigidSubunit}, Δr::V) where V <: AbstractVector = translate!(subunits, Δr[1], Δr[2], Δr[3])

function morph!(subunit::RigidSubunit, site_displacements::Vector{Vector{Float64}})
    n, d = subunit.body_axes
    dxn = cross(d, n)
    for (i, a) in enumerate(site_displacements)
        subunit.interaction_sites[i].position .+= a[1] .* dxn .+ a[2] .* d .+ a[3] .* n
    end
end