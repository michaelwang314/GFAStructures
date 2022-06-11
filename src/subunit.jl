export Subunit
export rotate!
export translate!

mutable struct Subunit
    position::MVector{3, Float64}
    binding_sites::Vector{BindingSite}

    body_axis::NTuple{2, MVector{3, Float64}}
end

function rotate!(subunit::Subunit, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)
    sin, cos = sincos(θ)
    Rxx, Rxy, Rxz = cos + axis_x^2 * (1 - cos), axis_x * axis_y * (1 - cos) - axis_z * sin, axis_x * axis_z * (1 - cos) + axis_y * sin
    Ryy, Ryz = cos + axis_y^2 * (1 - cos), axis_y * axis_z * (1 - cos) - axis_x * sin
    Rzz = cos + axis_z^2 * (1 - cos)

    for site in subunit.binding_sites
        pos = site.position
        pos[1], pos[2], pos[3] = Rxx * pos[1] + Rxy * pos[2] + Rxz * pos[3], Rxy * pos[1] + Ryy * pos[2] + Ryz * pos[3], Rxz * pos[1] + Ryz * pos[2] + Rzz * pos[3]
    end

    b1, b2 = subunit.body_axis
    b1[1], b1[2], b1[3] = Rxx * b1[1] + Rxy * b1[2] + Rxz * b1[3], Rxy * b1[1] + Ryy * b1[2] + Ryz * b1[3], Rxz * b1[1] + Ryz * b1[2] + Rzz * b1[3]
    b2[1], b2[2], b2[3] = Rxx * b2[1] + Rxy * b2[2] + Rxz * b2[3], Rxy * b2[1] + Ryy * b2[2] + Ryz * b2[3], Rxz * b2[1] + Ryz * b2[2] + Rzz * b2[3]
end
rotate!(subunit::Subunit, axis::V, θ::Float64) where V <: AbstractVector = rotate!(subunit, axis[1], axis[2], axis[3], θ)

function translate!(subunit::Subunit, Δx::Float64, Δy::Float64, Δz::Float64)
    pos = subunit.position
    pos[1] += Δx
    pos[2] += Δy
    pos[3] += Δz

    for site in subunit.binding_sites
        pos = site.position
        pos[1] += Δx
        pos[2] += Δy
        pos[3] += Δz
    end
end
translate!(subunit::Subunit, Δr::V) where V <: AbstractVector = translate!(subunit, Δr[1], Δr[2], Δr[3])