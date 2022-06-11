export Subunit

mutable struct Subunit
    position::SVector{3, Float64}

    binding_sites::Vector{BindingSite}
end

function rotate!(subunit::Subunit, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)

end
rotate!(subunit::Subunit, axis::Vector{Float64}, θ::Float64) = rotate!(subunit, axis[1], axis[2], axis[3], θ)

function translate!(subunit::Subunit, Δx::Float64, Δy::Float64, Δz::Float64)

end
translate!(subunit::Subunit, r::Vector{Float64}) = translate!(subunit, r[1], r[2], r[3])