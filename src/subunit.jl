export BindingSite
export Subunit
export Bond

mutable struct BindingSite
    position::SVector{3, Float64}
end

mutable struct Subunit
    position::SVector{3, Float64}

    binding_sites::Vector{BindingSite}
end

mutable struct Bond{P <: Function, F <: Function}
    pair::Tuple{BindingSite, BindingSite}

    potential::P
    force::F
end

function compute_forces!(bonds::Vector{Bond})

end