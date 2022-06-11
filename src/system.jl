export System

mutable struct System
    subunits::Vector{Subunit}
    bonds::Vector{Bond}

    integrator::TemperatureQuench
end

function initialize_lattice(unit_cell::Vector{Subunit}, lattice_vectors::Vector{Vector{Float64}}, lattice_dimensions::Vector{Int64})
    subunits = Vector{Subunit}()

    a1, a2, a3 = lattice_vectors
    N1, N2, N3 = lattice_dimensions
    for i = 0 : N1 - 1, j = 0 : N2 - 1, k = 0 : N3 - 1
        for subunit in deepcopy(unit_cell)
            translate!(subunit, i * a1 .+ j * a2 .+ k * a3)
            push!(subunits, subunit)
        end
    end

    return subunits
end

function link(subunits::Vector{Subunit}, interactions::Vector{Tuple{Int64, Int64, I}}, neighbor_cutoff::Float64, bond_cutoff::Float64) where I <: Interaction
    bonds = Vector{Bond}()

    return bonds
end

function run_simulation!(system::System; num_steps::Int64 = 1, message_interval::Float64 = 10.0)

end

function format_for_mathematica!(system::System; params = [], file::String = "TEMP.txt")
    if !isdir(dirname(file))
        mkpath(dirname(file))
    end

    
end