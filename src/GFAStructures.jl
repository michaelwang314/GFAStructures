module GFAStructures
    using StaticArrays
    using LinearAlgebra
    using Serialization

    export RigidSubunit, InteractionSite, rotate!, translate!, morph!
    export NeighborList, CellList, FixedPairList
    export Interaction, HarmonicBond, LennardJones, HertzianSphere
    export ExternalForce, AnisotropicSpring, ConstantForce
    export GradientDescent
    export System, initialize_lattice, create_interaction_matrix, find_neighbors, sort_by_id, get_neighbor_statistics, run_simulation!, format_for_mathematica, save!, load

    """
    Flush output so that jobs can be monitored on a cluster.
    """
    @inline println(args...) = println(stdout, args...)
    @inline function println(io::IO, args...)
        Base.println(io, args...)
        flush(io)
    end

    include("subunit.jl")
    include("neighbor_lists.jl")
    include("interactions.jl")
    include("external_forces.jl")
    include("integrators.jl")
    include("system.jl")
end