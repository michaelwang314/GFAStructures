module GFAStructures
    using StaticArrays
    using LinearAlgebra
    using Serialization

    export RigidSubunit, InteractionSite, rotate!, translate!
    export NeighborList, CellList, PairList
    export Interaction, HarmonicBond, LennardJones, HertzianSphere
    export System, initialize_lattice, find_neighbors, run_simulation!, save!, load

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
    include("integrators.jl")
    include("system.jl")
end