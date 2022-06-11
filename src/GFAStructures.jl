module GFAStructures
    using StaticArrays

    """
    Flush output so that jobs can be monitored on a cluster.
    """
    @inline println(args...) = println(stdout, args...)
    @inline function println(io::IO, args...)
        Base.println(io, args...)
        flush(io)
    end

    include("subunit.jl")
    include("integrator.jl")
    include("system.jl")
end