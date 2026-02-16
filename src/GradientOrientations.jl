"A package for representing and operation on orientation-related types."
module GradientOrientations

using LinearAlgebra
using Random
using StaticArrays

"""
All subtypes of AbstractOrientation are expected to support:

* reframe
* compose
* difference
* distance
* interpolate
* other
* smallest

* Base.zero
* Base.inv
* Random.rand

To enable a type to support numerical integration, implement `rates` for the type.

It's helpful to write conversions for each abstract type, especially to and from
EulerRodriguesParameters. That provides a way to convert to all other types. Where more
direct conversions exist, those are helpful to implement as well.
"""
abstract type AbstractOrientation{T} end

# Necessary things
export reframe, compose, difference, distance, interpolate

# May not exist for all orientation types
export other, smallest, rate

# Users can `convert` from one type to another, but this lets them specify their intention
# more clearly and is easier to write.
export aa2dcm, aa2erp, aa2rpy, aa2rv
export dcm2aa, dcm2erp, dcm2rpy, dcm2rv
export erp2aa, erp2dcm, erp2rpy, erp2rv
export rv2aa, rv2dcm, rv2erp, rv2rpy
export rpy2aa, rpy2dcm, rpy2erp, rpy2rv

# Composition operator
export ⊗

include("utilities.jl")

# The primary types
include("erp.jl")
include("aa.jl")
include("dcm.jl")
include("rpy.jl")
include("rv.jl")
include("conversions.jl")

# Implement some fallback methods for AbstractOrientation that just convert to ERP, do the
# operation, and then convert back.
function reframe(a::AbstractOrientation, v)
    return reframe(convert(ERP, a), v)
end
function compose(a::T, b::T) where {T <: AbstractOrientation}
    return convert(T, compose(convert(ERP, a), convert(ERP, b)))
end
function difference(a::T, b::T) where {T <: AbstractOrientation}
    return convert(T, difference(convert(ERP, a), convert(ERP, b)))
end
function distance(a::T, b::T) where {T <: AbstractOrientation}
    return distance(convert(ERP, a), convert(ERP, b))
end
function interpolate(a::T, b::T, f) where {T <: AbstractOrientation}
    return convert(T, interpolate(convert(ERP, a), convert(ERP, b), f))
end
Base.eltype(::AbstractOrientation{T}) where {T} = T
Base.zero(x::AbstractOrientation) = zero(typeof(x))
Base.zero(type::Type{<:AbstractOrientation}) = zero(ERP{eltype(type)})
Base.inv(x::AbstractOrientation) = convert(typeof(x), inv(convert(ERP, x)))
Random.rand(x::AbstractOrientation) = convert(typeof(x), rand(typeof(convert(ERP, x))))

"Composition operator, with the same interface as `compose`."
⊗(a, b) = compose(a, b)

# Allow a user to compose multiple orientations.
compose(a, b, c, args...) = compose(compose(a, b), c, args...)

end # module GradientOrientations
