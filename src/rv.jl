export RotationVector, RV, RV_F64

"""
Represents an orientation using a rotation vector -- essentially the angle of rotation times
the axis of rotation between the two frames. It has a single field, `vector`.
"""
@kwdef struct RotationVector{T} <: AbstractOrientation{T}
    vector::SVector{3, T}
end
const RV = RotationVector
const RV_F64 = RotationVector{Float64}

################
# Constructors #
################

Base.convert(type::Type{<:RV}, v::AbstractVector) = type(v)

Base.zero(rv::RV{T}) where {T} = zero(typeof(rv))
Base.zero(::Type{RV}) = zero(RV_F64)
Base.zero(::Type{RV{T}}) where {T} = RV(zero(SVector{3, T}))

"Constructs a random RotationVector following a uniform distribution on SO(3)."
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{RV{T}}) where {T}
    return erp2rv(rand(rng, ERP{T}))
end

# TODO: We could do the deg2rad/rad2deg thing here too.

##############
# Operations #
##############

Base.inv(rv::RV) = RV(-rv.vector)
compose(a::RV, b::RV) = erp2rv(compose(rv2erp(a), rv2erp(b)))
reframe(rv::RV, v) = reframe(rv2aa(rv), v)
difference(a::RV, b::RV) = erp2rv(difference(rv2erp(a), rv2erp(b)))
distance(rv::RV) = norm(rv.vector)
interpolate(a::RV, b::RV, f) = erp2rv(interpolate(rv2erp(a), rv2erp(b), f))

#############
# Iteration #
#############

# TODO?

#################
# Miscellaneous #
#################

# The default `show` method is fine.

########################################
# Conversions to Non-Orientation Types #
########################################

# Get a Vector of the internal values.
function Base.convert(type::Type{<:AbstractVector}, rv::RotationVector)
    return convert(type, rv.vector)
end
