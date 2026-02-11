export RollPitchYaw, RPY, RPY_F64, RPYDeg, RPYDeg_F64

"""
Represents an orientation as a roll, pitch, and yaw (rad) from a reference. Yaw is the first
rotation, followed by pitch, followed by roll. That is, this describes a frame oriented
`roll` around the x axis of a frame that's oriented `pitch` around the y axis of a frame
that's oriented `yaw` around the z axis of the reference frame.
"""
@kwdef struct RollPitchYaw{T} <: AbstractOrientation{T}
    roll::T
    pitch::T
    yaw::T
end
const RPY = RollPitchYaw
const RPY_F64 = RollPitchYaw{Float64}

"""
Represents an orientation as a roll, pitch, and yaw (deg) from a reference. Yaw is the first
rotation, followed by pitch, followed by roll. That is, this describes a frame oriented
`roll` around the x axis of a frame that's oriented `pitch` around the y axis of a frame
that's oriented `yaw` around the z axis of the reference frame.

This type can be converted to RollPitchYaw via `deg2rad`.
"""
@kwdef struct RollPitchYawDeg{T} <: AbstractOrientation{T}
    roll::T
    pitch::T
    yaw::T
end
const RPYDeg = RollPitchYawDeg
const RPYDeg_F64 = RollPitchYawDeg{Float64}

##############################
# Conversion To/From Degrees #
##############################

Base.deg2rad(rpy_deg::RPYDeg{T}) where {T} = RPY{T}(
    deg2rad(rpy_deg.roll),
    deg2rad(rpy_deg.pitch),
    deg2rad(rpy_deg.yaw),
)

Base.rad2deg(rpy::RPY{T}) where {T} = RPYDeg{T}(
    rad2deg(rpy.roll),
    rad2deg(rpy.pitch),
    rad2deg(rpy.yaw),
)

Base.convert(::Type{<:RPY}, rpy_deg::RPYDeg) = deg2rad(rpy_deg)
Base.convert(::Type{<:RPYDeg}, rpy::RPY) = rad2deg(rpy)

################
# Constructors #
################

Base.zero(rv::RPY{T}) where {T} = zero(typeof(rv))
Base.zero(::Type{RPY}) = zero(RPY_F64)
Base.zero(::Type{RPY{T}}) where {T} = RPY{T}(zero(T), zero(T), zero(T))

"Constructs a random RollPitchYaw following a uniform distribution on SO(3)."
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{RPY{T}}) where {T}
    return erp2rpy(rand(rng, ERP{T}))
end

##############
# Operations #
##############

# We use the abstract fallback for all operations. RPYs aren't really meant for operations.
