export RollPitchYawZXY, RPYZXY, RPYZXY_F64, RPYZXYDeg, RPYZXYDeg_F64

"""
Represents an orientation as a roll, pitch, and yaw (rad) from a reference, using the 3-1-2
Euler angle rotation order.

Yaw is the first rotation (about the z axis of the reference), followed by roll (about the x
axis of the intermediate frame), followed by pitch (about the y axis of the second
intermediate frame). Equivalently, the corresponding direction cosine matrix is:

    Ry(pitch) * Rx(roll) * Rz(yaw)

This convention can be useful for tailsitter / VTOL vehicles where `roll` approaches ±90° in
hover, avoiding the classic 3-2-1 pitch gimbal lock.
"""
@kwdef struct RollPitchYawZXY{T} <: AbstractOrientation{T}
    roll::T
    pitch::T
    yaw::T
end
const RPYZXY = RollPitchYawZXY
const RPYZXY_F64 = RollPitchYawZXY{Float64}

"""
Represents an orientation as a roll, pitch, and yaw (deg) from a reference, using the 3-1-2
Euler angle rotation order.

This type can be converted to `RollPitchYawZXY` via `deg2rad`.
"""
@kwdef struct RollPitchYawZXYDeg{T} <: AbstractOrientation{T}
    roll::T
    pitch::T
    yaw::T
end
const RPYZXYDeg = RollPitchYawZXYDeg
const RPYZXYDeg_F64 = RollPitchYawZXYDeg{Float64}

##############################
# Conversion To/From Degrees #
##############################

Base.deg2rad(rpy_deg::RPYZXYDeg{T}) where {T} = RPYZXY{T}(
    deg2rad(rpy_deg.roll),
    deg2rad(rpy_deg.pitch),
    deg2rad(rpy_deg.yaw),
)

Base.rad2deg(rpy::RPYZXY{T}) where {T} = RPYZXYDeg{T}(
    rad2deg(rpy.roll),
    rad2deg(rpy.pitch),
    rad2deg(rpy.yaw),
)

Base.convert(::Type{<:RPYZXY}, rpy_deg::RPYZXYDeg) = deg2rad(rpy_deg)
Base.convert(::Type{<:RPYZXYDeg}, rpy::RPYZXY) = rad2deg(rpy)

################
# Constructors #
################

Base.zero(rpy::RPYZXY{T}) where {T} = zero(typeof(rpy))
Base.zero(::Type{RPYZXY}) = zero(RPYZXY_F64)
Base.zero(::Type{RPYZXY{T}}) where {T} = RPYZXY{T}(zero(T), zero(T), zero(T))

##############
# Operations #
##############

# We use the abstract fallback for all operations. RPYZXYs aren't really meant for
# operations.
