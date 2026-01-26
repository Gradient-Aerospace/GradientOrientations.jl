############
# AA to... #
############

"Returns the DirectionCosineMatrix for the given AxisAngle."
function aa2dcm(aa::AA{T}) where {T}
    s, c = sincos(aa.angle)
    r = aa.axis
    R = diagm(SVector{3, T}(c, c, c)) + (one(T) - c) .* (r * r') - s .* crs3(r)
    return DCM{T}(R)
end
Base.convert(::Type{DCM}, aa::AA) = aa2dcm(erp) # Type is not specified on the LHS.
Base.convert(::Type{DCM{T}}, aa::AA{T}) where {T} = aa2dcm(erp) # Type is not specified on the LHS.

"Returns the EulerRodriguesParameters for the given AxisAngle."
function aa2erp(aa::AA{T}) where {T}
    r = normalize(aa.axis)
    s, c = sincos(aa.angle/2)
    return ERP{T}(s * r[1], s * r[2], s * r[3], c)
end
Base.convert(::Type{ERP}, aa::AA) = aa2erp(aa) # Type is not specified on the LHS.
Base.convert(::Type{ERP{T}}, aa::AA{T}) where {T} = aa2erp(aa) # Type is not specified on the LHS.

"Returns the RotationVector for the given AxisAngle."
function aa2rv(aa::AA{T}) where {T}
    return RV{T}(aa.angle .* aa.axis)
end
Base.convert(::Type{RV}, aa::AA) = aa2rv(aa) # Type is not specified on the LHS.
Base.convert(::Type{RV{T}}, aa::AA{T}) where {T} = aa2rv(aa) # Type is not specified on the LHS.

"Returns the RollPitchYaw for the given AxisAngle."
function aa2rpy(aa::AA)
    return erp2rpy(aa2erp(aa)) # TODO: Be more direct.
end
Base.convert(::Type{RPY}, aa::AA) = aa2rpy(aa) # Type is not specified on the LHS.
Base.convert(::Type{RPY{T}}, aa::AA{T}) where {T} = aa2rpy(aa) # Type is not specified on the LHS.

#############
# DCM to... #
#############

"Returns the AxisAngle for the given DirectionCosineMatrix."
function dcm2aa(dcm::DCM)
    return erp2aa(dcm2erp(dcm)) # TODO: Be more direct.
end
Base.convert(::Type{AA}, dcm::DCM) = dcm2aa(dcm) # Type is not specified on the LHS.
Base.convert(::Type{AA{T}}, dcm::DCM{T}) where {T} = dcm2aa(dcm) # Type is not specified on the LHS.

"Returns the EulerRodriguesParameters for the given DirectionCosineMatrix."
function dcm2erp(dcm::DCM{T}) where {T}

    R = dcm.matrix

    # Split the conversion so as to divide by the largest possible number.
    if R[1,1] + R[2,2] + R[3,3] >= 0
        η4 = 0.5 * √(1. + R[1,1] + R[2,2] + R[3,3])
        α  = 0.25 / η4 # Necessarily safe from above
        return EulerRodriguesParameters{T}(
            α * (R[2,3] - R[3,2]),
            α * (R[3,1] - R[1,3]),
            α * (R[1,2] - R[2,1]),
            η4,
        )
    elseif R[1,1] - R[2,2] - R[3,3] >= 0
        η1 = 0.5 * √(1. + R[1,1] - R[2,2] - R[3,3])
        α  = 0.25 / η1 # Necessarily safe from above
        return EulerRodriguesParameters{T}(
            η1,
            α * (R[1,2] + R[2,1]),
            α * (R[3,1] + R[1,3]),
            α * (R[2,3] - R[3,2]),
        )
    elseif - R[1,1] + R[2,2] - R[3,3] >= 0
        η2 = 0.5 * √(1. - R[1,1] + R[2,2] - R[3,3])
        α  = 0.25 / η2 # Necessarily safe from above
        return EulerRodriguesParameters{T}(
            α * (R[1,2] + R[2,1]),
            η2,
            α * (R[3,2] + R[2,3]),
            α * (R[3,1] - R[1,3]),
        )
    else
        η3 = 0.5 * √(1. - R[1,1] - R[2,2] + R[3,3])
        α  = 0.25 / η3 # Safe if R is a DCM
        return EulerRodriguesParameters{T}(
            α * (R[1,3] + R[3,1]),
            α * (R[3,2] + R[2,3]),
            η3,
            α * (R[1,2] - R[2,1]),
        )
    end

end
Base.convert(::Type{ERP}, dcm::DCM) = dcm2erp(dcm) # Type is not specified on the LHS.
Base.convert(::Type{ERP{T}}, dcm::DCM{T}) where {T} = dcm2erp(dcm) # Type is not specified on the LHS.

"Returns the RotationVector for the given DirectionCosineMatrix."
function dcm2rv(dcm::DCM)
    return erp2rv(dcm2erp(dcm)) # TODO: Be more direct.
end
Base.convert(::Type{RV}, dcm::DCM) = dcm2rv(dcm) # Type is not specified on the LHS.
Base.convert(::Type{RV{T}}, dcm::DCM{T}) where {T} = dcm2rv(dcm) # Type is not specified on the LHS.

"Returns the RollPitchYaw for the given DirectionCosineMatrix."
function dcm2rpy(dcm::DCM)
    return erp2rpy(dcm2erp(dcm)) # TODO: Be more direct.
end
Base.convert(::Type{RPY}, dcm::DCM) = dcm2rpy(dcm) # Type is not specified on the LHS.
Base.convert(::Type{RPY{T}}, dcm::DCM{T}) where {T} = dcm2rpy(dcm) # Type is not specified on the LHS.

#############
# ERP to... #
#############

"Returns the AxisAngle for the given EulerRodriguesParameters."
function erp2aa(erp::EulerRodriguesParameters{T}) where {T}
    θ = 2 * acos(clamp(erp.s, -one(T), one(T)))
    m = sqrt(erp.x^2 + erp.y^2 + erp.z^2)
    if iszero(m)
        r = SVector{3,T}(one(T), zero(T), zero(T))
    else
        r = SVector{3,T}(erp.x/m, erp.y/m, erp.z/m)
    end
    return AA{T}(r, θ)
end
Base.convert(::Type{AA}, erp::ERP) = erp2aa(erp) # Type is not specified on the LHS.
Base.convert(::Type{AA{T}}, erp::ERP{T}) where {T} = erp2aa(erp) # Types are the same.
# Base.convert(::Type{AA{T}}, erp::ERP) = AA{T}(erp2aa(erp)) # Different types? Weird.

"Returns the DirectionCosineMatrix for the given EulerRodriguesParameters."
function erp2dcm(e::EulerRodriguesParameters{T}) where {T}
    x2 = e.x^2
    y2 = e.y^2
    z2 = e.z^2
    s2 = e.s^2
    m = @SMatrix [
        (x2 - y2 - z2 + s2)         2 * (e.x * e.y + e.s * e.z) 2 * (e.x * e.z - e.s * e.y);
        2 * (e.y * e.x - e.s * e.z)        (-x2 + y2 - z2 + s2) 2 * (e.y * e.z + e.s * e.x);
        2 * (e.z * e.x + e.s * e.y) 2 * (e.z * e.y - e.s * e.x)        (-x2 - y2 + z2 + s2);
    ]
    return DCM{T}(m)
end
Base.convert(::Type{DCM}, erp::ERP) = erp2dcm(erp) # Type is not specified on the LHS.
Base.convert(::Type{DCM{T}}, erp::ERP{T}) where {T} = erp2dcm(erp) # Types are the same.

"Returns the RotationVector for the given EulerRodriguesParameters."
function erp2rv(erp::EulerRodriguesParameters{T}) where {T}
    aa = erp2aa(erp)
    return RV{T}(aa.angle .* aa.axis)
end
Base.convert(::Type{RV}, erp::ERP) = erp2rv(erp) # Type is not specified on the LHS.
Base.convert(::Type{RV{T}}, erp::ERP{T}) where {T} = erp2rv(erp) # Types are the same.

"Returns the RollPitchYaw for the given EulerRodriguesParameters."
function erp2rpy(erp::EulerRodriguesParameters{T}) where {T}

    w = erp.s
    x = erp.x
    y = erp.y
    z = erp.z

    # Roll (x-axis rotation)
    sinr_cosp = 2 * (w * x + y * z)
    cosr_cosp = 1 - 2 * (x^2 + y^2)
    roll = atan(sinr_cosp, cosr_cosp)

    # Pitch (y-axis rotation)
    sinp = 2 * (w * y - z * x)

    # Clamp value to [-1, 1] to handle potential numerical noise
    if abs(sinp) >= 1
        pitch = copysign(π / 2, sinp) # 90 degrees if out of range
    else
        pitch = asin(sinp)
    end

    # Yaw (z-axis rotation)
    siny_cosp = 2 * (w * z + x * y)
    cosy_cosp = 1 - 2 * (y^2 + z^2)
    yaw = atan(siny_cosp, cosy_cosp)

    return RPY{T}(roll, pitch, yaw)

end
Base.convert(::Type{RPY}, erp::ERP) = erp2rpy(erp) # Type is not specified on the LHS.
Base.convert(::Type{RPY{T}}, erp::ERP{T}) where {T} = erp2rpy(erp) # Types are the same.

############
# RV to... #
############

"Returns the AxisAngle for the given RotationVector."
function rv2aa(rv::RV{T}) where {T}
    θ = norm(rv.vector)
    if iszero(θ) # If no rotation...
        return zero(AA{T})
    end
    return AA{T}(rv.vector ./ θ, θ)
end
Base.convert(::Type{AA}, rv::RV) = rv2aa(rv) # Type is not specified on the LHS.
Base.convert(::Type{AA{T}}, rv::RV{T}) where {T} = rv2aa(rv) # Type is not specified on the LHS.

"Returns the DirectionCosineMatrix for the given RotationVector."
function rv2dcm(rv::RV)
    return erp2dcm(rv2erp(rv)) # TODO: Be more direct.
end
Base.convert(::Type{DCM}, rv::RV) = rv2dcm(rv) # Type is not specified on the LHS.
Base.convert(::Type{DCM{T}}, rv::RV{T}) where {T} = rv2dcm(rv) # Type is not specified on the LHS.

"Returns the EulerRodriguesParameters for the given RotationVector."
function rv2erp(rv::RV)
    return aa2erp(rv2aa(rv)) # This is already efficient.
end
Base.convert(::Type{ERP}, rv::RV) = rv2erp(rv) # Type is not specified on the LHS.
Base.convert(::Type{ERP{T}}, rv::RV{T}) where {T} = rv2erp(rv) # Type is not specified on the LHS.

"Returns the RollPitchYaw for the given RotationVector."
function rv2rpy(rv::RV)
    return erp2rpy(rv2erp(rv))
end
Base.convert(::Type{RPY}, rv::RV) = rv2rpy(rv) # Type is not specified on the LHS.
Base.convert(::Type{RPY{T}}, rv::RV{T}) where {T} = rv2rpy(rv) # Type is not specified on the LHS.

#############
# RPY to... #
#############

"Returns the AxisAngle for the given RollPitchYaw."
function rpy2aa(rpy::RPY)
    return erp2aa(rpy2erp(rpy)) # TODO: Be more direct.
end
Base.convert(::Type{AA}, rpy::RPY) = rpy2aa(rpy) # Type is not specified on the LHS.
Base.convert(::Type{AA{T}}, rpy::RPY{T}) where {T} = rpy2aa(rpy) # Type is not specified on the LHS.

"Returns the DirectionCosineMatrix for the given RollPitchYaw."
function rpy2dcm(rpy::RPY)
    return erp2dcm(rpy2erp(rpy)) # TODO: Be more direct.
end
Base.convert(::Type{DCM}, rpy::RPY) = rpy2dcm(rpy) # Type is not specified on the LHS.
Base.convert(::Type{DCM{T}}, rpy::RPY{T}) where {T} = rpy2dcm(rpy) # Type is not specified on the LHS.

"Returns the EulerRodriguesParameters for the given RollPitchYaw."
function rpy2erp(rpy::RPY)
    return compose( # TODO: Be more direct.
        erpx(rpy.roll),
        erpy(rpy.pitch),
        erpz(rpy.yaw),
    )
end
Base.convert(::Type{ERP}, rpy::RPY) = rpy2erp(rpy) # Type is not specified on the LHS.
Base.convert(::Type{ERP{T}}, rpy::RPY{T}) where {T} = rpy2erp(rpy) # Type is not specified on the LHS.

"Returns the RotationVector for the given RollPitchYaw."
function rpy2rv(rpy::RPY)
    return erp2rv(rpy2erp(rpy)) # TODO: Be more direct.
end
Base.convert(::Type{RV}, rpy::RPY) = rpy2rv(rpy) # Type is not specified on the LHS.
Base.convert(::Type{RV{T}}, rpy::RPY{T}) where {T} = rpy2rv(rpy) # Type is not specified on the LHS.
