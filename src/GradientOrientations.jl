"""
TODO
"""
module GradientOrientations

using LinearAlgebra
using Random
using StaticArrays

abstract type AbstractOrientation end

export reframe, compose, difference, distance, interpolate, other, smallest, rate

export aa2erp, erp2aa, rv2erp, erp2rv, dcm2erp, erp2dcm

# TODO: What interface do orientation types support?
#
# reframe
# compose
# difference
# distance
# interpolate
# other
# smallest
#
# Base.one
# Base.zero
# Random.rand
# Base.inv
# `*` to mean "compose" or "reframe"?
#
# Or do we even want "orientation types"? ERP is potentially the only custom type we need.
#
# DCM: SMatrix
# RotationVector: SVector
# AxisAngle: Float64, SVector
# MRP: SVector
#
# On the other hand, there's no real harm in having types for these things.

# TODO: Include this simply because it's so useful for orientations/rotations? Where should
# we put functions related to orientations/rotations in general? This seems like a great
# package for them.
function crs3(r)
    @assert length(r) == 3 "The length of the input to `crs3` was $(length(r)) but should be 3."
    return @SMatrix [0 -r[3] r[2]; r[3] 0 -r[1]; -r[2] r[1] 0]
end

include("erp.jl")

end # module GradientOrientations
