# GradientOrientations.jl

*NOTE*: This is a sandbox for now. If we like this package, we should release it as open source under Orientations.jl.

This package is useful for specifying the orientation of a thing relative to another thing and for operating on those orientations.

The available orientation types are:

* `EulerRodriguesParameters` (aka `ERP`)
* `DirectionCosineMatrix` (aka `DCM`)
* `RotationVector`
* `AxisAngle`
* `ModifiedRodriguesParameters` (MRP)

Euler-Rodrigues Parameters are equivalent to "quaternions" using the JPL/Shuster convention common in aerospace. The scalar part is last. The composition rules are such that the order of operations works like direction cosine matrices. That is if `erp_BA` is the orientation of B wrt A and `erp_CB` is the orientation of C wrt B, then the orientation of C wrt A, `erp_CA`, is `compose(erp_CB, erp_BA)`, which is the same as `ERP(DCM(erp_CB) * DCM(erp_BA))`.

Just like $\otimes$ is often used to mean "composition" for Euler-Rodrigues symmetric parameters, this package uses the `compose` function two compose two orientations rather than `*`. Similarly, it uses `reframe(erp, v)` to express `v` in a frame orientation according to `erp`, rather than `erp * v`.

# Comparison With Other Packages

Rotations.jl is an excellent package that implements _active rotations_. When constructing its rotations, you say, "This _does something_ to a vector. `RotX(pi/4) * v` rotates `v` by `pi/4` radians about the x axis. This is exactly the opposite of Orientations, where `erp_BA = erpx(pi/4)` says, "frame B is rotated pi/4 radians from frame A", and if `v_A` is a vector expressed in frame A, then `v_B = reframe(erp_BA, v_A)` is that same vector expressed in frame B.

## Documentation

### Operations on Orientations

* [`compose`](@ref): "Composes" two orientations, with the same order and meaning as the multiplication of two direction cosine matrices.
* [`difference`](@ref): Returns orientation corresponding to the difference in two orientations.
* [`interpolate`](@ref): Interpolates between two orientations.
* [`inv`](@ref): Inverts an orientation (`compose(η, inv(η))` means "no rotation").
* [`normalize`](@ref): Divides the Euler-Rodrigues parameters by their 2-norm for a proper unit-norm orientation. (Euler-Rodrigues parameters are _not_ automatically normalized in other operations.)
* [`other`](@ref): Returns an equivalent orientation, going around the opposite rotation axis -- going the "other way around".
* [`rate`](@ref): Returns the time-derivative orientation given a rotation rate expressed in the "this" frame (as opposed to the "wrt that" frame).
* [`reframe`](@ref): Given the orientation of frame B wrt A and a vector expressed in A, this returns the vector expressed in B (see above). Given a Vector for the vector input, this will return a Vector. Given a Tuple, it will return a Tuple. Given any other type, `T`, it will attemp to construct it as `T(v1, v2, v3)` (works for `SVectors` from `StaticArrays` for instance).
* [`smallest`](@ref): Returns the representation with the smaller rotation angle from the reference -- the "shortest way around".

### Conversions

TODO: Should we keep these or rely on constructors?

* [`aa2ep`](@ref): Converts axis-angle representation to Euler-Rodrigues parameters
* [`ep2aa`](@ref): Converts Euler-Rodrigues parameters to axis-angle representation (the axis will be a tuple)
* [`rv2ep`](@ref): Converts a rotation vector to Euler-Rodrigues parameters
* [`dcm2ep`](@ref): Converts a direction cosine matrix to Euler-Rodrigues parameters

## Euler-Rodrigues Parameters

This is an implementation of Euler-Rodrigues symmetric parameters, as described by Shuster in [A Survey of Attitude Representations](http://malcolmdshuster.com/Pub_1993h_J_Repsurv_scan.pdf). These are equivalent to unit quaternions used for rotations, but whereas quaternion operations are more general, have many competing conventions, and typically represent "vector" rotations rather than "frame" rotations, Euler parameters are used only for representations of orientation (frame rotations) and can be implemented directly from Shuster's paper without confusion.

The conventions here are that the rotation describes how the frame is rotated, the scalar part is last, and EP "composition" has the same order of arguments as rotation matrices for the equivalent functionalitly, all consistent with Shuster.

Note that for speed and other reasons, EP are not automatically normalized after calculations. Feel free to normalize when that makes sense for your application.

## TODO: Documentation of other types.

## Example

Let's consider two frames, frame A and frame B, with B rotated 45 degrees from frame A about its y-axis.

```julia
r = [0., 1., 0.] # Rotation axis
θ = π/4 # Rotation angle
η_BA = aa2erp(r, θ) # Euler parameters for the orientation of B wrt A, from axis-angle notation
# output
4-element EP{Float64}:
 0.0
 0.3826834323650898
 0.0
 0.9238795325112867
```

If something were 2 units along the A frame's x axis, that same vector expressed in B would be `[√2, 0, √2]`.

```julia
v_A = [2., 0., 0.]
v_B = reframe(η_BA, v_A)
# output
3-element Vector{Float64}:
 1.414213562373095
 0.0
 1.4142135623730951
```

Now consider frame C rotated 45 degree from B about B's z axis.

```julia
η_CB = aa2erp([0., 0., 1.], π/4)
v_C = reframe(η_CB, v_B)
# output
3-element Vector{Float64}:
  0.9999999999999998
 -1.0
  1.4142135623730951
```

We can "compose" the two rotations to represent C wrt A, in the same order as rotation matrices:

```julia
η_CA = compose(η_CB, η_BA)
v_C = reframe(η_CA, v_A)
# output
3-element Vector{Float64}:
  0.9999999999999998
 -0.9999999999999999
  1.414213562373095
```

Naturally, `v_C` is the same either way.
