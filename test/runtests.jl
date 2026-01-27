using Test
using GradientOrientations
using GradientOrientations: AbstractOrientation
using LinearAlgebra
using StaticArrays

import GradientOrientations: crs3

# By default, this is pretty tight.
Base.isapprox(a::ERP, b::ERP; atol = eps(4π)) = distance(a, b) <= atol

# Obviously, by relying on ERP conversion to tell all orientations, we should do a good job
# of making sure the ERP stuff is right on its own (separate from these tests).
function Base.isapprox(a::AbstractOrientation, b::AbstractOrientation; kwargs...)
    return isapprox(
        convert(ERP, a),
        convert(ERP, b);
        kwargs...
    )
end

include("test_erp.jl")

# Now that we've tested ERPs, let's test the other types by converting to/from ERPs.

@testset "AA" begin

    a = AxisAngle(normalize(SA[1., 2., 3.]), 1.1)
    b = AxisAngle(normalize(SA[1., 0., -1.]), 2.2)
    v = SA[5., 2., -4.]
    f = 0.1

    @test aa2erp(zero(AA_F64)) ≈ zero(ERP_F64)

    @test reframe(a, v) ≈ reframe(aa2erp(a), v)
    @test a ⊗ b ≈ erp2aa(aa2erp(a) ⊗ aa2erp(b))
    @test difference(a, b) ≈ erp2aa(difference(aa2erp(a), aa2erp(b)))
    @test distance(a, b) ≈ distance(aa2erp(a), aa2erp(b))
    @test interpolate(a, b, f) ≈ erp2aa(interpolate(aa2erp(a), aa2erp(b), f))
    @test inv(a) ≈ erp2aa(inv(aa2erp(a)))

    rand(AA_F64) # Just test that we can do it.

    # Test conversions to all other types.
    tol = 1e-7
    @test aa2dcm(a) ≈ a atol = tol
    @test aa2erp(a) ≈ a atol = tol
    @test aa2rpy(a) ≈ a atol = tol
    @test aa2rv(a) ≈ a atol = tol

end

@testset "ERP conversions" begin

    a = normalize(ERP(1., 2., 3., 4.))

    tol = 1e-7
    @test erp2aa(a) ≈ a atol = tol
    @test erp2dcm(a) ≈ a atol = tol
    @test erp2rpy(a) ≈ a atol = tol
    @test erp2rv(a) ≈ a atol = tol

end

@testset "RV" begin

    a = RV(SA[1., 2., 3.])
    b = RV(SA[1., 0., -1.])
    v = SA[5., 2., -4.]
    f = 0.1

    @test rv2erp(zero(RV_F64)) ≈ zero(ERP_F64)

    @test reframe(a, v) ≈ reframe(rv2erp(a), v)
    @test a ⊗ b ≈ erp2rv(rv2erp(a) ⊗ rv2erp(b))
    @test difference(a, b) ≈ erp2rv(difference(rv2erp(a), rv2erp(b)))
    @test distance(a, b) ≈ distance(rv2erp(a), rv2erp(b))
    @test interpolate(a, b, f) ≈ erp2rv(interpolate(rv2erp(a), rv2erp(b), f))
    @test inv(a) ≈ erp2rv(inv(rv2erp(a)))

    rand(RV_F64) # Just test that we can do it.

    # Test conversions to all other types.
    tol = 1e-7
    @test rv2aa(a) ≈ a atol = tol
    @test rv2dcm(a) ≈ a atol = tol
    @test rv2erp(a) ≈ a atol = tol
    @test rv2rpy(a) ≈ a atol = tol

end

@testset "DCM" begin

    a = rv2dcm(RV(SA[1., 2., 3.]))
    b = rv2dcm(RV(SA[1., 0., -1.]))
    v = SA[5., 2., -4.]
    f = 0.1

    @test dcm2erp(zero(DCM_F64)) ≈ zero(ERP_F64)

    # Triggy conversions have much worse tolerances than the non-trig types.
    dcm_tol = 1e-7

    @test reframe(a, v) ≈ reframe(dcm2erp(a), v)
    @test a ⊗ b ≈ erp2dcm(dcm2erp(a) ⊗ dcm2erp(b)) atol = dcm_tol
    @test difference(a, b) ≈ erp2dcm(difference(dcm2erp(a), dcm2erp(b)))
    @test distance(a, b) ≈ distance(dcm2erp(a), dcm2erp(b))
    @test interpolate(a, b, f) ≈ erp2dcm(interpolate(dcm2erp(a), dcm2erp(b), f)) atol = dcm_tol
    @test inv(a) ≈ erp2dcm(inv(dcm2erp(a))) atol = dcm_tol

    rand(DCM_F64) # Just test that we can do it.

    # Test conversions to all other types.
    tol = 1e-7
    @test dcm2aa(a) ≈ a atol = tol
    @test dcm2erp(a) ≈ a atol = tol
    @test dcm2rpy(a) ≈ a atol = tol
    @test dcm2rv(a) ≈ a atol = tol

end

@testset "RPY" begin

    a = rv2rpy(RV(SA[1., 2., 3.]))
    b = rv2rpy(RV(SA[1., 0., -1.]))
    v = SA[5., 2., -4.]
    f = 0.1

    @test rv2erp(zero(RV_F64)) ≈ zero(ERP_F64)

    rpy_tol = 1e-7

    # Note that *all* of these fall back to the ERP implementation anyway, so these aren't
    # tests so much as truisms.
    @test reframe(a, v) ≈ reframe(rpy2erp(a), v)
    @test a ⊗ b ≈ erp2rpy(rpy2erp(a) ⊗ rpy2erp(b)) atol = rpy_tol
    @test difference(a, b) ≈ erp2rpy(difference(rpy2erp(a), rpy2erp(b))) atol = rpy_tol
    @test distance(a, b) ≈ distance(rpy2erp(a), rpy2erp(b))
    @test interpolate(a, b, f) ≈ erp2rpy(interpolate(rpy2erp(a), rpy2erp(b), f)) atol = rpy_tol
    @test inv(a) ≈ erp2rpy(inv(rpy2erp(a))) atol = rpy_tol

    rand(RPY_F64) # Just test that we can do it.

    # Test conversions to all other types.
    tol = 1e-7
    @test rpy2aa(a) ≈ a atol = tol
    @test rpy2dcm(a) ≈ a atol = tol
    @test rpy2erp(a) ≈ a atol = tol
    @test rpy2rv(a) ≈ a atol = tol

end

@testset "lots of orientations" begin

    # For binary operations, we'll need some reference to use.
    ref_erp = normalize(ERP(1., 2., 3., 4.))
    ref_aa = erp2aa(ref_erp)
    ref_dcm = erp2dcm(ref_erp)
    ref_rpy = erp2rpy(ref_erp)
    ref_rv = erp2rv(ref_erp)

    v = SA[2., -3., 4.]

    v_tol = 1e-6
    erp_tol = 1e-6

    # Let's construct a bunch of different orientations and see that operations using them
    # Doing this with roll, pitch, and yaw let's us test a lot of annoying conditions
    # directly.
    # for yaw in LinRange(-2π, 2π, 21)
    #     for pitch in LinRange(-π/2, π/2, 11)
    #         for roll in LinRange(-π, π, 21)
    for yaw in LinRange(-2π, 2π, 21)
        for pitch in LinRange(-π/2, π/2, 11)
            for roll in LinRange(-π, π, 21)

                rpy = RPY(roll, pitch, yaw)
                erp = rpy2erp(rpy)
                aa = erp2aa(erp)
                dcm = erp2dcm(erp)
                rv = erp2rv(erp)

                v_expected = reframe(erp, v)
                @test reframe(aa, v) ≈ v_expected atol = v_tol
                @test reframe(dcm, v) ≈ v_expected atol = v_tol
                @test reframe(rpy, v) ≈ v_expected atol = v_tol
                @test reframe(rv, v) ≈ v_expected atol = v_tol

                erp_expected = compose(erp, ref_erp)
                @test compose(aa, ref_aa) ≈ erp_expected atol = erp_tol
                @test compose(dcm, ref_dcm) ≈ erp_expected atol = erp_tol
                @test compose(rpy, ref_rpy) ≈ erp_expected atol = erp_tol
                @test compose(rv, ref_rv) ≈ erp_expected atol = erp_tol

                erp_expected = difference(erp, ref_erp)
                @test difference(aa, ref_aa) ≈ erp_expected atol = erp_tol
                @test difference(dcm, ref_dcm) ≈ erp_expected atol = erp_tol
                @test difference(rpy, ref_rpy) ≈ erp_expected atol = erp_tol
                @test difference(rv, ref_rv) ≈ erp_expected atol = erp_tol

                angle_expected = distance(erp, ref_erp)
                @test distance(aa, ref_aa) ≈ angle_expected atol = erp_tol
                @test distance(dcm, ref_dcm) ≈ angle_expected atol = erp_tol
                @test distance(rpy, ref_rpy) ≈ angle_expected atol = erp_tol
                @test distance(rv, ref_rv) ≈ angle_expected atol = erp_tol

                for f in LinRange(0., 1., 5)
                    erp_expected = interpolate(erp, ref_erp, f)
                    @test interpolate(aa, ref_aa, f) ≈ erp_expected atol = erp_tol
                    @test interpolate(dcm, ref_dcm, f) ≈ erp_expected atol = erp_tol
                    @test interpolate(rpy, ref_rpy, f) ≈ erp_expected atol = erp_tol
                    if !isapprox(interpolate(rpy, ref_rpy, f), erp_expected; atol = erp_tol)
                        @show rpy
                        @show ref_rpy
                        @show erp
                        @show ref_erp
                        @show f
                    end
                    @test interpolate(rv, ref_rv, f) ≈ erp_expected atol = erp_tol
                end

                erp_expected = inv(erp)
                @test inv(aa) ≈ erp_expected atol = erp_tol
                @test inv(dcm) ≈ erp_expected atol = erp_tol
                @test inv(rpy) ≈ erp_expected atol = erp_tol
                @test inv(rv) ≈ erp_expected atol = erp_tol

                eltype(aa) == Float64
                eltype(dcm) == Float64
                eltype(erp) == Float64
                eltype(rpy) == Float64
                eltype(rv) == Float64

            end
        end
    end

end
