using GradientOrientations
using Test
using StaticArrays
using LinearAlgebra

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

@testset "dcm2aa" begin

    # Test 1: Identity matrix should give zero angle
    R_identity = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    dcm_identity = DCM(R_identity)
    aa_identity = dcm2aa(dcm_identity)
    @test aa_identity.angle ≈ 0.0 atol=1e-10
    @test norm(aa_identity.axis) == 1.

    # Test 2: 90 degree rotation about x-axis
    aa_90x = AxisAngle(SA[1.0, 0.0, 0.0], π/2)
    dcm_90x = aa2dcm(aa_90x)
    aa_from_dcm = dcm2aa(dcm_90x)
    @test aa_from_dcm.angle ≈ π/2 atol=1e-10
    @test norm(aa_from_dcm.axis .- SA[1.0, 0.0, 0.0]) < 1e-10

    # Test 3: 90 degree rotation about y-axis
    aa_90y = AxisAngle(SA[0.0, 1.0, 0.0], π/2)
    dcm_90y = aa2dcm(aa_90y)
    aa_from_dcm_y = dcm2aa(dcm_90y)
    @test aa_from_dcm_y.angle ≈ π/2 atol=1e-10
    @test norm(aa_from_dcm_y.axis .- SA[0.0, 1.0, 0.0]) < 1e-10

    # Test 4: 90 degree rotation about z-axis
    aa_90z = AxisAngle(SA[0.0, 0.0, 1.0], π/2)
    dcm_90z = aa2dcm(aa_90z)
    aa_from_dcm_z = dcm2aa(dcm_90z)
    @test aa_from_dcm_z.angle ≈ π/2 atol=1e-10
    @test norm(aa_from_dcm_z.axis .- SA[0.0, 0.0, 1.0]) < 1e-10

    # Test 5: Arbitrary rotation (45 degrees about [1,1,1] axis)
    axis_arbitrary = normalize(SA[1.0, 1.0, 1.0])
    angle_arbitrary = π/4
    aa_arbitrary = AxisAngle(axis_arbitrary, angle_arbitrary)
    dcm_arbitrary = aa2dcm(aa_arbitrary)
    aa_from_dcm_arb = dcm2aa(dcm_arbitrary)
    @test aa_from_dcm_arb.angle ≈ angle_arbitrary atol=1e-10
    @test norm(aa_from_dcm_arb.axis .- axis_arbitrary) < 1e-10

    # Test 6: Large angle (180 degrees)
    axis_180 = normalize(SA[1.0, 0.0, 1.0])
    angle_180 = 1π
    aa_180 = AxisAngle(axis_180, angle_180)
    dcm_180 = aa2dcm(aa_180)
    aa_from_dcm_180 = dcm2aa(dcm_180)
    # This rotation could go either way.
    if aa_from_dcm_180.angle < 0
        @test aa_from_dcm_180.angle ≈ angle_180 atol=1e-10
        @test norm(aa_from_dcm_180.axis .- -axis_180) < 1e-10
    else
        @test aa_from_dcm_180.angle ≈ angle_180 atol=1e-10
        @test norm(aa_from_dcm_180.axis .- axis_180) < 1e-10
    end

    # Test 7: Small angle (very small rotation)
    axis_small = normalize(SA[1.0, 2.0, 3.0])
    angle_small = 1e-5
    aa_small = AxisAngle(axis_small, angle_small)
    dcm_small = aa2dcm(aa_small)
    aa_from_dcm_small = dcm2aa(dcm_small)
    @test aa_from_dcm_small.angle ≈ angle_small atol=1e-8
    @test norm(aa_from_dcm_small.axis .- axis_small) < 1e-6

    # Test 8: Round-trip conversion - AA -> DCM -> AA
    aa_original = AxisAngle(normalize(SA[1.0, 2.0, 3.0]), 0.7)
    dcm_intermediate = aa2dcm(aa_original)
    aa_final = dcm2aa(dcm_intermediate)
    @test aa_final.angle ≈ aa_original.angle atol=1e-10
    @test norm(aa_final.axis .- aa_original.axis) < 1e-10

    # Test 9: Type preservation (Float64)
    aa_typed = AxisAngle(normalize(SA[1.0, 1.0, 0.0]), 0.5)
    dcm_typed = aa2dcm(aa_typed)
    aa_from_dcm_typed = dcm2aa(dcm_typed)
    @test typeof(aa_from_dcm_typed) <: AxisAngle{Float64}
    @test typeof(aa_from_dcm_typed.axis) <: SVector{3, Float64}
    @test typeof(aa_from_dcm_typed.angle) <: Float64

end
