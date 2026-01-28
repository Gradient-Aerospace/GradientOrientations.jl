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
