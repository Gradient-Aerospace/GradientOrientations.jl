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
