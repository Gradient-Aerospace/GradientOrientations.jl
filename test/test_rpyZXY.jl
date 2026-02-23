@testset "RPYZXY" begin

    a = rv2rpyZXY(RV(SA[1., 2., 3.]))
    b = rv2rpyZXY(RV(SA[1., 0., -1.]))
    v = SA[5., 2., -4.]
    f = 0.1

    rpy_tol = 1e-7

    # Note that *all* of these fall back to the ERP implementation anyway, so these aren't
    # tests so much as truisms.
    @test reframe(a, v) ≈ reframe(rpyZXY2erp(a), v)
    @test a ⊗ b ≈ erp2rpyZXY(rpyZXY2erp(a) ⊗ rpyZXY2erp(b)) atol = rpy_tol
    @test difference(a, b) ≈ erp2rpyZXY(difference(rpyZXY2erp(a), rpyZXY2erp(b))) atol = rpy_tol
    @test distance(a, b) ≈ distance(rpyZXY2erp(a), rpyZXY2erp(b))
    @test interpolate(a, b, f) ≈ erp2rpyZXY(interpolate(rpyZXY2erp(a), rpyZXY2erp(b), f)) atol = rpy_tol
    @test inv(a) ≈ erp2rpyZXY(inv(rpyZXY2erp(a))) atol = rpy_tol

    rand(RPYZXY_F64) # Just test that we can do it.

    # Test conversions to all other types.
    tol = 1e-7
    @test rpyZXY2aa(a) ≈ a atol = tol
    @test rpyZXY2dcm(a) ≈ a atol = tol
    @test rpyZXY2erp(a) ≈ a atol = tol
    @test rpyZXY2rv(a) ≈ a atol = tol

    # Conversion to specific types.
    rpy = RPYZXY(0.1, 0.2, 0.3)
    dcm = rpyZXY2dcm(rpy)
    dcm_expected = Ry(rpy.pitch) * Rx(rpy.roll) * Rz(rpy.yaw)
    @test all(dcm.matrix .≈ dcm_expected.matrix)

    erp = rpyZXY2erp(rpy)
    erp_expected = compose(erpy(rpy.pitch), erpx(rpy.roll), erpz(rpy.yaw))
    @test all(smallest(erp) .≈ smallest(erp_expected))

    # Check that the inverse conversions are consistent.
    @test dcm2rpyZXY(dcm) ≈ rpy atol = tol
    @test erp2rpyZXY(erp) ≈ rpy atol = tol

    # Gimbal lock checks (roll = ±90°). The yaw/pitch split is not unique here, but the
    # represented orientation should still match.
    for roll in (π/2, -π/2)
        rpy_sing = RPYZXY(roll, 0.4, -0.5)
        dcm_sing = rpyZXY2dcm(rpy_sing)
        erp_sing = rpyZXY2erp(rpy_sing)
        @test dcm2rpyZXY(dcm_sing) ≈ rpy_sing atol = tol
        @test erp2rpyZXY(erp_sing) ≈ rpy_sing atol = tol
    end

    # Test conversion to/from degrees.
    deg = [10., 20., -30.]
    rpy_deg = convert(RPYZXYDeg, RPYZXY(deg2rad.(deg)...))
    @test rpy_deg.roll ≈ deg[1]
    @test rpy_deg.pitch ≈ deg[2]
    @test rpy_deg.yaw ≈ deg[3]
    rpy = convert(RPYZXY, rpy_deg)
    @test rpy.roll ≈ deg2rad(deg[1])
    @test rpy.pitch ≈ deg2rad(deg[2])
    @test rpy.yaw ≈ deg2rad(deg[3])

end
