
@testset "aa2erp, erp2aa, and compose" begin

    r = [1., 0., 0.]
    theta1 = 1.
    theta2 = 2.
    ep1 = aa2erp(r, theta1)
    erp2 = aa2erp(r, theta2)
    epf = compose(erp2, ep1)
    rf, thetaf = erp2aa(epf)
    @test norm(rf .- r) < 1e-10
    @test thetaf ≈ theta1 + theta2 atol=1e-8

end

@testset "rv2erp, erp2rv, dcm2erp, erp2dcm" begin

    r = [1/sqrt(2.), 1/sqrt(2.), 0.]
    theta = -2.
    rv = theta * r
    ep_from_rv = rv2erp(rv)
    ep_from_aa = aa2erp(r, theta)
    R = cos(theta) * Matrix(I, 3, 3) + (1 - cos(theta)) * (r * r') - sin(theta) * crs3(r)
    ep_from_R = dcm2erp(R)
    @test [ep_from_rv...] ≈ [ep_from_aa...] atol=1e-15
    @test [ep_from_rv...] ≈ [ep_from_R...] atol=1e-15

    rv2 = erp2rv(ep_from_aa)
    @test rv2 ≈ rv atol=1e-12

    R2 = erp2dcm(ep_from_rv)
    @test R2 ≈ R atol=1e-12

    @test rv2erp([0., 0., 0.]) == zero(ERP{Float64})

    @test all(dcm2erp([ 1 0 0; 0 -1 0; 0 0 -1]) .≈ ERP(1., 0., 0., 0.))
    @test all(dcm2erp([-1 0 0; 0  1 0; 0 0 -1]) .≈ ERP(0., 1., 0., 0.))
    @test all(dcm2erp([-1 0 0; 0 -1 0; 0 0  1]) .≈ ERP(0., 0., 1., 0.))

end

@testset "reframe with Vector" begin
    v_A = [1., 0., 0.]
    η_BA = rv2erp(π/4 * [0., 1., 0.])
    v_B = reframe(η_BA, v_A)
    @test typeof(v_B) == Vector{Float64}
    @test v_B ≈ [1/sqrt(2), 0, 1/sqrt(2)]
end

@testset "reframe with SVector" begin
    v_A = SA[0., 2., 0.]
    η_BA = rv2erp(π/4 * SA[0., 0., 1.])
    v_B = reframe(η_BA, v_A)
    @test typeof(v_B) == SVector{3, Float64}
    @test reframe(η_BA, v_A) ≈ SA[2/sqrt(2), 2/sqrt(2), 0]
end

@testset "norm, squared_norm, normalize, inv, other, smallest" begin
    ep = normalize(ERP(1., 2., 3., 4.))
    @test norm(ep) ≈ 1.
    epi = inv(ep)
    epf = normalize(compose(epi, ep))
    @test [epf...] ≈ [0., 0., 0., 1.]
    _, angle = erp2aa(other(epf))
    @test angle ≈ 2π
    r, angle = erp2aa(smallest(epf))
    @test angle ≈ 0.
    r, angle = erp2aa(smallest(other(epf)))
    @test angle ≈ 0.
end

@testset "rate" begin

    r      = [1., 0., 0.]
    theta  = 0.
    ep     = aa2erp(r, theta)
    w      = 1. * r
    ep_dot = rate(ep, w)

    # While we're here, test that this results in zero rotation.
    @test all(ep .≈ ERP(0., 0., 0., 1.))
    @test all(ep .≈ ERP([0., 0., 0., 1.]))
    @test all(ep .≈ zero(ERP{Float64}))

    # These are only true for rotations about a single axis.
    @test ep_dot.x ≈ cos(theta/2) * 1/2 * w[1]
    @test ep_dot.y ≈ cos(theta/2) * 1/2 * w[2]
    @test ep_dot.z ≈ cos(theta/2) * 1/2 * w[3]
    @test ep_dot.s ≈ -sin(theta/2) * 1/2 * (w ⋅ r)

    r      = [0., 1., 0.]
    theta  = 1.
    ep     = aa2erp(r, theta)
    w      = 2. * r
    ep_dot = rate(ep, w)

    # These are only true for rotations about a single axis.
    @test ep_dot.x ≈ cos(theta/2) * 1/2 * w[1]
    @test ep_dot.y ≈ cos(theta/2) * 1/2 * w[2]
    @test ep_dot.z ≈ cos(theta/2) * 1/2 * w[3]
    @test ep_dot.s ≈ -sin(theta/2) * 1/2 * (w ⋅ r)

    r      = [0., 0., 1.]
    theta  = 2π - 0.001
    ep     = aa2erp(r, theta)
    w      = 3. * r
    ep_dot = rate(ep, w)

    # These are only true for rotations about a single axis.
    @test ep_dot.x ≈ cos(theta/2) * 1/2 * w[1]
    @test ep_dot.y ≈ cos(theta/2) * 1/2 * w[2]
    @test ep_dot.z ≈ cos(theta/2) * 1/2 * w[3]
    @test ep_dot.s ≈ -sin(theta/2) * 1/2 * (w ⋅ r)

    @test all(ep_dot .≈ 1/2 * compose(ERP(w..., 0.), ep))

end

@testset "rate norm correction" begin

    # Define a function to get the magnitude of the error between two sets of EPs.
    angle_from_vector(ep) = 2 * asin(norm([ep.x, ep.y, ep.z]))
    eperr(ep1, erp2) = angle_from_vector(compose(ep1, inv(erp2)))

    # Let's set up a little sim where we have a set of EulerParameters. One will be the
    # truth values, propagated numerical integration (Euler method). Another will start at a
    # non-norm value and be integrated without norm correction. The last will start at the
    # same non-norm value and be integrated _with_ norm correction. We'll first rotate a bit
    # about x, then about y, and then about z, and we'll repeat that sequence for a bunch of
    # samples. We're assuming the rotation is constant for each sample. We'll do this
    # for a few different starting points for the true EPs.
    for ep0 in (
        ERP(0., 0., 0., 1.),
        ERP(0., 0., 0., -1.),
        normalize(ERP(1., 1., 1., 1.)),
        other(normalize(ERP(1., 1., 1., 1.))), # The long way around version of the above
        inv(normalize(ERP(1., 1., 1., 1.))),
        other(inv(normalize(ERP(1., 1., 1., 1.)))),
    )
        # The non-norm error we'll add on
        Δep = ERP(0., 0., 0., 0.000001)

        ep  = ep0       # Corrected propagation from the right initial condition
        epu = ep0 + Δep # Uncorrected propagation from non-norm i.c.
        epc = ep0 + Δep # Corrected propagation from non-norm i.c.
        Δt  = 0.1       # Time step (s)
        for z = 1:1000  # Repeat the x, y, z rotation a whole bunch.
            for k = 1:3
                w = [
                    k == 1 ? 0.1 : 0., # First rotate about x, then y, then z.
                    k == 2 ? 0.2 : 0.,
                    k == 3 ? 0.3 : 0.,
                ]
                ep  += rate(ep,  w)     * Δt
                epu += rate(epu, w, 0.) * Δt
                epc += rate(epc, w)     * Δt
            end
        end

        # Check that the norm is better when corrected than without correction.
        @test abs(norm(epc) - 1.) < abs(norm(epu) - 1.)

        # Check that the propagation error with the correction term is better than the error
        # without.
        @test eperr(epc, ep) <= eperr(epu, ep)

    end

end

@testset "operations" begin

    x = [1., 2., 3., 4.]
    y = [5., 6., 7., 8.]
    @test ERP(x + y) == ERP(x) + ERP(y)
    @test ERP(2x) == 2. * ERP(x)
    @test ERP(2x) == ERP(x) * 2.
    @test all((x...,) .== convert(NTuple{4,Float64}, ERP(x)))
    @test x == convert(Vector{Float64}, ERP(x))

    ep = ERP(x)
    @test ep[1] == x[1]
    @test ep[4] == x[4]
    @test all((ep...,) .== (x...,))

end

@testset "misc" begin

    ep = ERP(0., 0., 0., 1.)

    # "Test" the `show` method.
    @show ep

    @test length(ep) == 4

    [(el for el in ep)...,]

    x = [1., 2., 3., 4.]
    @test all((x...,) .== (convert(ERP, x)...,))
    @test all((x...,) .== (convert(ERP{Float64}, x)...,))
    @test all((x...,) .== (ERP(x)...,))

    dict = convert(Dict{String, Float64}, ERP(x))
    @test dict["x"] == 1.
    @test dict["y"] == 2.
    @test dict["z"] == 3.
    @test dict["s"] == 4.

end

@testset "erpx, erpy, erpz" begin

    # This just tests that I didn't have a typo, basically.
    @test all(erpx(1.) .== aa2erp([1., 0., 0.], 1.))
    @test all(erpy(2.) .== aa2erp([0., 1., 0.], 2.))
    @test all(erpz(3.) .== aa2erp([0., 0., 1.], 3.))

    # Test that we can construct with an Irrational.
    @test erp2aa(erpx(π))[2] ≈ π
    @test erp2aa(erpy(π))[2] ≈ π
    @test erp2aa(erpz(π))[2] ≈ π
    @test erp2aa(aa2erp([1/√(2), 1/√(2), 0.], π))[2] ≈ π

end

@testset "difference: zero difference" begin

    # This ERP is not normalized on purpose to make sure a non-normal ERP is properly handled.
    a = ERP(1., 2., 3., -4.)

    b1 = a
    b2 = other(a)
    b3 = normalize(a)
    b4 = normalize(other(a))

    @test distance(difference(a, b1), zero(ERP{Float64})) ≈ 0 atol = eps(1.)
    @test distance(difference(a, b2), zero(ERP{Float64})) ≈ 0 atol = eps(1.)
    @test distance(difference(a, b3), zero(ERP{Float64})) ≈ 0 atol = eps(1.)
    @test distance(difference(a, b4), zero(ERP{Float64})) ≈ 0 atol = eps(1.)

end

@testset "difference: non-zero difference" begin

    a = normalize(ERP(1., 2., 3., -4.))

    η_diff_true = aa2erp([1., 0., 0.], π/8)

    b1 = compose(η_diff_true, a)
    b2 = compose(η_diff_true, other(a))

    @test distance(difference(b1, a), η_diff_true) ≈ 0. atol = eps(1.)
    @test distance(difference(b2, a), η_diff_true) ≈ 0. atol = eps(1.)

end

@testset "distance" begin

    θ = deg2rad(15.)

    @test distance(erpx(θ), zero(ERP{Float64})) ≈ θ
    @test distance(erpy(θ), zero(ERP{Float64})) ≈ θ

    # Other way around
    @test distance(erpx(2π + θ), zero(ERP{Float64})) ≈ θ
    @test distance(other(erpy(θ)), zero(ERP{Float64})) ≈ θ

end

include("test_erp_interp.jl")
