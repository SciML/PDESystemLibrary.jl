using PDESystemLibrary
using Test

@testset "Precompiled analytic solutions" begin
    systems = get_pdesys_with_tags(String[])
    system = only(filter(system -> system.name == :inviscid_burgers_monotonic, systems))
    analytic = only(values(system.analytic_func))

    @test analytic([], 0.0, 0.5) == 0.5
end

@testset "convection_diffusion analytic solutions" begin
    systems = get_pdesys_with_tags(String[])
    # Values at t = 2 from an independent Crank-Nicolson finite-difference solve
    # (dz = L / 40000, dt = 2.5e-4).
    cases = (
        (:convdiff1, [0.01, 1.0], 1.0, [(0.05, 0.99326203)]),
        (:convdiff2, [1.0, -1.0], 5.0, [(2.5, 0.0451322), (4.5, 0.574723)]),
        (:convdiff3, [0.01, -3.0], 10.0, Tuple{Float64, Float64}[]),
        (:convdiff4, [0.5, 3.0], 10.0, [(0.5, 0.0083057), (5.0, 0.79956466)]),
    )
    for (name, ps, L, fdvalues) in cases
        system = only(filter(system -> system.name == name, systems))
        analytic = only(values(system.analytic_func))
        for t in (0.0, 0.01, 2.0, 30.0), z in range(0.0, L, length = 41)
            @test isfinite(analytic(ps, t, z))
        end
        @test analytic(ps, 2.0, 0.0) ≈ 0.0 atol = 1.0e-10
        @test analytic(ps, 2.0, L) ≈ 1.0 atol = 1.0e-10
        for (z, expected) in fdvalues
            @test analytic(ps, 2.0, z) ≈ expected atol = 1.0e-5
        end
        # Long-time limit is the steady state (1 - e^{-vz/k}) / (1 - e^{-vL/k}).
        a = -ps[2] / ps[1]
        steady(z) = expm1(a * z) / expm1(a * L)
        z = 0.3L
        a < 0 && @test analytic(ps, 100.0, z) ≈ steady(z) atol = 1.0e-8
    end
end
