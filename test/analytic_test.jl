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
    # Same solve at small t, where the series needs many terms: (t, z, value).
    smallt = Dict(
        :convdiff1 => [(0.01, 0.05, 0.0), (0.01, 0.95, 0.0039782), (0.1, 0.95, 0.92730929)],
        :convdiff4 => [(0.01, 9.9, 0.4183456), (0.1, 9.0, 0.0213762)],
    )
    nearL = Dict(
        :convdiff1 => [(0.001, 0.86178896), (0.005, 0.33369459)],
        :convdiff2 => [(0.001, 0.98166495), (0.005, 0.90868539)],
        :convdiff3 => [(0.001, 0.68385034), (0.005, 0.10986312)],
        :convdiff4 => [(0.001, 0.97759218), (0.005, 0.88711224)],
    )
    for (name, ps, L, fdvalues) in cases
        system = only(filter(system -> system.name == name, systems))
        analytic = only(values(system.analytic_func))
        for t in (0.0, 0.01, 2.0, 30.0), z in range(0.0, L, length = 41)
            @test isfinite(analytic(ps, t, z))
        end
        # Initial condition, and a solution that stays within [0, 1] at small t.
        @test analytic(ps, 0.0, 0.5L) == 0.0
        @test analytic(ps, 0.0, L) == 1.0
        for t in (1.0e-4, 0.01, 0.1), z in range(0.0, L, length = 41)
            @test -1.0e-12 <= analytic(ps, t, z) <= 1 + 1.0e-12
        end
        @test analytic(ps, 2.0, 0.0) ≈ 0.0 atol = 1.0e-10
        @test analytic(ps, 2.0, L) ≈ 1.0 atol = 1.0e-10
        for (z, expected) in fdvalues
            @test analytic(ps, 2.0, z) ≈ expected atol = 1.0e-5
        end
        for (t, z, expected) in get(smallt, name, ())
            @test analytic(ps, t, z) ≈ expected atol = 1.0e-6
        end
        # Boundary layer at t = 1e-3 near z = L (Crank-Nicolson, dz = L / 160000).
        for (off, expected) in nearL[name]
            @test analytic(ps, 1.0e-3, L - off) ≈ expected atol = 2.0e-5
        end
        # Too small t is a loud failure, never a silently truncated series.
        @test_throws ArgumentError analytic(ps, 1.0e-40, L - 1.0e-5)
        # Long-time limit is the steady state (1 - e^{-vz/k}) / (1 - e^{-vL/k}).
        a = -ps[2] / ps[1]
        steady(z) = expm1(a * z) / expm1(a * L)
        z = 0.3L
        a < 0 && @test analytic(ps, 100.0, z) ≈ steady(z) atol = 1.0e-8
    end
end
