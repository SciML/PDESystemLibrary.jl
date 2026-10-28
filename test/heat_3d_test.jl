using PDESystemLibrary
using Test
using DomainSets

@testset "heat_3d analytic solution" begin
    system = only(get_pdesys_with_tags(["3D", "Heat", "Dirichlet"]))
    @test system.name == :heat_3d
    kval = only(system.ps)[2]
    analytic = only(values(system.analytic_func))
    u(t, x, y, z) = analytic([kval], t, x, y, z)

    # Residual of Dt(u) - k*(Dxx+Dyy+Dzz)(u) by central differences (O(h^2) error).
    h = 1.0e-3
    for (t, x, y, z) in [(0.1, 0.3, 0.4, 0.5), (0.2, 0.7, 0.1, 0.9), (0.4, 0.25, 0.5, 0.75)]
        ut = (u(t + h, x, y, z) - u(t - h, x, y, z)) / 2h
        lap = (
            u(t, x + h, y, z) + u(t, x - h, y, z) + u(t, x, y + h, z) +
                u(t, x, y - h, z) + u(t, x, y, z + h) + u(t, x, y, z - h) -
                6 * u(t, x, y, z)
        ) / h^2
        @test abs(ut - kval * lap) < 1.0e-4
    end

    @test u(0.0, 0.3, 0.6, 0.2) ≈ sin(pi * 0.3) * sin(pi * 0.6) * sin(pi * 0.2)
    for t in (0.0, 0.2, 0.5), a in (0.3, 0.8), b in (0.1, 0.6), face in (0.0, 1.0)
        @test abs(u(t, face, a, b)) < 1.0e-12
        @test abs(u(t, a, face, b)) < 1.0e-12
        @test abs(u(t, a, b, face)) < 1.0e-12
    end

    t_end = DomainSets.supremum(system.domain[1].domain)
    @test u(t_end, 0.5, 0.5, 0.5) > 0.1
end
