using PDESystemLibrary
using Test
using DomainSets

@testset "Precompiled analytic solutions" begin
    systems = get_pdesys_with_tags(String[])
    system = only(filter(system -> system.name == :inviscid_burgers_monotonic, systems))
    analytic = only(values(system.analytic_func))

    @test analytic([], 0.0, 0.5) == 0.5
end

@testset "Advection-dominated Gaussian" begin
    system = only(filter(s -> s.name == :advdom_gaussian, get_pdesys_with_tags(["AdvectionDominated"])))
    @test all(t -> t in system.metadata, ["1D", "Advection", "Diffusion", "AdvectionDominated"])

    pvals = Dict(string(first(p)) => last(p) for p in system.ps)
    v, ε = pvals["v"], pvals["eps"]
    xdom = only(filter(d -> string(d.variables) == "x", system.domain)).domain
    tdom = only(filter(d -> string(d.variables) == "t", system.domain)).domain
    L = DomainSets.supremum(xdom) - DomainSets.infimum(xdom)
    @test v * L / ε == 200  # documented Peclet number

    analytic = only(values(system.analytic_func))
    ps = [v, ε]
    ts = range(DomainSets.infimum(tdom), DomainSets.supremum(tdom); length = 11)
    xs = range(0.0, 1.0; length = 201)
    @test all(isfinite, [analytic(ps, t, x) for t in ts, x in xs])

    # PDE residual through central differences on the exact solution.
    h = 1.0e-4
    for t in ts[2:(end - 1)], x in xs[10:10:(end - 10)]
        ut = (analytic(ps, t + h, x) - analytic(ps, t - h, x)) / 2h
        ux = (analytic(ps, t, x + h) - analytic(ps, t, x - h)) / 2h
        uxx = (analytic(ps, t, x + h) - 2analytic(ps, t, x) + analytic(ps, t, x - h)) / h^2
        @test abs(ut + v * ux - ε * uxx) < 1.0e-4
    end

    # Initial and boundary conditions agree with the exact solution; the pulse starts at x0 = 0.2.
    @test analytic(ps, 0.0, 0.2) ≈ 1.0
    @test analytic(ps, 0.0, 0.0) < 1.0e-6
    @test analytic(ps, 0.5, 1.0) < 1.0e-2
    @test analytic(ps, 0.5, 0.7) ≈ sqrt(0.1 / 0.6)
end
