using PDESystemLibrary
using PDEBase
using Test

@testset "MeshDomain heat equation with trace conditions" begin
    system = only(get_pdesys_with_tags(["MeshDomain"]))
    @test system.name == :heat_2d_meshdomain
    @test all(tag -> tag in system.metadata, ["2D", "Dirichlet", "Neumann", "Robin"])

    domain = only(filter(d -> d.domain isa MeshDomain, system.domain)).domain
    @test PDEBase.meshdim(domain.mesh) == 2
    @test Set(PDEBase.facetset_names(domain.mesh)) == Set([:left, :right, :bottom, :top])
    @test [0.5, 0.5] ∈ domain
    @test [1.5, 0.5] ∉ domain

    @test length(system.bcs) == 5   # four trace conditions, then the initial condition
    @test PDEBase.bckind.(system.bcs[1:4]) ==
        [PDEBase.Essential(), PDEBase.Essential(), PDEBase.Natural(), PDEBase.Robin()]
    @test_throws ArgumentError PDEBase.bckind(system.bcs[5])

    f = only(values(system.analytic_func))
    u(t, x, y) = f([], t, x, y)
    lambda = pi / 4
    alpha, beta = lambda * tan(lambda), 1.0
    h = 1.0e-4
    ut(t, x, y) = (u(t + h, x, y) - u(t - h, x, y)) / 2h
    ux(t, x, y) = (u(t, x + h, y) - u(t, x - h, y)) / 2h
    uy(t, x, y) = (u(t, x, y + h) - u(t, x, y - h)) / 2h
    lap(t, x, y) = (
        u(t, x + h, y) + u(t, x - h, y) + u(t, x, y + h) + u(t, x, y - h) - 4u(t, x, y)
    ) / h^2

    for t in 0.1:0.2:0.9, x in 0.1:0.2:0.9, y in 0.1:0.2:0.9
        @test ut(t, x, y) ≈ lap(t, x, y) atol = 1.0e-5
    end
    for t in 0.0:0.25:1.0, s in 0.0:0.25:1.0
        @test -ux(t, 0.0, s) ≈ 0 atol = 1.0e-8                      # left: Dn u = 0
        @test ux(t, 1.0, s) + alpha * u(t, 1.0, s) ≈ beta atol = 1.0e-8   # right: Robin
        @test u(t, s, 0.0) ≈ beta / alpha + exp(-2lambda^2 * t) * cos(lambda * s)   # bottom
        @test u(t, s, 1.0) ≈ beta / alpha + exp(-2lambda^2 * t) * cos(lambda * s) * cos(lambda)  # top
    end
end
