using PDESystemLibrary
using Test

@testset "Robin heat equation analytic solution" begin
    system = only(get_pdesys_with_tags(["1D", "Robin"]))
    @test system.name == :heat_1d_robin
    f = only(values(system.analytic_func))
    u(t, x) = f([], t, x)

    alpha, beta = pi / 4 * tan(pi / 4), 1.0
    h = 1.0e-4
    ut(t, x) = (u(t + h, x) - u(t - h, x)) / (2h)
    ux(t, x) = (u(t, x + h) - u(t, x - h)) / (2h)
    uxx(t, x) = (u(t, x + h) - 2u(t, x) + u(t, x - h)) / h^2

    for t in 0.1:0.2:0.9, x in 0.1:0.2:0.9
        @test ut(t, x) ≈ uxx(t, x) atol = 1.0e-5
    end
    for t in 0.0:0.25:1.0
        @test ux(t, 0.0) ≈ 0 atol = 1.0e-8
        @test ux(t, 1.0) + alpha * u(t, 1.0) ≈ beta atol = 1.0e-8
    end
    for x in 0.0:0.25:1.0
        @test u(0.0, x) ≈ beta / alpha + cos(pi / 4 * x)
    end
end
