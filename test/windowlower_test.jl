using PDESystemLibrary
using Test

@testset "windowlower and dconv analytic solutions" begin
    windowlower = PDESystemLibrary.windowlower
    @test windowlower(5.0, -0.1) == 0.0
    @test windowlower(5.0, 0.0) == 5.0
    @test windowlower(5.0, 0.5) == 5.0
    @test windowlower(5.0, 1.0) == 5.0
    @test windowlower(5.0, 1.1) == 0.0

    systems = get_pdesys_with_tags(String[])
    system = first(filter(system -> system.name == :dconvcos, systems))
    analytic = only(values(system.analytic_func))
    @test analytic([2], 0.0, 0.25) ≈ cos(2pi * 0.25) atol = 1.0e-12
    @test analytic([2], 0.0, 0.1) ≈ cos(2pi * 0.1)
    @test analytic([2], 0.25, 0.8) ≈ cos(2pi * 0.3)  # x - v*t = 0.3
    @test analytic([2], 0.25, 0.3) == 0.0            # x - v*t < 0: inflow Dirichlet BC
end
