using PDESystemLibrary
using Test
import ModelingToolkit

@testset "Spherical laplacian is finite at r = 0" begin
    system = only(filter(system -> system.name == :sph, get_pdesys_with_tags(["Spherical"])))
    analytic = only(values(system.analytic_func))
    initial_condition = system.bcs[1].rhs
    r = system.ivs[2]

    @test isfinite(analytic([], 0.0, 0.0))
    @test analytic([], 0.0, 0.0) == 1.0
    @test analytic([], 0.5, 0.0) ≈ exp(-2.0)
    @test analytic([], 0.5, 1.0) ≈ exp(-2.0) * sin(1.0)
    ic_func = ModelingToolkit.Symbolics.build_function(initial_condition, r; expression = Val{false})
    @test ic_func(0.0) == 1
    @test ic_func(1.0) ≈ sin(1.0)
end
