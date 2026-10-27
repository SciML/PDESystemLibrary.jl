using PDESystemLibrary
using Test

@testset "Precompiled analytic solutions" begin
    systems = get_pdesys_with_tags(String[])
    system = only(filter(system -> system.name == :inviscid_burgers_monotonic, systems))
    analytic = only(values(system.analytic_func))

    @test analytic([], 0.0, 0.5) == 0.5
end

@testset "funcconv velocities are deterministic" begin
    build() = begin
        systems = []
        PDESystemLibrary.add_systems!(
            systems, PDESystemLibrary.funcs,
            PDESystemLibrary.linear_convection_dirichlet3, "funcconv"
        )
        systems
    end
    velocity(system) = last(only(system.ps))
    @test length(build()) == 4 * length(PDESystemLibrary.funcs)^2
    @test map(velocity, build()) == map(velocity, build())
    @test velocity(first(build())) == mod(0.6180339887498949, 1.0)

    library = get_pdesys_with_tags(String[])
    funcconv = filter(s -> startswith(string(s.name), "funcconv"), library)
    @test length(funcconv) == 2 * 4 * length(PDESystemLibrary.funcs)^2
end
