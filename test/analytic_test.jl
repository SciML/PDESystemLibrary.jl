using PDESystemLibrary
using ModelingToolkit: Num, build_function
using Test

@testset "Precompiled analytic solutions" begin
    systems = get_pdesys_with_tags(String[])
    system = only(filter(system -> system.name == :inviscid_burgers_monotonic, systems))
    analytic = only(values(system.analytic_func))

    @test analytic([], 0.0, 0.5) == 0.5
end

@testset "heat_1d1 solution does not decay to zero" begin
    system = only(filter(s -> s.name == :heat_1d1, get_pdesys_with_tags(String[])))
    D, D_value = only(getfield(system, :ps))
    t, x = getfield(system, :ivs)
    t_end = 1.0
    rhs = Num(only(getfield(system, :analytic)).rhs)
    exact = build_function(rhs, D, t, x; expression = Val{false})
    u_end(xi) = exact(D_value, t_end, xi)
    @test maximum(abs, u_end.(range(0, 1, length = 101))) > 1.0e-3
    @test u_end(0.25) ≈ exp(-4pi^2 * D_value * t_end)
end
