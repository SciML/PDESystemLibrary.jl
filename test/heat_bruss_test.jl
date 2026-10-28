using PDESystemLibrary
using ModelingToolkit: Num, build_function
using DomainSets: supremum
using Test

@testset "heat_1d1 solution is not negligible at the final time" begin
    system = only(filter(s -> s.name == :heat_1d1, get_pdesys_with_tags(String[])))
    D, D_value = only(getfield(system, :ps))
    t, x = getfield(system, :ivs)
    t_end = supremum(only(filter(d -> isequal(d.variables, t), system.domain)).domain)
    rhs = Num(only(getfield(system, :analytic)).rhs)
    exact = build_function(rhs, D, t, x; expression = Val{false})
    u_end(xi) = exact(D_value, t_end, xi)
    @test maximum(abs, u_end.(range(0, 1, length = 101))) > 1.0e-3
    @test t_end == 0.05
    @test u_end(0.25) ≈ exp(-4pi^2 * D_value * t_end)
end

@testset "bruss does not claim an analytic solution" begin
    system = only(filter(s -> s.name == :bruss, get_pdesys_with_tags(String[])))
    @test getfield(system, :analytic) === nothing
    @test "NumericalReference" in system.metadata
    @test !isempty(system.analytic_func)
end
