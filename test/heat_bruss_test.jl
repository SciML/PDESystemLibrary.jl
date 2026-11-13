using PDESystemLibrary
using ModelingToolkit: Num, build_function
using DomainSets: supremum
using Test

# Returns the final time and a function evaluating the exact solution at that time.
function final_time_solution(name)
    system = only(filter(s -> s.name == name, get_pdesys_with_tags(String[])))
    D, D_value = only(getfield(system, :ps))
    t, x = getfield(system, :ivs)
    t_end = supremum(only(filter(d -> isequal(d.variables, t), system.domain)).domain)
    rhs = Num(only(getfield(system, :analytic)).rhs)
    exact = build_function(rhs, D, t, x; expression = Val{false})
    return t_end, xi -> exact(D_value, t_end, xi), D_value
end

@testset "heat_1d1 keeps its original interval and exact solution" begin
    t_end, u_end, D_value = final_time_solution(:heat_1d1)
    @test t_end == 1.0
    @test D_value == 1.0
    @test u_end(0.25) ≈ exp(-4pi^2 * D_value * t_end)
end

@testset "heat_1d1_short solution is not negligible at the final time" begin
    t_end, u_end, D_value = final_time_solution(:heat_1d1_short)
    @test t_end == 0.05
    @test maximum(abs, u_end.(range(0, 1, length = 101))) > 1.0e-3
    @test u_end(0.25) ≈ exp(-4pi^2 * D_value * t_end)
    @test "ShortTime" in only(filter(s -> s.name == :heat_1d1_short, get_pdesys_with_tags(String[]))).metadata
end

@testset "bruss does not claim an analytic solution" begin
    system = only(filter(s -> s.name == :bruss, get_pdesys_with_tags(String[])))
    @test getfield(system, :analytic) === nothing
    @test "NumericalReference" in system.metadata
    @test !isempty(system.analytic_func)
end
