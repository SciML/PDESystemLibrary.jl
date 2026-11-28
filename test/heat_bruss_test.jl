using PDESystemLibrary
using Symbolics: Num, build_function
using DomainSets: supremum
using Test

heat_system(tags; kwargs...) = only(get_pdesys_with_tags(tags; kwargs...))

# Checks the default diffusivity and the stored analytic callback against exp(-4π² t) sin(2π x).
function check_heat_solution(system, t_end)
    t, _ = getfield(system, :ivs)
    t_domain = only(filter(d -> isequal(d.variables, t), system.domain)).domain
    @test supremum(t_domain) == t_end

    D_value = last(only(getfield(system, :ps)))
    @test D_value == 1.0
    callback = only(values(system.analytic_func))
    for time in (0.0, t_end), xi in (0.0, 0.25, 0.5, 0.75, 1.0)
        @test callback([D_value], time, xi) ≈ exp(-4pi^2 * time) * sinpi(2 * xi) atol = 1.0e-12
    end
    return callback, D_value
end

@testset "heat_1d1 keeps its original interval and exact solution" begin
    system = heat_system(["Heat", "Dirichlet", "1D"]; without = ["ShortTime"])
    @test system.name == :heat_1d1
    check_heat_solution(system, 1.0)
end

@testset "heat_1d1_short solution is not negligible at the final time" begin
    system = heat_system(["Heat", "ShortTime"])
    @test system.name == :heat_1d1_short
    callback, D_value = check_heat_solution(system, 0.05)
    @test maximum(abs, [callback([D_value], 0.05, xi) for xi in range(0, 1, length = 101)]) > 1.0e-3
    @test callback([D_value], 0.05, 0.25) ≈ exp(-4pi^2 * 0.05) atol = 1.0e-12

    D, _ = only(getfield(system, :ps))
    t, x = getfield(system, :ivs)
    rhs = Num(only(getfield(system, :analytic)).rhs)
    exact = build_function(rhs, D, t, x; expression = Val{false})
    @test exact(D_value, 0.05, 0.25) ≈ callback([D_value], 0.05, 0.25)
end

@testset "heat_1d1 and heat_1d1_short do not share tag vectors" begin
    long = only(filter(s -> s.name == :heat_1d1, get_pdesys_with_tags(String[])))
    short = only(filter(s -> s.name == :heat_1d1_short, get_pdesys_with_tags(String[])))
    @test long.metadata !== short.metadata
    @test "ShortTime" ∉ long.metadata
    @test "ShortTime" ∈ short.metadata
end

@testset "bruss does not claim an analytic solution" begin
    system = only(filter(s -> s.name == :bruss, get_pdesys_with_tags(String[])))
    @test getfield(system, :analytic) === nothing
    @test "NumericalReference" in system.metadata
    @test !isempty(system.analytic_func)
end
