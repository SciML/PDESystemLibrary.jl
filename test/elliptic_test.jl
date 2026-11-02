using PDESystemLibrary
using Test

@testset "Steady elliptic analytic solution: laplace_2d_mixed" begin
    using Symbolics: Differential, expand_derivatives, substitute, unwrap, arguments, operation, build_function
    system = only(filter(s -> s.name == :laplace_2d_mixed, get_pdesys_with_tags(String[])))
    @test all(t -> t in system.metadata, ["2D", "Steady", "Elliptic"])

    x, y = system.ivs
    exact = only(system.analytic).rhs
    points = [(0.0, 0.0), (0.3, 0.7), (0.5, 0.5), (1.0, 1.0), (0.9, 0.1), (0.25, 1.0)]
    function vanishes(residual)
        f = build_function(residual, x, y; expression = Val{false})
        return all(((xv, yv),) -> abs(f(xv, yv)) < 1.0e-12, points)
    end

    # PDE: the exact solution's derivatives satisfy the equation.
    eq = only(system.eqs)
    derivative(expr, D) = expand_derivatives(D(expr))
    laplacian = sum(D -> derivative(exact, D), [Differential(x)^2, Differential(y)^2])
    @test vanishes(laplacian - eq.rhs)

    # BCs: the exact solution (or its derivative, for Neumann edges) matches each edge.
    for bc in system.bcs
        lhs = unwrap(bc.lhs)
        if operation(lhs) isa Differential
            call = only(arguments(lhs))
            value_at_edge = derivative(exact, operation(lhs))
        else
            call = lhs
            value_at_edge = exact
        end
        edge = Dict(zip([x, y], arguments(call)))
        @test vanishes(substitute(value_at_edge, edge) - bc.rhs)
    end

    # The precompiled analytic function agrees with the closed form.
    analytic = only(values(system.analytic_func))
    for (xv, yv) in points
        @test analytic([], xv, yv) ≈ sin(pi * xv) * sinh(pi * yv) / sinh(pi)
    end
end
