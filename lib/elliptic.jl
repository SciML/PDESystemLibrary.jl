"""
# The Laplace Equation in 2D with Mixed Boundary Conditions.

Steady 2D Laplace equation on the unit square with Dirichlet conditions on three edges and a
Neumann condition on the top edge. It has no time variable.
The equation is given by:

```math
\\frac{\\partial^2 u}{\\partial x^2} + \\frac{\\partial^2 u}{\\partial y^2} = 0
```

The exact solution is ``u = \\sin(\\pi x) \\sinh(\\pi y) / \\sinh(\\pi)``.
"""
function laplace_2d_mixed()
    @parameters x y
    @variables u(..)

    Dy = Differential(y)
    Dxx = Differential(x)^2
    Dyy = Differential(y)^2

    eq = Dxx(u(x, y)) + Dyy(u(x, y)) ~ 0.0
    bcs = [
        u(0, y) ~ 0.0,
        u(1, y) ~ 0.0,
        u(x, 0) ~ 0.0,
        Dy(u(x, 1)) ~ pi * sin(pi * x) * cosh(pi) / sinh(pi),
    ]

    domains = [
        x ∈ Interval(0.0, 1.0),
        y ∈ Interval(0.0, 1.0),
    ]

    analytic = [u(x, y) ~ sin(pi * x) * sinh(pi * y) / sinh(pi)]

    tags = ["2D", "Steady", "Elliptic", "Linear", "Dirichlet", "Neumann", "Laplace"]

    @named laplace_2d_mixed = PDESystem(
        eq, bcs, domains, [x, y], [u(x, y)],
        analytic = analytic, metadata = tags, eval_module = @__MODULE__
    )

    return laplace_2d_mixed
end
push!(all_systems, laplace_2d_mixed())
