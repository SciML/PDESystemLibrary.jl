"""
# The Heat Equation in 1D with Dirichlet Boundary Conditions.

1D heat equation with Dirichlet boundary conditions.
This models the temperature of a rod over time, where the ends are held at a constant temperature.

It is initialized with a sinusoidal profile.
The equation is given by:

```math
\\frac{\\partial u}{\\partial t} = D \\frac{\\partial^2 u}{\\partial x^2}
```
"""
function heat_1d1()
    @variables x t u(..)
    @parameters D

    Dxx = Differential(x)
    Dt = Differential(t)

    eqs = [Dt(u(t, x)) ~ D * Dxx(u(t, x))]
    bcs = [u(0, x) ~ sin(2pi * x),
        u(t, 0) ~ 0.0, u(t, 1) ~ 0.0]

    domains = [t ∈ Interval(0.0, 1.0),
        x ∈ Interval(0.0, 1.0)]

    analytic = [u(t, x) ~ exp(-4pi^2 * D * t) * sin(2pi * x)]

    tags = ["1D", "Dirichlet", "Linear", "Diffusion", "Heat"]

    @named heat_1d1 = PDESystem(eqs, bcs, domains, [t, x], [u(t, x)], [D]; defaults = Dict(D => 1.0),
                                analytic = analytic, metadata = tags)

    heat_1d1
end
push!(all_systems, heat_1d1())

"""
# The Heat Equation in 1D with Neumann Boundary Conditions.

1D heat equation with Neumann boundary conditions.
This models the temperature of a rod over time, where the ends are held at a constant temperature.

It is initialized with a sinusoidal profile.
The equation is given by:

```math
\\frac{\\partial u}{\\partial t} = D \\frac{\\partial^2 u}{\\partial x^2}
```
"""
function heat_1d_neumann()
    # Method of Manufactured Solutions

    # Parameters, variables, and derivatives
    @parameters t x
    @variables u(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2

    # 1D PDE and boundary conditions
    eq = Dt(u(t, x)) ~ Dxx(u(t, x))
    bcs = [u(0, x) ~ cos(x),
        Dx(u(t, 0)) ~ 0,
        Dx(u(t, Float64(pi))) ~ 0]

    # Space and time domains
    domains = [t ∈ Interval(0.0, 1.0),
        x ∈ Interval(0.0, Float64(pi))]

    analytic = [u(t, x) ~ exp(-t) * cos(x)]

    tags = ["1D", "Neumann", "Linear", "Diffusion", "Heat"]
    # PDE system
    @named pdesys = PDESystem(eq, bcs, domains, [t, x], [u(t, x)], analytic = analytic,
                              metadata = tags)
end
push!(all_systems, heat_1d_neumann())
