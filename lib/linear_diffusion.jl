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
begin
    @variables x t
    @parameters D

    Dxx = Differential(x)
    Dt = Differential(t)

    eqs = [Dt(u(t, x)) ~ D * Dxx(u(t, x))]
    bcs = [u(0, x) ~ sin(2pi * x),
           u(t, 0) ~ 0.0, u(t, 1) ~ 0.0]

    domains = [t ∈ IntervalDomain(0.0, 1.0),
               x ∈ IntervalDomain(0.0, 1.0)]

    analytic = [u(t, x) ~ exp(-4pi^2 * D * t) * sin(2pi * x)]

    tags = ["1D", "Dirichlet", "Linear", "Diffusion", "Heat"]

    @named heat_1d1 = PDESystem(eqs, bcs, domains, [t, x], [u(t, x)], [D => 1.0],
                               analytic = analytic, metadata = tags)

    heat_1d
end

push!(all_systems, heat_1d)
