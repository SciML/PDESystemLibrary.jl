"""
# The Heat Equation in 3D with Dirichlet Boundary Conditions.

3D heat equation on the unit cube with homogeneous Dirichlet boundary conditions.
This models the temperature of a cube whose faces are held at zero temperature.

It is initialized with the lowest sinusoidal mode, so the exact solution decays as
`exp(-3pi^2 k t)`. With `k = 0.1` and `t ∈ [0, 0.5]` the solution at the final time is
about 23% of the initial amplitude.
The equation is given by:

```math
\\frac{\\partial u}{\\partial t} = k \\left(\\frac{\\partial^2 u}{\\partial x^2} +
\\frac{\\partial^2 u}{\\partial y^2} + \\frac{\\partial^2 u}{\\partial z^2}\\right)
```
"""
function heat_3d()
    @parameters x y z t k
    @variables u(..)

    Dxx = Differential(x)^2
    Dyy = Differential(y)^2
    Dzz = Differential(z)^2
    Dt = Differential(t)

    eqs = [Dt(u(t, x, y, z)) ~ k * (Dxx(u(t, x, y, z)) + Dyy(u(t, x, y, z)) + Dzz(u(t, x, y, z)))]
    bcs = [
        u(0, x, y, z) ~ sin(pi * x) * sin(pi * y) * sin(pi * z),
        u(t, 0, y, z) ~ 0.0, u(t, 1, y, z) ~ 0.0,
        u(t, x, 0, z) ~ 0.0, u(t, x, 1, z) ~ 0.0,
        u(t, x, y, 0) ~ 0.0, u(t, x, y, 1) ~ 0.0,
    ]

    domains = [
        t ∈ Interval(0.0, 0.5),
        x ∈ Interval(0.0, 1.0),
        y ∈ Interval(0.0, 1.0),
        z ∈ Interval(0.0, 1.0),
    ]

    analytic = [u(t, x, y, z) ~ exp(-3pi^2 * k * t) * sin(pi * x) * sin(pi * y) * sin(pi * z)]

    tags = ["3D", "Dirichlet", "Linear", "Diffusion", "Heat"]

    @named heat_3d = PDESystem(
        eqs, bcs, domains, [t, x, y, z], [u(t, x, y, z)], [k => 0.1],
        analytic = analytic, metadata = tags, eval_module = @__MODULE__
    )

    return heat_3d
end
push!(all_systems, heat_3d())
