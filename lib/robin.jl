"""
# The Heat Equation in 1D with a Robin Boundary Condition.

1D heat equation on ``[0, 1]`` with a Neumann condition at ``x = 0`` and a Robin
(convective) condition at ``x = 1``, which models a rod whose end exchanges heat with
an environment:

```math
\\frac{\\partial u}{\\partial t} = \\frac{\\partial^2 u}{\\partial x^2}, \\quad
\\frac{\\partial u}{\\partial x}(t, 0) = 0, \\quad
\\frac{\\partial u}{\\partial x}(t, 1) + \\alpha u(t, 1) = \\beta
```

with ``\\alpha = \\pi/4`` and ``\\beta = 1``. The exact solution is the steady state
``\\beta/\\alpha`` plus a decaying mode ``\\cos(\\lambda x)`` with ``\\lambda = \\pi/4``;
``\\alpha = \\lambda \\tan\\lambda`` is exactly the condition that makes the mode satisfy the
homogeneous Robin condition:

```math
u(t, x) = \\frac{\\beta}{\\alpha} + e^{-\\lambda^2 t} \\cos(\\lambda x)
```
"""
function heat_1d_robin()
    @parameters t x
    @variables u(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2

    lambda = pi / 4
    alpha = lambda * tan(lambda)
    beta = 1.0

    eq = Dt(u(t, x)) ~ Dxx(u(t, x))
    bcs = [
        u(0, x) ~ beta / alpha + cos(lambda * x),
        Dx(u(t, 0)) ~ 0,
        Dx(u(t, 1)) + alpha * u(t, 1) ~ beta,
    ]

    domains = [
        t ∈ Interval(0.0, 1.0),
        x ∈ Interval(0.0, 1.0),
    ]

    analytic = [u(t, x) ~ beta / alpha + exp(-lambda^2 * t) * cos(lambda * x)]

    tags = ["1D", "Neumann", "Robin", "Linear", "Diffusion", "Heat"]

    return @named heat_1d_robin = PDESystem(
        eq, bcs, domains, [t, x], [u(t, x)],
        analytic = analytic, metadata = tags, eval_module = @__MODULE__
    )
end
push!(all_systems, heat_1d_robin())
