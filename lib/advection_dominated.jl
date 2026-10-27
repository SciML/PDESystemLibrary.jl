# Advection-dominated systems
"""
# Advection-Dominated Convection-Diffusion in 1D

Transient convection-diffusion of a travelling Gaussian pulse with Dirichlet boundary
conditions, a benchmark for stabilised (SUPG) finite element methods.

With velocity ``v = 1`` and diffusivity ``\\varepsilon = 0.005`` on a domain of length ``L = 1`` the
Peclet number is ``v L / \\varepsilon = 200``. The pulse is sharp enough that central
differences or an unstabilised Galerkin method oscillate on coarse meshes (the mesh
Peclet number ``v h / 2\\varepsilon`` exceeds 1 for ``h > 0.01``), yet the pulse remains resolvable.

```math
\\frac{\\partial u}{\\partial t} + v \\frac{\\partial u}{\\partial x} = \\varepsilon \\frac{\\partial^2 u}{\\partial x^2}
```

The exact solution is a Gaussian that translates at speed ``v`` and spreads diffusively,

```math
u(t, x) = \\sqrt{\\frac{t_0}{t + t_0}} \\exp\\left(-\\frac{(x - x_0 - v t)^2}{4 \\varepsilon (t + t_0)}\\right)
```

with ``x_0 = 0.2`` and ``t_0 = 0.1``. The initial and boundary data are taken from it.
"""
function advection_dominated_gaussian(ps, name = :advdom_gaussian; x0 = 0.2, t0 = 0.1)
    @variables u(..)
    @parameters t x v eps
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2

    eq = Dt(u(t, x)) + v * Dx(u(t, x)) ~ eps * Dxx(u(t, x))

    u_exact(t, x) = sqrt(t0 / (t + t0)) * exp(-(x - x0 - v * t)^2 / (4 * eps * (t + t0)))

    bcs = [
        u(0, x) ~ u_exact(0, x),
        u(t, 0) ~ u_exact(t, 0),
        u(t, 1) ~ u_exact(t, 1),
    ]

    domains = [
        t ∈ Interval(0.0, 0.5),
        x ∈ Interval(0.0, 1.0),
    ]

    analytic = [u(t, x) ~ u_exact(t, x)]

    tags = ["1D", "Dirichlet", "Linear", "Advection", "Diffusion", "AdvectionDominated"]

    return PDESystem(
        eq, bcs, domains, [t, x], [u(t, x)], [v => ps[1], eps => ps[2]],
        analytic = analytic, metadata = tags, name = name, eval_module = @__MODULE__
    )
end

# v = 1.0, eps = 0.005: Peclet number 200
push!(all_systems, advection_dominated_gaussian([1.0, 0.005], :advdom_gaussian))
