"""
# Spherical Laplacian 1D

The laplacian in spherical coordinates, with Dirichlet and Neumann boundary conditions.
Symmetrical in azimuth and elevation, so the solution is a function of only one spatial
variable.

```math
\\frac{\\partial^2 u}{\\partial t^2} - \\frac{4}{r^2} \\frac{\\partial}{\\partial r}(r^2 \\frac{\\partial}{\\partial r}(u) = 0
```
"""
function spherical_laplacian()
    # Parameters, variables, and derivatives
    @parameters t r
    @variables u(..)
    Dt = Differential(t)
    Dr = Differential(r)

    # 1D PDE and boundary conditions

    eq = Dt(u(t, r)) ~ 4 / r^2 * Dr(r^2 * Dr(u(t, r)))
    bcs = [u(0, r) ~ sin(r) / r,
        Dr(u(t, 0)) ~ 0,
        u(t, 1) ~ exp(-4t) * sin(1)]

    # Space and time domains
    domains = [t ∈ Interval(0.0, 1.0),
        r ∈ Interval(0.0, 1.0)]

    u_exact = [u(t, r) ~ exp.(-4t) * sin.(r) ./ r]

    tags = ["1D", "Dirichlet", "Neumann", "Spherical", "Diffusion"]

    # PDE system
    @named sph = PDESystem(eq, bcs, domains, [t, r], [u(t, r)], analytic = u_exact,
        metadata = tags)

    sph
end

push!(all_systems, spherical_laplacian())
