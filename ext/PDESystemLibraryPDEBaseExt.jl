module PDESystemLibraryPDEBaseExt

using IntervalSets: Interval
using ModelingToolkitBase: @named, @parameters, PDESystem
using PDEBase: MeshDomain, NormalDerivative, boundarypiece
using PDESystemLibrary: PDESystemLibrary, all_systems
using Symbolics: @variables, Differential
import PDEBase

"""
    UnitSquareMesh()

The smallest mesh implementing the PDEBase mesh-domain interface: the unit square
``[0, 1]^2`` with the named boundary facet sets `:left`, `:right`, `:bottom` and `:top`.
It stores no cells, so it only supports declaring a `PDESystem` and evaluating its exact
solution; a discretizer needs a backend mesh.
"""
struct UnitSquareMesh end

PDEBase.meshdim(::UnitSquareMesh) = 2
PDEBase.facetset_names(::UnitSquareMesh) = (:left, :right, :bottom, :top)
PDEBase.cellset_names(::UnitSquareMesh) = (:all,)
PDEBase.contains_point(::UnitSquareMesh, p) = all(c -> 0 <= c <= 1, p)

"""
# The Heat Equation on a `MeshDomain` with Trace Boundary Conditions.

2D heat equation on the unit square, declared on a PDEBase `MeshDomain` whose boundary
conditions are written as conditions on the traces of `u` on named boundary pieces:
a Neumann condition on `left`, a Robin condition on `right`, and Dirichlet conditions on
`bottom` and `top`:

```math
\\frac{\\partial u}{\\partial t} = \\Delta u, \\quad
\\partial_n u|_{\\mathrm{left}} = 0, \\quad
\\partial_n u + \\alpha u|_{\\mathrm{right}} = \\beta, \\quad
u|_{\\mathrm{bottom}} = g_b(t, x), \\quad u|_{\\mathrm{top}} = g_t(t, x)
```

where ``\\partial_n`` is the outward normal derivative, ``\\lambda = \\pi/4``,
``\\alpha = \\lambda \\tan\\lambda`` and ``\\beta = 1``. The exact solution is the steady
state ``\\beta/\\alpha`` plus a decaying product mode with decay rate ``2\\lambda^2``;
``\\cos(\\lambda x)`` has zero slope at ``x = 0`` and satisfies the homogeneous Robin condition
at ``x = 1`` because ``\\alpha = \\lambda \\tan\\lambda``:

```math
u(t, x, y) = \\frac{\\beta}{\\alpha} + e^{-2\\lambda^2 t} \\cos(\\lambda x) \\cos(\\lambda y)
```

so ``g_b = \\beta/\\alpha + e^{-2\\lambda^2 t}\\cos(\\lambda x)`` and
``g_t = \\beta/\\alpha + e^{-2\\lambda^2 t}\\cos(\\lambda x)\\cos\\lambda``.
"""
function heat_2d_meshdomain(mesh = UnitSquareMesh())
    @parameters t x y
    @variables u(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)

    Ω = MeshDomain(mesh)
    Γ_left, Γ_right, Γ_bottom, Γ_top = boundarypiece.(Ref(Ω), (:left, :right, :bottom, :top))
    Dn_left, Dn_right = NormalDerivative(Γ_left), NormalDerivative(Γ_right)

    lambda = pi / 4
    alpha = lambda * tan(lambda)
    beta = 1.0

    eq = Dt(u(t, x, y)) ~ Dx(Dx(u(t, x, y))) + Dy(Dy(u(t, x, y)))
    bcs = [
        u(t, Γ_bottom) ~ beta / alpha + exp(-2 * lambda^2 * t) * cos(lambda * x),
        u(t, Γ_top) ~ beta / alpha + exp(-2 * lambda^2 * t) * cos(lambda * x) * cos(lambda),
        Dn_left(u(t, Γ_left)) ~ 0,
        Dn_right(u(t, Γ_right)) + alpha * u(t, Γ_right) ~ beta,
        u(0, x, y) ~ beta / alpha + cos(lambda * x) * cos(lambda * y),
    ]

    domains = [(x, y) ∈ Ω, t ∈ Interval(0.0, 1.0)]

    analytic = [u(t, x, y) ~ beta / alpha + exp(-2 * lambda^2 * t) * cos(lambda * x) * cos(lambda * y)]

    tags = ["2D", "MeshDomain", "Dirichlet", "Neumann", "Robin", "Linear", "Diffusion", "Heat"]

    return @named heat_2d_meshdomain = PDESystem(
        eq, bcs, domains, [t, x, y], [u(t, x, y)],
        analytic = analytic, metadata = tags, eval_module = PDESystemLibrary
    )
end

# Registered at load time: mutations of `all_systems` made while precompiling this
# extension would not survive into the session.
__init__() = push!(all_systems, heat_2d_meshdomain())

end
