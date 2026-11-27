using PDESystemLibrary
using PDEBase
using ModelingToolkit: Differential, Symbolics, unwrap, iscall, operation, arguments
using .Symbolics: wrap
using DomainSets: dimension
using Test

const MD_LAMBDA = pi / 4
const MD_ALPHA = MD_LAMBDA * tan(MD_LAMBDA)
const MD_BETA = 1.0
const MD_VARS = ("t", "x", "y")

# Partial derivative d^a/dt^a d^b/dx^b d^c/dy^c of the exact solution, written in closed form:
# the steady state beta/alpha plus exp(-2 lambda^2 t) cos(lambda x) cos(lambda y).
function exact_u(orders, t, x, y)
    a, b, c = orders
    steady = a + b + c == 0 ? MD_BETA / MD_ALPHA : 0.0
    return steady + (-2MD_LAMBDA^2)^a * exp(-2MD_LAMBDA^2 * t) *
        MD_LAMBDA^(b + c) * cos(MD_LAMBDA * x + b * pi / 2) * cos(MD_LAMBDA * y + c * pi / 2)
end

# Outward normal derivative on each named facet of the unit square.
const MD_FACETS = Dict(
    :left => (s, t) -> (0.0, s, -1, 0),
    :right => (s, t) -> (1.0, s, 1, 0),
    :bottom => (s, t) -> (s, 0.0, 0, -1),
    :top => (s, t) -> (s, 1.0, 0, 1),
)
facetpoint(name, s) = MD_FACETS[name](s, nothing)

# Evaluates a stored symbolic expression with the exact solution substituted for `u`.
# `s` parametrizes the facet and `t, x, y` give the point for the variables.
function evaluate(expr, env; orders = (0, 0, 0))
    e = unwrap(expr)
    Symbolics.value(e) isa Number && return Symbolics.value(e)
    iscall(e) || return env[string(e)]
    op, args = operation(e), arguments(e)
    if op isa Differential
        i = findfirst(==(string(op.x)), MD_VARS)
        return evaluate(args[1], env; orders = Tuple(orders .+ op.order .* (1:3 .== i)))
    elseif op isa PDEBase.NormalDerivative
        name = facetname(args[1])
        op.piece.name == name || error("normal derivative on $(op.piece.name) applied to a trace on $name")
        x, y, nx, ny = facetpoint(name, env["s"])
        t = env["t"]
        return nx * exact_u((orders[1], 1, 0), t, x, y) + ny * exact_u((orders[1], 0, 1), t, x, y)
    elseif string(op) == "u"
        facet = length(args) == 2 ? facetname(e) : nothing
        x, y = isnothing(facet) ? (evaluate(args[2], env), evaluate(args[3], env)) :
            facetpoint(facet, env["s"])[1:2]
        return exact_u(orders, evaluate(args[1], env), x, y)
    end
    return op((evaluate(a, env) for a in args)...)
end

# `u(t, Γ_name)` is the trace of `u` on the boundary piece named `name`.
facetname(trace) = Symbol(replace(string(arguments(unwrap(trace))[2]), "Γ_" => ""))

residual(eq, env) = evaluate(eq.lhs, env) - evaluate(eq.rhs, env)

# Largest residual of the PDE, each facet condition and the initial condition at asymmetric
# sample points. Returns residuals keyed by condition.
function residuals(eqs, bcs)
    samples = [0.3, 0.71]
    pde = maximum(
        abs(residual(eqs[1], Dict("t" => t, "x" => x, "y" => y)))
            for t in (0.13, 0.62), x in samples, y in (0.17, 0.89)
    )
    facets = [:bottom, :top, :left, :right]
    out = Dict{Symbol, Float64}(:pde => pde)
    for (name, bc) in zip(facets, bcs[1:4])
        out[name] = maximum(
            abs(residual(bc, Dict("t" => t, "s" => s, "x" => facetpoint(name, s)[1], "y" => facetpoint(name, s)[2])))
                for t in (0.13, 0.62), s in samples
        )
    end
    out[:initial] = maximum(
        abs(residual(bcs[5], Dict("t" => 0.0, "x" => x, "y" => y))) for x in samples, y in (0.17, 0.89)
    )
    return out
end

# Whether every call of the dependent variable `u` in a facet condition is a trace on the facet
# `name` at the free time variable, and every normal derivative is taken on that same facet.
function on_facet(domain, t, name, bc)
    piece = Symbolics.value(PDEBase.boundarypiece(domain, name))
    stored = only(filter(p -> p.name == name, PDEBase.boundarypieces(domain)))
    found = Ref(0)
    function check(expr)
        e = unwrap(expr)
        iscall(e) || return true
        op, args = operation(e), arguments(e)
        ok = all(check, args)
        if string(op) == "u"
            found[] += 1
            ok &= length(args) == 2 && isequal(args[1], unwrap(t)) && isequal(Symbolics.value(args[2]), piece)
        elseif op isa PDEBase.NormalDerivative
            found[] += 1
            ok &= op.piece == stored
        end
        return ok
    end
    return check(bc.lhs) && check(bc.rhs) && found[] >= 1
end

# Whether the initial condition is posed at time zero, as `u(0, x, y)`.
initial_at_zero(uf, x, y, bc) = isequal(arguments(unwrap(bc.lhs)), arguments(unwrap(uf(0, x, y))))

const MD_TOL = 1.0e-12

@testset "MeshDomain heat equation with trace conditions" begin
    system = only(get_pdesys_with_tags(["MeshDomain"]))
    @test system.name == :heat_2d_meshdomain
    @test all(tag -> tag in system.metadata, ["2D", "Dirichlet", "Neumann", "Robin"])

    pairing = only(filter(d -> d.domain isa MeshDomain, system.domain))
    domain = pairing.domain
    @test dimension(domain) == 2
    @test Set(p.name for p in PDEBase.boundarypieces(domain)) ==
        Set([:left, :right, :bottom, :top])
    @test length(PDEBase.boundarypieces(domain)) == 4
    @test [0.5, 0.5] ∈ domain
    @test [1.5, 0.5] ∉ domain
    mesh = domain.mesh
    @test PDEBase.contains_point(mesh, [0.5, 0.5])
    for p in (Float64[], [0.5], [0.5, 0.5, 0.5])
        @test !PDEBase.contains_point(mesh, p)
    end

    @test length(system.eqs) == 1
    @test length(system.bcs) == 5   # four trace conditions, then the initial condition
    @test PDEBase.bckind.(system.bcs[1:4]) ==
        [PDEBase.Essential(), PDEBase.Essential(), PDEBase.Natural(), PDEBase.Robin()]
    @test_throws ArgumentError PDEBase.bckind(system.bcs[5])

    f = only(values(system.analytic_func))
    for t in (0.13, 0.62), x in (0.3, 0.71), y in (0.17, 0.89)
        @test f([], t, x, y) ≈ exact_u((0, 0, 0), t, x, y) rtol = 1.0e-12
    end

    # The stored equations hold for the exact solution.
    res = residuals(system.eqs, system.bcs)
    for name in (:pde, :bottom, :top, :left, :right, :initial)
        @test res[name] < MD_TOL
    end

    # The residual check trusts that each condition sits at its facet and the initial time, so
    # verify the arguments of the dependent variable first.
    t, x, y = system.ivs
    u = only(system.dvs)
    uf = wrap(operation(unwrap(u)))
    # The facet geometry assumes the mesh coordinates are (x, y), in that order.
    mesh_order(p) = isequal(collect(arguments(p.variables)), unwrap.([x, y]))
    @test mesh_order(pairing)
    @test !mesh_order((y, x) ∈ domain)
    for (name, bc) in zip((:bottom, :top, :left, :right), system.bcs[1:4])
        @test on_facet(domain, t, name, bc)
    end
    @test initial_at_zero(uf, x, y, system.bcs[5])

    # The residual check cannot tell the diffusion coefficients apart, because the exact
    # solution has equal second derivatives in x and y, so compare the operator structurally.
    Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
    is_heat(eq) = isequal(eq.lhs, Dt(u)) && isequal(eq.rhs, Dx(Dx(u)) + Dy(Dy(u)))
    @test is_heat(only(system.eqs))

    # Each mutation of the stored model must be caught by its own residual.
    bcs = copy(system.bcs)
    robin_terms = arguments(unwrap(bcs[4].lhs))
    robin_u = wrap(only(filter(a -> string(operation(a)) == "*", robin_terms)))
    robin_dn = wrap(only(filter(a -> operation(a) isa PDEBase.NormalDerivative, robin_terms)))
    mutants = [
        (:pde, [Dt(u) ~ Dx(Dx(u)) - Dy(Dy(u))], bcs),
        (:top, system.eqs, [bcs[1:1]; bcs[1].rhs |> r -> bcs[2].lhs ~ r; bcs[3:5]]),
        (:right, system.eqs, [bcs[1:3]; robin_dn - robin_u ~ bcs[4].rhs; bcs[5:5]]),
        (:initial, system.eqs, [bcs[1:4]; bcs[5].lhs ~ 2MD_BETA / MD_ALPHA - bcs[5].rhs]),
    ]
    for (name, eqs, mutated) in mutants
        @test residuals(eqs, mutated)[name] > 1.0e-3
    end

    @test !is_heat(Dt(u) ~ 2 * Dx(Dx(u)))

    # A Robin datum that depends on y differs from the stored one away from y = 0.37.
    robin_y = [bcs[1:3]; bcs[4].lhs ~ 1 + y - 0.37; bcs[5:5]]
    @test residuals(system.eqs, robin_y)[:right] > 0.3

    # Location mutants have the right residuals but sit on the wrong facet or time.
    Γ_bottom = PDEBase.boundarypiece(domain, :bottom)
    Dn_bottom = PDEBase.NormalDerivative(Γ_bottom)
    neumann_on_bottom = Dn_bottom(uf(t, Γ_bottom)) ~ 0
    @test !on_facet(domain, t, :left, neumann_on_bottom)
    @test !initial_at_zero(uf, x, y, uf(t, x, y) ~ bcs[5].rhs)
end
