using PDESystemLibrary
using Test
using DomainSets
using Symbolics

@testset "heat_3d analytic solution" begin
    system = only(get_pdesys_with_tags(["3D", "Heat", "Dirichlet"]))
    @test system.name == :heat_3d
    kval = only(system.ps)[2]
    @test kval == 0.1  # documented default diffusivity
    analytic = only(values(system.analytic_func))
    u(t, x, y, z) = analytic([kval], t, x, y, z)

    # Residual of Dt(u) - k*(Dxx+Dyy+Dzz)(u) by central differences (O(h^2) error).
    h = 1.0e-3
    for (t, x, y, z) in [(0.1, 0.3, 0.4, 0.5), (0.2, 0.7, 0.1, 0.9), (0.4, 0.25, 0.5, 0.75)]
        ut = (u(t + h, x, y, z) - u(t - h, x, y, z)) / 2h
        lap = (
            u(t, x + h, y, z) + u(t, x - h, y, z) + u(t, x, y + h, z) +
                u(t, x, y - h, z) + u(t, x, y, z + h) + u(t, x, y, z - h) -
                6 * u(t, x, y, z)
        ) / h^2
        @test abs(ut - kval * lap) < 1.0e-4
    end

    @test u(0.0, 0.3, 0.6, 0.2) ≈ sin(pi * 0.3) * sin(pi * 0.6) * sin(pi * 0.2)
    for t in (0.0, 0.2, 0.5), a in (0.3, 0.8), b in (0.1, 0.6), face in (0.0, 1.0)
        @test abs(u(t, face, a, b)) < 1.0e-12
        @test abs(u(t, a, face, b)) < 1.0e-12
        @test abs(u(t, a, b, face)) < 1.0e-12
    end

    t_end = DomainSets.supremum(system.domain[1].domain)
    @test t_end == 0.5  # documented final time
    @test analytic([0.2], t_end, 0.5, 0.5, 0.5) ≈ exp(-3pi^2 * 0.2 * t_end) rtol = 1.0e-12  # k is read from the parameter
    # Centre value is exp(-3pi^2 k t_end) * sin(pi/2)^3, the documented lowest mode decay.
    @test u(t_end, 0.5, 0.5, 0.5) ≈ exp(-3pi^2 * 0.1 * t_end) rtol = 1.0e-12

    # The stored equation and boundary conditions are satisfied by the exact solution.
    t, x, y, z = system.ivs
    pmap = Dict(first(p) => last(p) for p in system.ps)
    u_sym = only(system.dvs)
    exact = only(system.analytic).rhs
    pde = only(system.eqs)
    Dt, Dx, Dy, Dz = Differential(t), Differential(x), Differential(y), Differential(z)
    derivs = Dict(
        d(u_sym) => Symbolics.expand_derivatives(d(exact))
            for d in (Dt, Dx^2, Dy^2, Dz^2)
    )
    pde_res = substitute(pde.lhs - pde.rhs, derivs)
    samples = [(0.1, 0.3, 0.4, 0.5), (0.2, 0.7, 0.1, 0.9), (0.4, 0.25, 0.5, 0.75)]
    for (tv, xv, yv, zv) in samples
        vals = merge(pmap, Dict(t => tv, x => xv, y => yv, z => zv))
        @test abs(Symbolics.value(substitute(pde_res, vals; fold = Val(true)))) < 1.0e-8
    end

    # Trial fields v = exp(-t) sin(a pi x) sin(b pi y) sin(c pi z) have Dt(v) = -v and
    # (Dxx + Dyy + Dzz)(v) = -(a^2 + b^2 + c^2) pi^2 v, so the stored residual is
    # (-1 + (a^2 + b^2 + c^2) pi^2 k) v only if each axis carries weight exactly 1. The three
    # frequency triples have linearly independent squares, which pins all three weights. It is
    # checked at the default diffusivity and at a second one supplied through the parameter, so a
    # literal k in the stored equation cannot pass.
    k_sym = only(system.ps)[1]
    for (a, b, c) in ((1, 2, 3), (2, 1, 1), (1, 1, 2))
        trial = exp(-t) * sin(a * pi * x) * sin(b * pi * y) * sin(c * pi * z)
        trial_derivs = Dict(d(u_sym) => Symbolics.expand_derivatives(d(trial)) for d in (Dt, Dx^2, Dy^2, Dz^2))
        trial_res = substitute(pde.lhs - pde.rhs, trial_derivs)
        for ksup in (kval, 0.2), (tv, xv, yv, zv) in samples
            vals = merge(pmap, Dict(k_sym => ksup, t => tv, x => xv, y => yv, z => zv))
            vval = exp(-tv) * sin(a * pi * xv) * sin(b * pi * yv) * sin(c * pi * zv)
            @test Symbolics.value(substitute(trial_res, vals; fold = Val(true))) ≈ (-1 + (a^2 + b^2 + c^2) * pi^2 * ksup) * vval atol = 1.0e-10
        end
    end

    # Each boundary condition sits on the initial time or one end of an x, y or z interval,
    # and its right-hand side matches the exact solution there.
    doms = Dict(string(d.variables) => d.domain for d in system.domain)
    ends = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]  # the unit cube, stated rather than read from the stored domains
    @test all(v -> (DomainSets.infimum(doms[v]), DomainSets.supremum(doms[v])) == (0.0, 1.0), ("x", "y", "z"))
    @test DomainSets.infimum(doms["t"]) == 0.0
    @test length(system.bcs) == 7
    signatures = Vector{Vector{Tuple{Int, Float64}}}()
    for bc in system.bcs
        args = Symbolics.arguments(Symbolics.unwrap(bc.lhs))
        fixed = [isempty(Symbolics.get_variables(a)) ? Symbolics.value(a) : nothing for a in args]
        # The whole left-hand side is a call of the registered dependent variable on these arguments.
        @test isequal(Symbolics.unwrap(bc.lhs), Symbolics.operation(Symbolics.unwrap(u_sym))(args...))
        push!(signatures, [(i, Float64(f)) for (i, f) in enumerate(fixed) if f !== nothing])
        # Every argument that is not fixed is its own independent variable, so the whole face is constrained.
        @test isequal(args[isnothing.(fixed)], Symbolics.unwrap.(system.ivs)[isnothing.(fixed)])
        @test fixed[1] == DomainSets.infimum(doms["t"]) ||
            any(i -> fixed[i + 1] in ends[i], 1:3)
        for (tv, xv, yv, zv) in samples
            coords = [Symbolics.value(substitute(a, Dict(t => tv, x => xv, y => yv, z => zv); fold = Val(true))) for a in args]
            rhs = Symbolics.value(substitute(bc.rhs, merge(pmap, Dict(t => coords[1], x => coords[2], y => coords[3], z => coords[4])); fold = Val(true)))
            @test abs(rhs - u(coords...)) < 1.0e-12
        end
    end
    # The initial plane and each of the six faces appear exactly once.
    t0 = DomainSets.infimum(doms["t"])
    expected = vcat([[(1, Float64(t0))]], [[(i + 1, Float64(e))] for i in 1:3 for e in ends[i]])
    @test Set(signatures) == Set(expected)
end
