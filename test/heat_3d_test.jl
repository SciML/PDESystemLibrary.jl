using PDESystemLibrary
using Test
using DomainSets
using Symbolics

@testset "heat_3d analytic solution" begin
    system = only(get_pdesys_with_tags(["3D", "Heat", "Dirichlet"]))
    @test system.name == :heat_3d
    kval = only(system.ps)[2]
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
    @test u(t_end, 0.5, 0.5, 0.5) > 0.1

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

    # Each boundary condition sits on the initial time or one end of an x, y or z interval,
    # and its right-hand side matches the exact solution there.
    doms = Dict(string(d.variables) => d.domain for d in system.domain)
    ends = [(DomainSets.infimum(doms[v]), DomainSets.supremum(doms[v])) for v in ("x", "y", "z")]
    @test length(system.bcs) == 7
    for bc in system.bcs
        args = Symbolics.arguments(Symbolics.unwrap(bc.lhs))
        fixed = [isempty(Symbolics.get_variables(a)) ? Symbolics.value(a) : nothing for a in args]
        @test fixed[1] == DomainSets.infimum(doms["t"]) ||
            any(i -> fixed[i + 1] in ends[i], 1:3)
        for (tv, xv, yv, zv) in samples
            coords = [Symbolics.value(substitute(a, Dict(t => tv, x => xv, y => yv, z => zv); fold = Val(true))) for a in args]
            rhs = Symbolics.value(substitute(bc.rhs, merge(pmap, Dict(t => coords[1], x => coords[2], y => coords[3], z => coords[4])); fold = Val(true)))
            @test abs(rhs - u(coords...)) < 1.0e-12
        end
    end
end
