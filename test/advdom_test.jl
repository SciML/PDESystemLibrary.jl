using PDESystemLibrary
using Test
using DomainSets
using Symbolics

@testset "Advection-dominated Gaussian" begin
    system = only(filter(s -> s.name == :advdom_gaussian, get_pdesys_with_tags(["AdvectionDominated"])))
    @test all(t -> t in system.metadata, ["1D", "Advection", "Diffusion", "AdvectionDominated"])

    pvals = Dict(string(first(p)) => last(p) for p in system.ps)
    v, ε = pvals["v"], pvals["eps"]
    xdom = only(filter(d -> string(d.variables) == "x", system.domain)).domain
    tdom = only(filter(d -> string(d.variables) == "t", system.domain)).domain
    L = DomainSets.supremum(xdom) - DomainSets.infimum(xdom)
    @test v * L / ε == 200  # documented Peclet number

    analytic = only(values(system.analytic_func))
    ps = [v, ε]
    ts = range(DomainSets.infimum(tdom), DomainSets.supremum(tdom); length = 11)
    xs = range(0.0, 1.0; length = 201)
    @test all(isfinite, [analytic(ps, t, x) for t in ts, x in xs])

    # PDE residual through central differences on the exact solution.
    h = 1.0e-4
    for t in ts[2:(end - 1)], x in xs[10:10:(end - 10)]
        ut = (analytic(ps, t + h, x) - analytic(ps, t - h, x)) / 2h
        ux = (analytic(ps, t, x + h) - analytic(ps, t, x - h)) / 2h
        uxx = (analytic(ps, t, x + h) - 2analytic(ps, t, x) + analytic(ps, t, x - h)) / h^2
        @test abs(ut + v * ux - ε * uxx) < 1.0e-4
    end

    # Initial and boundary conditions agree with the exact solution; the pulse starts at x0 = 0.2.
    @test analytic(ps, 0.0, 0.2) ≈ 1.0
    @test analytic(ps, 0.0, 0.0) < 1.0e-6
    @test analytic(ps, 0.5, 1.0) < 1.0e-2
    @test analytic(ps, 0.5, 0.7) ≈ sqrt(0.1 / 0.6)

    # The stored equation and boundary conditions are satisfied by the exact solution.
    function stored_residuals(sys)
        t, x = sys.ivs
        pmap = Dict(first(p) => last(p) for p in sys.ps)
        u_sym = only(sys.dvs)
        exact = only(sys.analytic).rhs
        pde = only(sys.eqs)
        Dt, Dx = Differential(t), Differential(x)
        derivs = Dict(
            d(u_sym) => Symbolics.expand_derivatives(d(exact))
                for d in (Dt, Dx, Dx^2)
        )
        pde_res = substitute(pde.lhs - pde.rhs, derivs)
        pde_vals = [
            Symbolics.value(substitute(pde_res, merge(pmap, Dict(t => tv, x => xv)); fold = Val(true)))
                for tv in ts[2:end], xv in xs[10:10:(end - 10)]
        ]
        bc_vals = Float64[]
        bc_points = Tuple[]
        for bc in sys.bcs
            args = Symbolics.arguments(Symbolics.unwrap(bc.lhs))
            # A coordinate is fixed (a number) or free (still a variable); free ones are swept below.
            fixed = [isempty(Symbolics.get_variables(a)) ? Symbolics.value(a) : nothing for a in args]
            push!(bc_points, Tuple(fixed))
            for s in range(0.0, 0.5; length = 5)  # each boundary coordinate is either fixed or swept
                tv, xv = [Symbolics.value(substitute(a, Dict(t => s, x => s); fold = Val(true))) for a in args]
                rhs = Symbolics.value(substitute(bc.rhs, merge(pmap, Dict(t => tv, x => xv)); fold = Val(true)))
                push!(bc_vals, rhs - analytic(ps, tv, xv))
            end
        end
        return pde_vals, bc_vals, bc_points
    end
    pde_vals, bc_vals, bc_points = stored_residuals(system)
    @test maximum(abs, pde_vals) < 1.0e-8
    @test maximum(abs, bc_vals) < 1.0e-12
    # Each boundary condition sits on the initial time or one end of the x domain.
    @test all(bc_points) do (tv, xv)
        tv == DomainSets.infimum(tdom) || xv in (DomainSets.infimum(xdom), DomainSets.supremum(xdom))
    end
end
