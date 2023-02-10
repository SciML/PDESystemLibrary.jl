begin
    @variables x t
    @parameters D

    Dxx = Differential(x)
    Dt = Differential(t)

    eqs = [Dt(u(t, x)) ~ D * Dxx(u(t, x))]
    bcs = [u(0, x) ~ sin(2pi * x), u(t, 0) ~ 0.0, u(t, 1) ~ 0.0]
    domains = [t ∈ IntervalDomain(0.0, 1.0), x ∈ IntervalDomain(0.0, 1.0)]

    function u_exact(t, x; ps = (1.0,))
        exp(-4 * (pi^2) * ps[D] * t) * sin(2pi * x)
    end

    analytic = [u(t, x) => u_exact]

    @named heat_1d = PDESystem(eqs, bcs, domains, [t, x], [u(t, x)],
                               parameters = [D => 1.0])

    heat_1d
end

push!(all_systems, heat_1d)
push!(linear_diffusion_systems, heat_1d)
