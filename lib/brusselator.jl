"""
# The Brusselator Equation in 2D

A 2D version of the Brusselator equation, which is a model for the formation of a pattern in a chemical reaction.
The reaction is between two chemicals, A and B, which are produced and consumed by the reaction.
The reaction is autocatalytic, meaning that the production of A and B is dependent on the concentration of A and B.
The reaction is also bistable, meaning that there are two stable states for the system.
The system is initialized with a small perturbation in the concentration of A, which causes the system to evolve.
The system is then perturbed again, eventually forming a bistable state.
This is a classic example of a pattern forming chemical reaction.

The Brusselator equation is given by:

```math
\\frac{\\partial u}{\\partial t} = 1 + v u^2 - 4.4 u + \\alpha \\Delta u + f(x, y, t)
```

```math
\\frac{\\partial v}{\\partial t} = 3.4 - vu^2 + u^2 + \\alpha \\Delta v
```

where ``\\alpha`` is a parameter that controls the diffusion of the system, and ``f(x, y, t)`` is a forcing term.

The initial conditions are given by:

```math
u(x, y, 0) = 22(y (1 - y))^{\\frac{3}{2}}
```

```math
v(x, y, 0) = 27(x (1 - x))^{\\frac{3}{2}}
```

The boundary conditions are periodic in both ``x`` and ``y``.
"""
bruss = begin
    @parameters x y t
    @variables u(..) v(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)
    Dxx = Differential(x)^2
    Dyy = Differential(y)^2

    ∇²(u) = Dxx(u) + Dyy(u)

    brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0

    x_min = y_min = t_min = 0.0
    x_max = y_max = 1.0
    t_max = 11.5

    α = 10.0

    u0(x, y, t) = 22(y * (1 - y))^(3 / 2)
    v0(x, y, t) = 27(x * (1 - x))^(3 / 2)

    eq = [
        Dt(u(x, y, t)) ~ 1.0 + v(x, y, t) * u(x, y, t)^2 - 4.4 * u(x, y, t) +
                         α * ∇²(u(x, y, t)) + brusselator_f(x, y, t),
        Dt(v(x, y, t)) ~ 3.4 * u(x, y, t) - v(x, y, t) * u(x, y, t)^2 + α * ∇²(v(x, y, t))]

    domains = [x ∈ Interval(x_min, x_max),
        y ∈ Interval(y_min, y_max),
        t ∈ Interval(t_min, t_max)]

    bcs = [u(x, y, 0) ~ u0(x, y, 0),
        u(0, y, t) ~ u(1, y, t),
        u(x, 0, t) ~ u(x, 1, t), v(x, y, 0) ~ v0(x, y, 0),
        v(0, y, t) ~ v(1, y, t),
        v(x, 0, t) ~ v(x, 1, t)]

    # Solve reference problem
    xyd_brusselator = range(0, stop = 1, length = N)
    brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
    limit(a, N) = a == N + 1 ? 1 : a == 0 ? N : a
    function brusselator_2d_loop(du, u, p, t)
        A, B, alpha, dx = p
        alpha = alpha / dx^2
        @inbounds for I in CartesianIndices((N, N))
            i, j = Tuple(I)
            x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
            ip1, im1, jp1, jm1 = limit(i + 1, N), limit(i - 1, N), limit(j + 1, N),
                                 limit(j - 1, N)
            du[i, j, 1] = alpha *
                          (u[im1, j, 1] + u[ip1, j, 1] + u[i, jp1, 1] + u[i, jm1, 1] -
                           4u[i, j, 1]) +
                          B + u[i, j, 1]^2 * u[i, j, 2] - (A + 1) * u[i, j, 1] +
                          brusselator_f(x, y, t)
            du[i, j, 2] = alpha *
                          (u[im1, j, 2] + u[ip1, j, 2] + u[i, jp1, 2] + u[i, jm1, 2] -
                           4u[i, j, 2]) +
                          A * u[i, j, 1] - u[i, j, 1]^2 * u[i, j, 2]
        end
    end
    p = (3.4, 1.0, 10.0, step(xyd_brusselator))

    function init_brusselator_2d(xyd)
        N = length(xyd)
        u = zeros(N, N, 2)
        for I in CartesianIndices((N, N))
            x = xyd[I[1]]
            y = xyd[I[2]]
            u[I, 1] = 22 * (y * (1 - y))^(3 / 2)
            u[I, 2] = 27 * (x * (1 - x))^(3 / 2)
        end
        u
    end
    u0_manual = init_brusselator_2d(xyd_brusselator)
    prob = ODEProblem(brusselator_2d_loop, u0_manual, (0.0, 11.5), p)

    msol = solve(prob, TRBDF2(), saveat = 0.01) # 2.771 s (5452 allocations: 65.73 MiB)

    # Create Interpolation for reference solution
    u_ref = zeros(N, N, length(msol.t))
    v_ref = zeros(N, N, length(msol.t))
    for (i, u) in enumerate(msol.u)
        u_ref[:, :, i] = u[:, :, 1]
        v_ref[:, :, i] = u[:, :, 2]
    end

    nodes = (xyd_brusselator, xyd_brusselator, msol.t)

    u_itp = interpolate(nodes, u_ref, Gridded(Linear()))
    v_itp = interpolate(nodes, v_ref, Gridded(Linear()))

    # periodic extrapolation

    u_extp = extrapolate(u_itp, Periodic())
    v_extp = extrapolate(v_itp, Periodic())

    u_func = (ps, x, y, t) -> u_extp(x, y, t)
    v_func = (ps, x, y, t) -> v_extp(x, y, t)

    analytic_func = [u(x, y, t) => u_func, v(x, y, t) => v_func]

    tags = ["2D", "Brusselator", "Periodic", "Nonlinear", "Reaction", "Diffusion"]

    @named bruss = PDESystem(eq, bcs, domains, [x, y, t], [u(x, y, t), v(x, y, t)],
                             analytic_func = analytic_func)

    bruss
end

push!(all_systems, bruss)
push!(nonlinear_systems, bruss)
