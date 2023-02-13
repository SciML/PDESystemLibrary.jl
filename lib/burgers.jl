"""
# The Inviscid Burgers Equation in 1D

The Inviscid Burgers equation is a model for the evolution of a fluid.
The fluid is assumed to be incompressible and inviscid, meaning that the fluid is not viscous and does not change in volume.
The fluid is also assumed to be one-dimensional, meaning that the fluid is only moving in one axis.
"""
function inviscid_burgers_monotonic()
    @parameters x t
    @variables u(..)
    Dx = Differential(x)
    Dt = Differential(t)
    x_min = 0.0
    x_max = 1.0
    t_min = 0.0
    t_max = 6.0

    analytic = u(t, x) ~ x / (t + 1)

    analytic_u(t, x) = x / (t + 1)

    eq = Dt(u(t, x)) ~ -u(t, x) * Dx(u(t, x))

    bcs = [u(0, x) ~ x,
        u(t, x_min) ~ analytic_u(t, x_min),
        u(t, x_max) ~ analytic_u(t, x_max)]

    domains = [t ∈ Interval(t_min, t_max),
        x ∈ Interval(x_min, x_max)]

    dx = 0.05

    tags = ["1D", "Monotonic", "Inviscid", "Burgers", "Advection", "Dirichlet"]

    @named inviscid_burgers_monotonic = PDESystem(eq, bcs, domains, [t, x], [u(t, x)];
                                                  analytic = analytic, metadata = tags)

    inviscid_burgers_monotonic
end

"""
# The Burgers Equation in 2D

The Burgers equation is a model for the evolution of a fluid.
This time the model has a viscosity term, which means that the fluid is viscous. The fluid is also assumed to be two-dimensional.
"""
function burgers_2d()
    @parameters x y t
    @variables u(..) v(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)
    Dxx = Differential(x)^2
    Dyy = Differential(y)^2

    R = 80.0

    x_min = y_min = t_min = 0.0
    x_max = y_max = 1.0
    t_max = 1.0

    #Exact solutions from: https://www.sciencedirect.com/science/article/pii/S0898122110003883

    u_exact(x, y, t) = 3 / 4 - 1 / (4 * (1 + exp(R * (-t - 4x + 4y) / 32)))
    v_exact(x, y, t) = 3 / 4 + 1 / (4 * (1 + exp(R * (-t - 4x + 4y) / 32)))
    analytic = [u(x, y, t) ~ u_exact(x, y, t),
        v(x, y, t) ~ v_exact(x, y, t)]

    eq = [
        Dt(u(x, y, t)) + u(x, y, t) * Dx(u(x, y, t)) + v(x, y, t) * Dy(u(x, y, t)) ~ (1 / R) *
                                                                                     (Dxx(u(x,
                                                                                            y,
                                                                                            t)) +
                                                                                      Dyy(u(x,
                                                                                            y,
                                                                                            t))),
        Dt(v(x, y, t)) + u(x, y, t) * Dx(v(x, y, t)) + v(x, y, t) * Dy(v(x, y, t)) ~ (1 / R) *
                                                                                     (Dxx(v(x,
                                                                                            y,
                                                                                            t)) +
                                                                                      Dyy(v(x,
                                                                                            y,
                                                                                            t))),
    ]

    domains = [x ∈ Interval(x_min, x_max),
        y ∈ Interval(y_min, y_max),
        t ∈ Interval(t_min, t_max)]

    bcs = [u(x, y, 0) ~ u_exact(x, y, 0),
        u(0, y, t) ~ u_exact(0, y, t),
        u(x, 0, t) ~ u_exact(x, 0, t),
        u(1, y, t) ~ u_exact(1, y, t),
        u(x, 1, t) ~ u_exact(x, 1, t), v(x, y, 0) ~ v_exact(x, y, 0),
        v(0, y, t) ~ v_exact(0, y, t),
        v(x, 0, t) ~ v_exact(x, 0, t),
        v(1, y, t) ~ v_exact(1, y, t),
        v(x, 1, t) ~ v_exact(x, 1, t)]

    tags = ["2D", "Non-Monotonic", "Viscous", "Burgers", "Advection", "Dirichlet"]

    @named burgers_2d = PDESystem(eq, bcs, domains, [t, x, y], [u(x, y, t), v(x, y, t)],
                                  analytic = analytic, metadata = tags)

    burgers_2d
end

all_systems = vcat(all_systems, [inviscid_burgers_monotonic(), burgers_2d()])
