# Linear Convection Systems
"""
# Linear Convection Equation

1D linear convection equation with periodic boundary conditions.

This equation models the flow of a substance, where the substance is moving at a constant velocity.

It can take any initial condition profile, as long as it is periodic.

```math
\\frac{\\partial u}{\\partial t} + c \\frac{\\partial u}{\\partial x} = 0
```
"""
function linear_convection(f, ps, name = :linear_convection)
    @variables u(..)
    @parameters t x v
    Dt = Differential(t)
    Dx = Differential(x)

    # 1D PDE and boundary conditions
    eq = Dt(u(t, x)) ~ -v * Dx(u(t, x))

    bcs = [u(0, x) ~ f(x),
        u(t, 0) ~ u(t, 1)]
    # Space and time domains
    domains = [t ∈ Interval(0.0, 10.0),
        x ∈ Interval(0.0, 1.0)]

    # Analytic solution
    u_exact = [u(t, x) ~ f(x - v * t)]

    tags = ["1D", "Periodic", "Linear", "Advection"]

    # PDE system
    lin_conv = PDESystem(eq, bcs, domains, [t, x], [u(t, x)], [v => ps[1]],
                         analytic = u_exact, metadata = tags, name = name)

    lin_conv
end

# sinusoidal input
convsin = linear_convection(x -> sin(2pi * x), [1], :convsin)
push!(convsin.metadata, "Sinusoidal")
push!(all_systems, convsin)

#cosine input
convcos = linear_convection(x -> cos(2pi * x), [-2], :convcos)
push!(convcos.metadata, "Sinusoidal")
push!(all_systems, convcos)

# triangular input
convtri = linear_convection(x -> 1.0 - abs(x - floor(x + 0.5)), [0.6], :convtri)
push!(convtri.metadata, "Triangular")
push!(all_systems, convtri)

# square wave
f = (x) -> IfElse.ifelse(x - floor(x) < 0.5, 1.0, -1.0)

convsquare = linear_convection(f, [-1.1], :convsquare)
push!(convsquare.metadata, "Square")
push!(convsquare.metadata, "Discontinuous")
push!(all_systems, convsquare)

convsquare = linear_convection(f, [2.1], :convsquare2)
push!(convsquare.metadata, "Square")
push!(convsquare.metadata, "Discontinuous")
push!(all_systems, convsquare)

"""
# Linear Convection Equation with Dirichlet Boundary Conditions 1

1D linear convection equation with Dirichlet boundary conditions.

This equation models the flow of a substance, where the substance is moving at a constant velocity.

v must be positive.

It can take any initial condition profile, as long as it is windowed over the domain.

```math
\\frac{\\partial u}{\\partial t} + v \\frac{\\partial u}{\\partial x} = 0
```

"""
function linear_convection_dirichlet1(f, ps, name = :linear_convection)
    @variables u(..)
    @parameters t x v
    Dt = Differential(t)
    Dx = Differential(x)

    # 1D PDE and boundary conditions
    eq = Dt(u(t, x)) ~ -v * Dx(u(t, x))

    bcs = [u(0, x) ~ f(x),
        u(t, 0) ~ 0.0]
    # Space and time domains
    domains = [t ∈ Interval(0.0, 2 / ps[1]),
        x ∈ Interval(0.0, 1.0)]

    # Analytic solution
    u_exact = [u(t, x) ~ windowlower(f(x - v * t), x - v * t)]

    tags = ["1D", "Dirichlet", "Linear", "Advection"]

    # PDE system
    lin_conv = PDESystem(eq, bcs, domains, [t, x], [u(t, x)], [v => ps[1]],
                         analytic = u_exact, metadata = tags, name = name)

    lin_conv
end

function windowlower(f, x, domain = (0.0, 1.0))
    IfElse.ifelse(x <= domain[1], IfElse.ifelse(x > domain[2], f, 0.0), 0.0)
end

# sinusoidal input
convsin = linear_convection_dirichlet1(x -> sin(2pi * x), [1], :dconvsin)
push!(convsin.metadata, "Sinusoidal")
push!(all_systems, convsin)

convsin = linear_convection_dirichlet1(x -> sin(2pi * x), [0.5], :dconvsin2)
push!(convsin.metadata, "Sinusoidal")
push!(all_systems, convsin)

#cosine input
convcos = linear_convection_dirichlet1(x -> cos(2pi * x), [2], :dconvcos)
push!(convcos.metadata, "Sinusoidal")
push!(all_systems, convcos)

# triangular input
convtri = linear_convection_dirichlet1(x -> 1.0 - abs(x - floor(x + 0.5)), [0.6], :dconvtri)
push!(convtri.metadata, "Triangular")
push!(all_systems, convtri)

# square wave
f = (x) -> IfElse.ifelse(x - floor(x) < 0.5, 1.0, -1.0)

convsquare = linear_convection_dirichlet1(f, [1.1], :dconvsquare)
push!(convsquare.metadata, "Square")
push!(convsquare.metadata, "Discontinuous")
push!(all_systems, convsquare)

"""
# Linear Convection Equation with Dirichlet Boundary Conditions 2

1D linear convection equation with Dirichlet boundary conditions.

This equation models the flow of a substance, where the substance is moving at a constant velocity.

v must be positive.

It can take any initial condition profile, as long as it is windowed over the domain.

```math
\\frac{\\partial u}{\\partial t} - v \\frac{\\partial u}{\\partial x} = 0
```

"""
function linear_convection_dirichlet2(f, ps, name = :linear_convection)
    @variables u(..)
    @parameters t x v
    Dt = Differential(t)
    Dx = Differential(x)

    # 1D PDE and boundary conditions
    eq = Dt(u(t, x)) ~ v * Dx(u(t, x))

    bcs = [u(0, x) ~ f(x),
        u(t, 2 / ps[1]) ~ 0.0]
    # Space and time domains
    domains = [t ∈ Interval(0.0, 2 / ps[1]),
        x ∈ Interval(0.0, 1.0)]

    # Analytic solution
    u_exact = [u(t, x) ~ windowlower(f(x - v * t), x - v * t)]

    tags = ["1D", "Dirichlet", "Linear", "Advection"]

    # PDE system
    lin_conv = PDESystem(eq, bcs, domains, [t, x], [u(t, x)], [v => ps[1]],
                         analytic = u_exact, metadata = tags, name = name, eval_module = @__MODULE__)

    lin_conv
end

function windowupper(f, x, domain = (0.0, 1.0))
    IfElse.ifelse(x < domain[1], IfElse.ifelse(x >= domain[2], f, 0.0), 0.0)
end

# sinusoidal input
convsin = linear_convection_dirichlet1(x -> sin(2pi * x), [1], :ddconvsin)
push!(convsin.metadata, "Sinusoidal")
push!(all_systems, convsin)

convsin = linear_convection_dirichlet1(x -> sin(2pi * x), [0.5], :ddconvsin2)
push!(convsin.metadata, "Sinusoidal")
push!(all_systems, convsin)

#cosine input
convcos = linear_convection_dirichlet1(x -> cos(2pi * x), [0.7], :ddconvcos)
push!(convcos.metadata, "Sinusoidal")
push!(all_systems, convcos)

# triangular input
convtri = linear_convection_dirichlet1(x -> 1.0 - abs(x - floor(x + 0.5)), [0.1],
                                       :ddconvtri)
push!(convtri.metadata, "Triangular")
push!(all_systems, convtri)

# square wave
sq = f = (x) -> IfElse.ifelse(x - floor(x) < 0.5, 1.0, -1.0)

convsquare = linear_convection_dirichlet1(f, [3.1], :ddconvsquare)
push!(convsquare.metadata, "Square")
push!(convsquare.metadata, "Discontinuous")
push!(all_systems, convsquare)

"""
# Linear Convection Equation with Dirichlet Boundary Conditions 3

1D linear convection equation with Dirichlet boundary conditions. the left boundary is time dependent.

This equation models the flow of a substance, where the substance is moving at a constant velocity.

v must be positive.

It can take any initial condition profile, as long as it is windowed over the domain.

```math
\\frac{\\partial u}{\\partial t} + v \\frac{\\partial u}{\\partial x} = 0
```

"""
function linear_convection_dirichlet3(f, h, ps, name = :linear_convection)
    @variables u(..)
    @parameters t x v
    Dt = Differential(t)
    Dx = Differential(x)

    # 1D PDE and boundary conditions
    eq = Dt(u(t, x)) ~ -v * Dx(u(t, x))

    bcs = [u(0, x) ~ f(x),
        u(t, 0.0) ~ h(t)]
    # Space and time domains
    domains = [t ∈ Interval(0.0, 10.0),
        x ∈ Interval(0.0, 1.0)]

    # Analytic solution
    u_exact = [u(t, x) ~ IfElse.ifelse(x > v * t, f(x - v * t), h(t - x / v))]

    tags = ["1D", "Dirichlet", "Linear", "Advection"]

    # PDE system
    lin_conv = PDESystem(eq, bcs, domains, [t, x], [u(t, x)], [v => ps[1]],
                         analytic = u_exact, metadata = tags, name = name, eval_module = @__MODULE__)

    lin_conv
end

funcs = [x -> x, x -> x^2, x -> x^3, sinpi, cospi, x -> 1.0 - abs(x - floor(x + 0.5)), sq]
function add_systems!(all_systems, funcs, sysconstructor, name)
    i = 0
    for f in funcs
        for h in funcs
            conv = sysconstructor(f, h, [rand()], Symbol(name, i))
            push!(all_systems, conv)
            i += 1
            conv = sysconstructor(f, x -> -2 * h(x), [2 * rand()], Symbol(name, i))
            push!(all_systems, conv)
            i += 1
            conv = sysconstructor(x -> -10 * f(x), h, [rand()], Symbol(name, i))
            push!(all_systems, conv)
            i += 1
            conv = sysconstructor(x -> -6 * f(x), x -> -5 * h(x), [rand()],
                                  Symbol(name, i))
            push!(all_systems, conv)
            i += 1
        end
    end
end

add_systems!(all_systems, funcs, linear_convection_dirichlet3, "funcconv")

"""
# Linear Convection Equation with Dirichlet Boundary Conditions 4

1D linear convection equation with Dirichlet boundary conditions. the left boundary is time dependent.

This equation models the flow of a substance, where the substance is moving at a constant velocity.

v must be positive.

It can take any initial condition profile, as long as it is windowed over the domain.

```math
\\frac{\\partial u}{\\partial t} - v \\frac{\\partial u}{\\partial x} = 0
```

"""
function linear_convection_dirichlet4(f, h, ps, name = :linear_convection)
    @variables u(..)
    @parameters t x v
    Dt = Differential(t)
    Dx = Differential(x)

    # 1D PDE and boundary conditions
    eq = Dt(u(t, x)) ~ v * Dx(u(t, x))

    bcs = [u(0, x) ~ f(x),
        u(t, 1.0) ~ h(t)]
    # Space and time domains
    domains = [t ∈ Interval(0.0, 10.0),
        x ∈ Interval(0.0, 1.0)]

    # Analytic solution
    u_exact = [u(t, x) ~ IfElse.ifelse(x < v * t, f(x + v * t), h(t + x / v))]

    tags = ["1D", "Dirichlet", "Linear", "Advection"]

    # PDE system
    lin_conv = PDESystem(eq, bcs, domains, [t, x], [u(t, x)], [v => ps[1]],
                         analytic = u_exact, metadata = tags, name = name)

    lin_conv
end

add_systems!(all_systems, funcs, linear_convection_dirichlet4, "funcconvneg")

"""
# Convection Diffusion Equation in 1D

1D convection diffusion equation with Dirichlet boundary conditions.

This equation models the flow of a substance, where the substance is moving at a constant velocity and diffusing in a medium.

It is initialized with a zero profile, and the boundary conditions are set to 1 at the right boundary and 0 at the left boundary.

```math
\\frac{\\partial u}{\\partial t} = \\kappa \\frac{\\partial^2 u}{\\partial z^2} + v \\frac{\\partial u}{\\partial z} +
```

Reference Solution found here: https://www.12000.org/my_notes/diffusion_convection_PDE/diffusion_convection_PDE.htm
Author of reference solution: Nasser M. Abbasi

"""
function convection_diffusion(L, ps, name = :convection_diffusion)
    @variables f(..)
    @parameters z, t, k, v
    Dt = Differential(t)
    Dz = Differential(z)
    Dzz = Differential(z)^2

    # 1D PDE and boundary conditions
    eq = Dt(f(t, z)) ~ k * Dzz(f(t, z)) + v * Dz(f(t, z))

    f_0(z) = IfElse.ifelse(z == L, 1.0, 0.0)

    bcs = [f(0, z) ~ f_0(z),
        f(t, 0) ~ 0.0, f(t, L) ~ 1.0]
    # Space and time domains
    domains = [t ∈ Interval(0.0, 30.0),
        z ∈ Interval(0.0, L)]

    # Analytic/reference solution
    λ(n) = (n * π / L)^2

    maxiters = 1000

    function A(ps, t, z)
        k = ps[1]
        v = ps[2]
        exp(-((t * v^2) / (4k) + (v * z) / (2k)))
    end

    function u(ps, t, z)
        k = ps[1]
        v = ps[2]
        a = z / L * exp((t * v^2) / (4k) + (v * L) / (2k))
        s = mapreduce((+), 1:maxiters) do n
            acu = (2 * (-1)^n * v^2 * exp(-k * λ(n) * t + (v * L) / (2k)) *
                   (exp(k * λ(n) * t + (t * v^2) / (4k)) - 1)) /
                  (n * π * (4 * k^2 * λ(n) + v^2))
            acu += 2 / L * ((-1)^n) * exp(v * L / (2k)) * exp(-k * λ(n) * t) / sqrt(λ(n))
            acu *= sin(λ(n) * z)

            acu
        end
        a + s
    end

    ref = [f(t, z) => (ps, t, z) -> A(ps, t, z) * u(ps, t, z)]

    tags = ["1D", "Dirichlet", "Linear", "Advection", "Diffusion", "Monotonic"]

    # PDE system
    convdiff = PDESystem(eq, bcs, domains, [t, z], [f(t, z)],
                         [k => ps[1], v => ps[2]], analytic_func = ref, metadata = tags,
                         name = name, eval_module = @__MODULE__)

    convdiff
end

# L = 1.0, k = 0.01, v = 1.0
convdiff1 = convection_diffusion(1.0, [0.01, 1.0], :convdiff1)
push!(all_systems, convdiff1)

# L = 5.0, k = 1.0, v = -1.0
convdiff2 = convection_diffusion(5.0, [1.0, -1.0], :convdiff2)
push!(all_systems, convdiff2)

# L = 10.0, k = 0.01, v = -3.0
convdiff3 = convection_diffusion(10.0, [0.01, -3.0], :convdiff3)
push!(all_systems, convdiff3)

# L = 10.0, k = 0.5, v = 3.0
convdiff4 = convection_diffusion(10.0, [0.5, 3.0], :convdiff4)
push!(all_systems, convdiff4)

"""
# Transport Equation in 1D

1D transport equation without boundary conditions

This equation models the flow of a substance, where the substance is moving at a constant velocity.

it is initialized with a sinusoid, and has a sinusoidal source term.

```math
\\frac{\\partial u}{\\partial t} + v 2\\frac{\\partial u}{\\partial z} = \\sin(z)
```
"""
function trans_sin()
    @variables u(..)
    @parameters z, t
    Dt = Differential(t)
    Dz = Differential(z)

    # 1D PDE and boundary conditions
    eq = Dt(u(t, z)) + 2 * Dz(u(t, z)) ~ sin(z)

    u_exact(t, z) = sin(z - 2t) + 0.5 * cos(z - 2t) - 0.5 * cos(z)

    bcs = [u(0, z) ~ u_exact(0, z),
        u(t, 0) ~ u_exact(t, 0),
        (t, 2π) ~ u_exact(t, 2π)]

    # Space and time domains

    domains = [t ∈ Interval(0.0, 1.0),
        z ∈ Interval(0.0, 2π)]

    # Analytic/reference solution
    ref = [u(t, z) ~ u_exact(t, z)]

    tags = ["1D", "Transport", "Linear", "Sinusoidal", "Inhomogeneous", "Advection"]

    # PDESystem

    @named trans_sin = PDESystem(eqs, bcs, domains, [t, z], [u(t, z)], analytic = ref,
                                 metadata = tags, eval_module = @__MODULE__)

    trans_sin
end
