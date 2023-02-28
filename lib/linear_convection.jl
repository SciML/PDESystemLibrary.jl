"""
# Linear Convection Equation

1D linear convection equation with periodic boundary conditions.

This equation models the flow of a substance, where the substance is moving at a constant velocity.

It can take any initial condition profile, as long as it is periodic.

```math
\\frac{\\partial u}{\\partial t} + c \\frac{\\partial u}{\\partial x} = 0
```
"""
function linear_convection(f, name = :linear_convection)
    @variables x t u(..)
    Dt = Differential(t)
    Dx = Differential(x)

    # 1D PDE and boundary conditions
    eq = Dt(u(t, x)) ~ - Dx(u(t, x))

    @register_symbolic f(x)

    bcs = [u(0, x) ~ f(x),
           u(t, 0) ~ u(t, 1)]
    # Space and time domains
    domains = [t ∈ Interval(0.0, 10.0),
               x ∈ Interval(0.0, 1.0)]

    # Analytic solution
    u_exact = [u(t, x) ~ f(x - t)]

    tags = ["1D", "Periodic", "Linear", "Advection"]

    # PDE system
    lin_conv = PDESystem(eq, bcs, domains, [t, x], [u(t, x)], analytic = u_exact, metadata = tags, name = name)

    lin_conv
end

# sinusoidal input
convsin = linear_convection(x -> sin(2pi * x), :convsin)
push!(convsin.metadata, "Sinusoidal")
push!(all_systems, convsin)

#cosine input
convcos = linear_convection(x -> cos(2pi * x), :convcos)
push!(convcos.metadata, "Sinusoidal")
push!(all_systems, convcos)

# triangular input
convtri = linear_convection(x -> 1.0 - abs(x - floor(x + 0.5)), :convtri)
push!(convtri.metadata, "Triangular")
push!(all_systems, convtri)

# square wave
f = (x) -> ifelse(x - floor(x) < 0.5, 1.0, -1.0)

convsquare = linear_convection(f, :convsquare)
push!(convsquare.metadata, "Square")
push!(convsquare.metadata, "Discontinuous")
push!(all_systems, convsquare)

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
    @variables x t f(..)
    @parameters k, v
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2

    # 1D PDE and boundary conditions
    eq = Dt(f(t, x)) ~ k * Dxx(f(t, x)) + v * Dx(f(t, x))

    f_0(x) = x == L ? 1.0 : 0.0
    @register_symbolic f_0(x)

    bcs = [f(0, x) ~ f_0(x),
           f(t, 0) ~ 0.0, f(t, L) ~ 1.0]
    # Space and time domains
    domains = [t ∈ Interval(0.0, 30.0),
               x ∈ Interval(0.0, L)]

    # Analytic/reference solution
    λ(n) = (n*π/L)^2

    maxiters = 1000

    function A(ps, t, z)
        k = ps[1]
        v = ps[2]
        exp(-((t*v^2)/(4k) + (v*z)/(2k)))
    end

    function u(ps, t, z)
        k = ps[1]
        v = ps[2]
        a = z/L * exp((t*v^2)/(4k) + (v*L)/(2k))
        s = mapreduce((+), 1:maxiters) do n
            acu = (2 * (-1)^n *v^2 * exp(-k*λ(n)*t + (v*L)/(2k))*(exp(k*λ(n)*t + (t*v^2)/(4k)) - 1)) / (n*π * (4*k^2*λ(n) + v^2))
            acu += 2/L*((-1)^n)*exp(v*L/(2k))*exp(-k*λ(n)*t)/sqrt(λ(n))
            acu *= sin(λ(n)*z)

            acu
        end
        a + s
    end

    ref = [f(t, z) => (ps, t, z) -> A(ps, t, z) * u(ps, t, z)]

    tags = ["1D", "Dirichlet", "Advection", "Diffusion", "Monotonic"]

    # PDE system
    convdiff = PDESystem(eq, bcs, domains, [t, x], [f(t, x)], [k => ps[1], v => ps[2]] analytic = ref, metadata = tags, name = name)

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
