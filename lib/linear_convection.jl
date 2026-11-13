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

    bcs = [
        u(0, x) ~ f(x),
        u(t, 0) ~ u(t, 1),
    ]
    # Space and time domains
    domains = [
        t ∈ Interval(0.0, 10.0),
        x ∈ Interval(0.0, 1.0),
    ]

    # Analytic solution
    u_exact = [u(t, x) ~ f(x - v * t)]

    tags = ["1D", "Periodic", "Linear", "Advection"]

    # PDE system
    lin_conv = PDESystem(
        eq, bcs, domains, [t, x], [u(t, x)], [v => ps[1]],
        analytic = u_exact, analytic_func = analytic_funcs(u_exact, [v => ps[1]]), metadata = tags, name = name, eval_module = @__MODULE__
    )

    return lin_conv
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
f = (x) -> ifelse(x - floor(x) < 0.5, 1.0, -1.0)

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

    bcs = [
        u(0, x) ~ f(x),
        u(t, 0) ~ 0.0,
    ]
    # Space and time domains
    domains = [
        t ∈ Interval(0.0, 2 / ps[1]),
        x ∈ Interval(0.0, 1.0),
    ]

    # Analytic solution
    u_exact = [u(t, x) ~ windowlower(f(x - v * t), x - v * t)]

    tags = ["1D", "Dirichlet", "Linear", "Advection"]

    # PDE system
    lin_conv = PDESystem(
        eq, bcs, domains, [t, x], [u(t, x)], [v => ps[1]],
        analytic = u_exact, analytic_func = analytic_funcs(u_exact, [v => ps[1]]), metadata = tags, name = name, eval_module = @__MODULE__
    )

    return lin_conv
end

function windowlower(f, x, domain = (0.0, 1.0))
    return ifelse(x <= domain[1], ifelse(x > domain[2], f, 0.0), 0.0)
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
f = (x) -> ifelse(x - floor(x) < 0.5, 1.0, -1.0)

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

    bcs = [
        u(0, x) ~ f(x),
        u(t, 2 / ps[1]) ~ 0.0,
    ]
    # Space and time domains
    domains = [
        t ∈ Interval(0.0, 2 / ps[1]),
        x ∈ Interval(0.0, 1.0),
    ]

    # Analytic solution
    u_exact = [u(t, x) ~ windowlower(f(x - v * t), x - v * t)]

    tags = ["1D", "Dirichlet", "Linear", "Advection"]

    # PDE system
    lin_conv = PDESystem(
        eq, bcs, domains, [t, x], [u(t, x)], [v => ps[1]],
        analytic = u_exact, analytic_func = analytic_funcs(u_exact, [v => ps[1]]), metadata = tags, name = name, eval_module = @__MODULE__
    )

    return lin_conv
end

function windowupper(f, x, domain = (0.0, 1.0))
    return ifelse(x < domain[1], ifelse(x >= domain[2], f, 0.0), 0.0)
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
convtri = linear_convection_dirichlet1(
    x -> 1.0 - abs(x - floor(x + 0.5)), [0.1],
    :ddconvtri
)
push!(convtri.metadata, "Triangular")
push!(all_systems, convtri)

# square wave
sq = f = (x) -> ifelse(x - floor(x) < 0.5, 1.0, -1.0)

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

    bcs = [
        u(0, x) ~ f(x),
        u(t, 0.0) ~ h(t),
    ]
    # Space and time domains
    domains = [
        t ∈ Interval(0.0, 10.0),
        x ∈ Interval(0.0, 1.0),
    ]

    # Analytic solution
    u_exact = [u(t, x) ~ ifelse(x > v * t, f(x - v * t), h(t - x / v))]

    tags = ["1D", "Dirichlet", "Linear", "Advection"]

    # PDE system
    lin_conv = PDESystem(
        eq, bcs, domains, [t, x], [u(t, x)], [v => ps[1]],
        analytic = u_exact, analytic_func = analytic_funcs(u_exact, [v => ps[1]]), metadata = tags, name = name, eval_module = @__MODULE__
    )

    return lin_conv
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
            conv = sysconstructor(
                x -> -6 * f(x), x -> -5 * h(x), [rand()],
                Symbol(name, i)
            )
            push!(all_systems, conv)
            i += 1
        end
    end
    return
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

    bcs = [
        u(0, x) ~ f(x),
        u(t, 1.0) ~ h(t),
    ]
    # Space and time domains
    domains = [
        t ∈ Interval(0.0, 10.0),
        x ∈ Interval(0.0, 1.0),
    ]

    # Analytic solution
    u_exact = [u(t, x) ~ ifelse(x < v * t, f(x + v * t), h(t + x / v))]

    tags = ["1D", "Dirichlet", "Linear", "Advection"]

    # PDE system
    lin_conv = PDESystem(
        eq, bcs, domains, [t, x], [u(t, x)], [v => ps[1]],
        analytic = u_exact, analytic_func = analytic_funcs(u_exact, [v => ps[1]]), metadata = tags, name = name, eval_module = @__MODULE__
    )

    return lin_conv
end

add_systems!(all_systems, funcs, linear_convection_dirichlet4, "funcconvneg")

# Fixed-point BigInt arithmetic (value = n / 2^p) for the γL > 1 series, whose terms of size e^{γL}
# cancel to O(1). It uses only public BigInt operations, so it neither reads nor sets the global
# BigFloat precision, which is process-wide on Julia 1.10.
function fixed(x, p)
    m, e = frexp(Float64(x))
    n = BigInt(ldexp(m, 53))
    return p + e - 53 >= 0 ? n << (p + e - 53) : n >> -(p + e - 53)
end

function fixedexp(x, p)
    m = max(ndigits(abs(x); base = 2) - p + 16, 0)  # reduce below 2^-16, square back m times
    r = x >> m
    term = total = BigInt(1) << p
    j = 0
    while term != 0
        term = div((term * r) >> p, j += 1)
        total += term
    end
    for _ in 1:m
        total = (total * total) >> p
    end
    return total
end

function fixedsincos(x, p)
    sn = cs = BigInt(0)
    term = BigInt(1) << p
    j = 0
    while term != 0
        j % 4 == 0 ? (cs += term) : j % 4 == 1 ? (sn += term) : j % 4 == 2 ? (cs -= term) : (sn -= term)
        term = div((term * x) >> p, j += 1)
    end
    return sn, cs
end

function fixedpi(p)  # Machin: 16 atan(1/5) - 4 atan(1/239)
    function ataninv(n)
        x = (BigInt(1) << p) ÷ n
        total = x
        for k in 3:2:p
            x = x ÷ n^2
            x == 0 && break
            total += (k % 4 == 1 ? 1 : -1) * (x ÷ k)
        end
        return total
    end
    return 16 * ataninv(5) - 4 * ataninv(239)
end

# steady + transient of the convection_diffusion solution (see there) for γ > 0, in fixed point
# with 64 + ⌈3γL⌉ bits plus 64 guard bits; the series is summed with the recurrences
# e^{-k t μₙ²} = e^{-k t μ₁²n²} and sin(nθ) = 2cos(θ) sin((n-1)θ) - sin((n-2)θ).
function convdiff_fixed(L, k, γ, t, z, nterms)
    L, k, γ, t, z = Float64.((L, k, γ, t, z))
    p = 128 + ceil(Int, 3γ * L)
    one_ = BigInt(1) << p
    mul(a, b) = (a * b) >> p
    quo(a, b) = div(a << p, b)
    Lₚ, kₚ, γₚ, tₚ, zₚ = fixed.((L, k, γ, t, z), p)
    steady = quo(one_ - fixedexp(-2mul(γₚ, zₚ), p), one_ - fixedexp(-2mul(γₚ, Lₚ), p))
    μ₁ = quo(fixedpi(p), Lₚ)
    sinθ, cosθ = fixedsincos(mul(μ₁, zₚ), p)
    γ² = mul(γₚ, γₚ)
    kt = mul(kₚ, tₚ)
    decay = fixedexp(-mul(kt, mul(μ₁, μ₁)), p)  # e^{-k t μ₁²}
    q, ratio, step = one_, decay, mul(decay, decay)
    sₙ₋₁, sₙ = BigInt(0), sinθ
    series = BigInt(0)
    for n in 1:nterms
        q = mul(q, ratio)
        q == 0 && break
        ratio = mul(ratio, step)
        μ = n * μ₁
        term = mul(quo(2μ, mul(Lₚ, γ² + mul(μ, μ))), mul(q, sₙ))
        series += iseven(n) ? term : -term
        sₙ₋₁, sₙ = sₙ, mul(2cosθ, sₙ) - sₙ₋₁
    end
    transient = mul(mul(fixedexp(mul(γₚ, Lₚ - zₚ), p), fixedexp(-mul(kt, γ²), p)), series)
    total = steady + transient
    shift = max(ndigits(abs(total); base = 2) - 64, 0)
    return ldexp(Float64(total >> shift), shift - p)
end

"""
# Convection Diffusion Equation in 1D

1D convection diffusion equation with Dirichlet boundary conditions.

This equation models the flow of a substance, where the substance is moving at a constant velocity and diffusing in a medium.

It is initialized with a zero profile, and the boundary conditions are set to 1 at the right boundary and 0 at the left boundary.

```math
\\frac{\\partial u}{\\partial t} = \\kappa \\frac{\\partial^2 u}{\\partial z^2} + v \\frac{\\partial u}{\\partial z} +
```

The analytic solution sums at most 50000 sine-series terms and throws an `ArgumentError` naming the minimum supported `t` when more would be needed (about `t ≥ 1.6e-7` for `convdiff1` and `1.6e-5` for `convdiff3`).

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

    f_0(z) = ifelse(z == L, 1.0, 0.0)

    bcs = [
        f(0, z) ~ f_0(z),
        f(t, 0) ~ 0.0, f(t, L) ~ 1.0,
    ]
    # Space and time domains
    domains = [
        t ∈ Interval(0.0, 30.0),
        z ∈ Interval(0.0, L),
    ]

    # Analytic/reference solution. With γ = v/2k and μₙ = nπ/L, f = fₛ + w where
    # fₛ = expm1(-2γz)/expm1(-2γL) is the steady state and
    # w = Σₙ 2μₙ(-1)ⁿ e^{γ(L-z) - k(γ²+μₙ²)t} sin(μₙz) / (L(γ²+μₙ²))
    # solves the homogeneous problem with initial data -fₛ.
    maxiters = 50_000

    function steady(γ, z)
        # Rewritten for γ < 0 so no exponential exceeds 1.
        return γ < 0 ? exp(-2γ * (z - L)) * expm1(2γ * z) / expm1(2γ * L) :
            γ == 0 ? z / L : expm1(-2γ * z) / expm1(-2γ * L)
    end

    # The terms are bounded by e^{γ(L-z) - kμₙ²t}, so truncate where that is below e^{-40}.
    function seriesterms(k, γ, t, z)
        reach = 40 + max(Float64(γ) * (L - Float64(z)), 0)
        nterms = L / π * sqrt(reach / (Float64(k) * Float64(t)))
        nterms <= maxiters || throw(
            ArgumentError(
                "convection_diffusion analytic solution needs ≈$(round(Float64(nterms))) series terms at " *
                    "t = $(Float64(t)), z = $(Float64(z)) but is limited to $maxiters; " *
                    "use t ≥ $(L^2 * reach / (π^2 * maxiters^2 * Float64(k)))"
            )
        )
        return ceil(Int, nterms)
    end

    function transient(k, γ, t, z)
        T = typeof(γ)
        Lₜ = T(L)
        return sum(1:seriesterms(k, γ, t, z)) do n
            μ = n * T(π) / Lₜ
            2μ * (-1)^n * exp(γ * (Lₜ - z) - k * (γ^2 + μ^2) * t) * sin(μ * z) /
                (Lₜ * (γ^2 + μ^2))
        end
    end

    function u(ps, t, z)
        k = ps[1]
        γ = ps[2] / 2k
        t == 0 && return f_0(z)
        # For γ > 0 the terms are of size e^{γL} while the sum is O(1), so
        # carry enough extra bits to survive the cancellation.
        γ * L > 1 && return convdiff_fixed(L, k, γ, t, z, seriesterms(k, γ, t, z))
        return steady(γ, z) + transient(k, γ, t, z)
    end

    ref = [f(t, z) => u]

    tags = ["1D", "Dirichlet", "Linear", "Advection", "Diffusion", "Monotonic"]

    # PDE system
    convdiff = PDESystem(
        eq, bcs, domains, [t, z], [f(t, z)],
        [k => ps[1], v => ps[2]], analytic_func = ref, metadata = tags,
        name = name
    )

    return convdiff
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

    bcs = [
        u(0, z) ~ u_exact(0, z),
        u(t, 0) ~ u_exact(t, 0),
        u(t, 2π) ~ u_exact(t, 2π),
    ]

    # Space and time domains

    domains = [
        t ∈ Interval(0.0, 1.0),
        z ∈ Interval(0.0, 2π),
    ]

    # Analytic/reference solution
    ref = [u(t, z) ~ u_exact(t, z)]

    tags = ["1D", "Transport", "Linear", "Sinusoidal", "Inhomogeneous", "Advection"]

    # PDESystem

    @named trans_sin = PDESystem(
        eq, bcs, domains, [t, z], [u(t, z)], analytic = ref,
        metadata = tags, eval_module = @__MODULE__
    )

    return trans_sin
end
