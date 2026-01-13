function general_system(α, β, γ, δ, name = :sys)
    @parameters t, x
    @variables u(..)

    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    Dxxx = Differential(x)^3
    Dxxxx = Differential(x)^4

    eq = Dt(u(t, x)) + α * Dx(u(t, x)) ~
        β * Dxx(u(t, x)) + γ * Dxxx(u(t, x)) -
        δ * Dxxxx(u(t, x))
    domain = [
        x ∈ Interval(0.0, 2π),
        t ∈ Interval(0.0, 3.0),
    ]

    ic_bc = [
        u(0.0, x) ~ cos(x)^2,
        u(t, 0.0) ~ u(t, 2π),
    ]

    analytic = [u(t, x) ~ 0.5 * (exp(-t * 4(β + 4δ)) * cos(t * (-8γ - 2α) + 2x) + 1)]

    return (eq, ic_bc, domain, [t, x], [u(t, x)], analytic, name)
end

function build_with_tags(sys, tags)
    eq, ic_bc, domain, indvars, depvars, analytic, name = sys
    sys = PDESystem(
        eq, ic_bc, domain, indvars, depvars, name = name,
        analytic = analytic, metadata = tags
    )
    return sys
end

adv = build_with_tags(
    general_system(1.0, 0.0, 0.0, 0.0, :adv),
    ["1D", "Advection", "Linear", "Periodic"]
)
diff = build_with_tags(
    general_system(0.0, 1.0, 0.0, 0.0, :diff),
    ["1D", "Diffusion", "Linear", "Periodic"]
)
adv3 = build_with_tags(
    general_system(0.0, 0.0, 1.0, 0.0, :adv3),
    ["1D", "3rd Order", "Linear", "Periodic"]
)
diff4 = build_with_tags(
    general_system(0.0, 0.0, 0.0, 1.0, :diff4),
    ["1D", "4th Order", "Linear", "Periodic"]
)

advdiff = build_with_tags(
    general_system(1.0, 1.0, 0.0, 0.0, :advdiff),
    ["1D", "Advection", "Diffusion", "Linear", "Periodic"]
)

advdiff3 = build_with_tags(
    general_system(1.0, 1.0, 1.0, 0.0, :advdiff3),
    [
        "1D",
        "Advection",
        "Diffusion",
        "3rd Order",
        "Linear",
        "Periodic",
    ]
)
advdiffno3 = build_with_tags(
    general_system(1.0, 1.0, 0.0, 1.0, :advdiffno3),
    [
        "1D",
        "Advection",
        "Diffusion",
        "4th Order",
        "Linear",
        "Periodic",
    ]
)

gen1 = build_with_tags(
    general_system(rand(4)..., :gen1),
    [
        "1D",
        "Advection",
        "Diffusion",
        "3rd Order",
        "4th Order",
        "Linear",
        "Periodic",
    ]
)
gen2 = build_with_tags(
    general_system(rand(4)..., :gen2),
    [
        "1D",
        "Advection",
        "Diffusion",
        "3rd Order",
        "4th Order",
        "Linear",
        "Periodic",
    ]
)
gen3 = build_with_tags(
    general_system(rand(4)..., :gen3),
    [
        "1D",
        "Advection",
        "Diffusion",
        "3rd Order",
        "4th Order",
        "Linear",
        "Periodic",
    ]
)

# Add to lists

gensys = [adv, diff, adv3, diff4, advdiff, advdiff3, advdiffno3, gen1, gen2, gen3]

all_systems = vcat(all_systems, gensys)
