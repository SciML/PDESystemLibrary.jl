function general_system(α, β, γ, δ, name = :sys)
    @parameters t, x
    @variables u(..)

    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    Dxxx = Differential(x)^3
    Dxxxx = Differential(x)^4

    eq = Dt(u(t, x)) + α * Dx(u(t, x)) ~ β * Dxx(u(t, x)) + γ * Dxxx(u(t, x)) -
                                         δ * Dxxxx(u(t, x))
    domain = [x ∈ Interval(0.0, 2π),
        t ∈ Interval(0.0, 3.0)]

    ic_bc = [u(0.0, x) ~ cos(x)^2,
        u(t, 0.0) ~ u(t, 2π)]

    function analytic(t, x; ps = nothing)
        0.5 * (exp(-t * 4(β + 4δ)) * cos(t * (-8γ - 2α) + 2x) + 1)
    end

    sys = PDESystem(eq, ic_bc, domain, [t, x], [u(t, x)],
                    analytic = analytic, name = name)
end

adv = general_system(1.0, 0.0, 0.0, 0.0, :adv)
diff = general_system(0.0, 1.0, 0.0, 0.0, :diff)
adv3 = general_system(0.0, 0.0, 1.0, 0.0, :adv3)
diff4 = general_system(0.0, 0.0, 0.0, 1.0, :diff4)

advdiff = general_system(1.0, 1.0, 0.0, 0.0, :advdiff)
advdiff3 = general_system(1.0, 1.0, 1.0, 0.0, :advdiff3)
advdiffno3 = general_system(1.0, 1.0, 0.0, 1.0, :advdiffno3)

gen1 = general_system(rand(4)..., :gen1)
gen2 = general_system(rand(4)..., :gen2)
gen3 = general_system(rand(4)..., :gen3)

# Add to lists

gensys = [adv, diff, adv3, diff4, advdiff, advdiff3, advdiffno3, gen1, gen2, gen3]

all_systems = vcat(all_systems, gensys)

push!(linear_diffusion_systems, diff)
push!(linear_diffusion_systems, diff4)

push!(convection_systems, adv)

push!(advection_systems, adv)
