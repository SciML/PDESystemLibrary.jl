using PDESystemLibrary
using Test

@testset "Precompiled analytic solutions" begin
    systems = get_pdesys_with_tags(String[])
    system = only(filter(system -> system.name == :inviscid_burgers_monotonic, systems))
    analytic = only(values(system.analytic_func))

    @test analytic([], 0.0, 0.5) == 0.5
end

@testset "windowlower and dconv analytic solutions" begin
    windowlower = PDESystemLibrary.windowlower
    @test windowlower(5.0, -0.1) == 0.0
    @test windowlower(5.0, 0.0) == 5.0
    @test windowlower(5.0, 0.5) == 5.0
    @test windowlower(5.0, 1.0) == 5.0
    @test windowlower(5.0, 1.1) == 0.0

    systems = get_pdesys_with_tags(String[])
    system = first(filter(system -> system.name == :dconvcos, systems))
    analytic = only(values(system.analytic_func))
    @test analytic([2], 0.0, 0.25) ≈ cos(2pi * 0.25) atol = 1.0e-12
    @test analytic([2], 0.0, 0.1) ≈ cos(2pi * 0.1)
    @test analytic([2], 0.25, 0.8) ≈ cos(2pi * 0.3)  # x - v*t = 0.3
    @test analytic([2], 0.25, 0.3) == 0.0            # x - v*t < 0: inflow Dirichlet BC
end

@testset "ddconv systems use linear_convection_dirichlet2" begin
    windowupper = PDESystemLibrary.windowupper
    @test windowupper(5.0, -0.1) == 0.0
    @test windowupper(5.0, 0.0) == 5.0
    @test windowupper(5.0, 1.0) == 5.0
    @test windowupper(5.0, 1.1) == 0.0

    systems = get_pdesys_with_tags(String[])
    byname = Dict(system.name => system for system in systems)
    fs = Dict(
        "sin" => x -> sin(2pi * x), "sin2" => x -> sin(2pi * x),
        "cos" => x -> cos(2pi * x),
        "tri" => x -> 1.0 - abs(x - floor(x + 0.5)),
        "square" => x -> ifelse(x - floor(x) < 0.5, 1.0, -1.0),
    )
    h = 1.0e-6
    for (suffix, f) in fs
        d = byname[Symbol("dconv", suffix)]
        dd = byname[Symbol("ddconv", suffix)]
        u = only(values(dd.analytic_func))
        p = [last(only(dd.ps))]
        v = only(p)

        # Dt u = +v Dx u (dconv has -v Dx u), Dirichlet boundary on the inflow side x = 1
        @test !isequal(dd.eqs[1], d.eqs[1])
        @test !occursin("-Differential", string(dd.eqs[1]))
        @test occursin("1.0", string(dd.bcs[2]))

        for x in 0.05:0.1:0.95
            @test u(p, 0.0, x) ≈ f(x) atol = 1.0e-12
        end
        for t in (0.1, 0.5) ./ v
            @test u(p, t, 1.0 + 1.0e-3) == 0.0
        end

        # u = f(x + v t) inside the domain, and it satisfies Dt u = v Dx u
        for (t, x) in ((0.05, 0.3), (0.1, 0.43), (0.2, 0.2))
            t /= v
            @test u(p, t, x) ≈ f(x + v * t) atol = 1.0e-12
            if suffix != "square" && suffix != "tri"
                residual = (u(p, t + h, x) - u(p, t - h, x)) / 2h -
                    v * (u(p, t, x + h) - u(p, t, x - h)) / 2h
                @test abs(residual) < 1.0e-4
            end
        end

        # the wave moves left, so it differs from the dconv solution
        t, x = 0.2 / v, 0.3
        @test u(p, t, x) != only(values(d.analytic_func))(p, t, x)
    end
end
