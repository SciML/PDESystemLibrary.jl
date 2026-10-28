using PDESystemLibrary
using Test

@testset "Precompiled analytic solutions" begin
    systems = get_pdesys_with_tags(String[])
    system = only(filter(system -> system.name == :inviscid_burgers_monotonic, systems))
    analytic = only(values(system.analytic_func))

    @test analytic([], 0.0, 0.5) == 0.5
end

@testset "analytic_func does not throw for any system" begin
    systems = PDESystemLibrary.all_systems
    @test length(systems) > 400
    for system in systems
        ps = system.ps isa AbstractVector ? [p isa Pair ? last(p) : 0.0 for p in system.ps] : []
        coords = ntuple(_ -> 0.5, length(system.ivs))
        for analytic in values(system.analytic_func)
            @test (analytic(ps, coords...); true)
        end
    end
end
