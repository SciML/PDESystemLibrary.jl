using PDESystemLibrary
using Test

@testset "Precompiled analytic solutions" begin
    systems = get_pdesys_with_tags(String[])
    system = only(filter(system -> system.name == :inviscid_burgers_monotonic, systems))
    analytic = only(values(system.analytic_func))

    @test analytic([], 0.0, 0.5) == 0.5
end

@testset "all_systems registry" begin
    names = [system.name for system in PDESystemLibrary.all_systems]
    @test :trans_sin in names
    @test allunique(names)
end
