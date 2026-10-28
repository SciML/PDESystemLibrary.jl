using PDESystemLibrary
using Test

@testset "all_systems registry" begin
    names = [system.name for system in PDESystemLibrary.all_systems]
    @test :trans_sin in names
    @test allunique(names)
end
