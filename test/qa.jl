using PDESystemLibrary, Aqua, JET, Test

@testset "Aqua" begin
    Aqua.test_all(PDESystemLibrary)
end

@testset "JET" begin
    JET.test_package(PDESystemLibrary; target_defined_modules = true)
end
