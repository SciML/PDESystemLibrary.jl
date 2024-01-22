# The test is simply that all of the examples build!
using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "MOL"
    @time @safetestset "Test against MethodOfLines.jl" begin include("mol_test.jl") end
end

if GROUP == "All" || GROUP == "NeuralPDE"
    @time @safetestset "Test against NeuralPDE.jl" begin include("neuralpde_test.jl") end
end
