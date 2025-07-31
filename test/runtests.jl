# The test is simply that all of the examples build!
using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")
const is_APPVEYOR = Sys.iswindows() && haskey(ENV, "APPVEYOR")
const is_TRAVIS = haskey(ENV, "TRAVIS")

if GROUP == "All" || GROUP == "MOL"
    @time @safetestset "Test against MethodOfLines.jl" begin
        include("mol_test.jl")
    end
end

# TODO: fix this when NeuralPDE.jl can be added to the test environment. (compat with Symbolics.jl 5)

# if GROUP == "All" || GROUP == "NeuralPDE"
#     @time @safetestset "Test against NeuralPDE.jl" begin include("neuralpde_test.jl") end
# end
