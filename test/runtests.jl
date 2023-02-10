 # The test is simply that all of the examples build!
 using SafeTestSets

const GROUP = get(ENV, "GROUP", "All")
const is_APPVEYOR = Sys.iswindows() && haskey(ENV, "APPVEYOR")
const is_TRAVIS = haskey(ENV, "TRAVIS")

if GROUP == "All" || GROUP == "MOL"
    @time @safetestsets begin include("mol_test.jl") end
end

if GROUP == "All" || GROUP == "NeuralPDE"
    @time @safetestsets begin include("neuralpde_test.jl") end
end
