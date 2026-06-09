# The test is simply that all of the examples build!
const GROUP = get(ENV, "GROUP", "All")

if GROUP == "QA"
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.develop(path = joinpath(@__DIR__, ".."))
    Pkg.instantiate()
    include("qa.jl")
else
    using SafeTestsets

    if GROUP == "All" || GROUP == "Core"
        @time @safetestset "Test against MethodOfLines.jl" begin
            include("mol_test.jl")
        end
    end

    # TODO: fix this when NeuralPDE.jl can be added to the test environment. (compat with Symbolics.jl 5)

    # if GROUP == "All" || GROUP == "NeuralPDE"
    #     @time @safetestset "Test against NeuralPDE.jl" begin include("neuralpde_test.jl") end
    # end
end
