module PDESystemLibrary
using ModelingToolkit, DomainSets

using Markdown
using Random

Random.seed!(100)

all_systems = []

include("../lib/burgers.jl")
include("../lib/linear_diffusion.jl")
include("../lib/general_linear_system.jl")
include("../lib/brusselator.jl")

# Don't export anything, just add your systems to the lists and import it downstream.
end # module PDESystemLibrary
