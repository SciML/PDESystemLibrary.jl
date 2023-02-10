module PDESystemLibrary
using ModelingToolkit, DomainSets

using Markdown
using Random

Random.seed!(100)

all_systems = []
linear_diffusion_systems = []
convection_systems = []
advection_systems = []
nonlinear_systems = []

include("../lib/burgers.jl")
include("../lib/linear_diffusion.jl")

# Don't export anything
end # module PDESystemLibrary
