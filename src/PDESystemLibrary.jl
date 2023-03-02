module PDESystemLibrary
using ModelingToolkit, DomainSets
using OrdinaryDiffEq
using Interpolations
using Symbolics
using Symbolics: unwrap

using IfElse
using Markdown
using Random

Random.seed!(100)

all_systems = []

include("../lib/burgers.jl")
include("../lib/linear_diffusion.jl")
include("../lib/linear_convection.jl")
include("../lib/nonlinear_diffusion.jl")
include("../lib/general_linear_system.jl")
include("../lib/brusselator.jl")

function get_pdesys_with_tags(tags...; f = all)
    filter(all_systems) do ex
        f(t -> t in ex.metadata, tags)
    end
end

export get_pdesys_with_tags
end # module PDESystemLibrary
