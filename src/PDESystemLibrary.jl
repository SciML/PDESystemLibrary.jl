module PDESystemLibrary
using ModelingToolkit, DomainSets
using OrdinaryDiffEq
using Interpolations

import SciMLBase

using IfElse
using IfElse: ifelse
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

function get_pdesys_with_tags(withtags; without = [], f = all)
    return filter(all_systems) do ex
        b = f(t -> t ∈ ex.metadata, withtags)
        b && all(t -> t ∉ ex.metadata, without)
    end
end

export get_pdesys_with_tags
end # module PDESystemLibrary
