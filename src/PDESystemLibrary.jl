module PDESystemLibrary
using CommonSolve: solve
using Interpolations: Gridded, Linear, Periodic, extrapolate, interpolate
using IntervalSets: (..), Interval
using ModelingToolkitBase: @named, @parameters, PDESystem
using OrdinaryDiffEqSDIRK: TRBDF2
import RuntimeGeneratedFunctions
using SciMLBase: ODEProblem
using Symbolics: @variables, Differential, arguments, build_function, unwrap

RuntimeGeneratedFunctions.init(@__MODULE__)

all_systems = []

# ModelingToolkitBase < 1.78.1 splices `Pair` parameters into a destructuring assignment when it
# builds `analytic_func` from `analytic`, so for such systems build the functions here.
function analytic_funcs(analytic, ps)
    return map(analytic) do eq
        ex = build_function(eq.rhs, first.(ps)..., arguments(unwrap(eq.lhs))...; expression = Val{true})
        f = RuntimeGeneratedFunctions.RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, ex)
        eq.lhs => (p, args...) -> f(p..., args...)
    end
end

include("../lib/burgers.jl")
include("../lib/linear_diffusion.jl")
include("../lib/linear_convection.jl")
include("../lib/nonlinear_diffusion.jl")
include("../lib/general_linear_system.jl")
include("../lib/brusselator.jl")
include("../lib/advection_dominated.jl")

"""
    get_pdesys_with_tags(withtags; without = [], f = all)

Return PDE systems whose `metadata` tags match `withtags`.

## Arguments

- `withtags`: An iterable of tags to include in the result.

## Keywords

- `without = []`: Tags that exclude a system from the result.
- `f = all`: Predicate that combines matches for `withtags`. Use `any` to select
  systems matching at least one requested tag.

# Examples

```julia
diffusion_systems = get_pdesys_with_tags(["Diffusion"])
one_dimensional_or_heat = get_pdesys_with_tags(["1D", "Heat"]; f = any)
smooth_diffusion = get_pdesys_with_tags(["Diffusion"]; without = ["Discontinuous"])
```
"""
function get_pdesys_with_tags(withtags; without = [], f = all)
    return filter(all_systems) do ex
        b = f(t -> t ∈ ex.metadata, withtags)
        b && all(t -> t ∉ ex.metadata, without)
    end
end

export get_pdesys_with_tags

include("precompile.jl")
end # module PDESystemLibrary
