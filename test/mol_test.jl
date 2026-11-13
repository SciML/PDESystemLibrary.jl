using PDESystemLibrary
PSL = PDESystemLibrary

using ModelingToolkit, MethodOfLines, DomainSets, OrdinaryDiffEq, OrdinaryDiffEqSDIRK,
    NonlinearSolve, Test
using DomainSets: supremum, infimum

N = 100

# `:adv3` is pure dispersive, so its central-difference discretization has a
# near-imaginary spectrum that needs TRBDF2's implicit stability.
# `:advdiff3` needs a Rosenbrock method for the same reason despite its diffusion.
const EXAMPLE_ALGORITHMS = Dict(:adv3 => TRBDF2(), :advdiff3 => Rodas5P())

for ex in PSL.all_systems
    @testset "Example: $(ex.name)" begin
        # MethodOfLines' stationary 2D path currently returns `Stalled` from NewtonRaphson, even
        # for its own `u = x*y` Laplace example, so steady systems are only checked analytically.
        "Steady" in ex.metadata && continue
        ivs = filter(x -> !isequal(Symbol(x), :t), ex.ivs)
        dxs = map(ivs) do x
            xdomain = ex.domain[findfirst(d -> isequal(x, d.variables), ex.domain)]
            x => (supremum(xdomain.domain) - infimum(xdomain.domain)) /
                (floor(N^(1 / length(ivs))) - 1)
        end
        if length(ivs) == 0
            continue
        elseif length(ivs) == length(ex.ivs)
            disc = MOLFiniteDifference(dxs)
            prob = discretize(ex, disc)
            sol = NonlinearSolve.solve(prob, NewtonRaphson())
            @test sol.retcode == SciMLBase.ReturnCode.Success
        else
            @parameters t
            disc = MOLFiniteDifference(dxs, t)
            sys, tspan = symbolic_discretize(ex, disc)
            prob = ODEProblem(mtkcompile(sys), nothing, tspan)
            sol = solve(prob, get(EXAMPLE_ALGORITHMS, ex.name, FBDF()))
            @test sol.retcode == SciMLBase.ReturnCode.Success
        end
    end
end
