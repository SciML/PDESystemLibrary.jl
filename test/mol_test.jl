using PDESystemLibrary
PSL = PDESystemLibrary

using ModelingToolkit, MethodOfLines, DomainSets, OrdinaryDiffEq, NonlinearSolve, Test
using DomainSets: supremum, infimum

N = 100

# Examples that are known to be numerically unstable with the default FBDF
# solver and central-difference spatial discretization. These are typically
# pure-dispersive equations (e.g. `Dt(u) = Dxxx(u)`) where the spatial
# discretization yields a near-imaginary spectrum that BDF methods handle
# poorly. These are tracked as broken so the test suite remains green; fixing
# them requires either an upwind/finite-volume discretization or an explicit
# solver such as `Vern9()`.
const BROKEN_EXAMPLES = Set([:adv3])

for ex in PSL.all_systems
    @testset "Example: $(ex.name)" begin
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
            if ex.name in BROKEN_EXAMPLES
                @test_broken sol.retcode == SciMLBase.ReturnCode.Success
            else
                @test sol.retcode == SciMLBase.ReturnCode.Success
            end
        else
            @parameters t
            disc = MOLFiniteDifference(dxs, t)
            prob = discretize(ex, disc)
            sol = solve(prob, FBDF())
            if ex.name in BROKEN_EXAMPLES
                @test_broken sol.retcode == SciMLBase.ReturnCode.Success
            else
                @test sol.retcode == SciMLBase.ReturnCode.Success
            end
        end
    end
end
