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
#
# `:advdiff3` is also unstable, but only on Julia >= 1.13: on earlier versions
# the integrator's adaptive step keeps it under control, while changes to
# floating-point tie-breaking / lowering on 1.13 push the solver over the
# stiffness threshold and FBDF aborts with `dt < eps`. Tracked conditionally
# so we keep regression coverage on lts/1.x and don't false-fail on 1.13+.
const BROKEN_EXAMPLES = if VERSION >= v"1.13-"
    Set([:adv3, :advdiff3])
else
    Set([:adv3])
end

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
            if ex.name in BROKEN_EXAMPLES
                # solve itself may throw on numerically-unstable examples;
                # wrap in try/catch so @test_broken sees a `false` instead of
                # an Error when that happens.
                @test_broken try
                    NonlinearSolve.solve(prob, NewtonRaphson()).retcode ==
                        SciMLBase.ReturnCode.Success
                catch
                    false
                end
            else
                sol = NonlinearSolve.solve(prob, NewtonRaphson())
                @test sol.retcode == SciMLBase.ReturnCode.Success
            end
        else
            @parameters t
            disc = MOLFiniteDifference(dxs, t)
            prob = discretize(ex, disc)
            if ex.name in BROKEN_EXAMPLES
                @test_broken try
                    solve(prob, FBDF()).retcode == SciMLBase.ReturnCode.Success
                catch
                    false
                end
            else
                sol = solve(prob, FBDF())
                @test sol.retcode == SciMLBase.ReturnCode.Success
            end
        end
    end
end
