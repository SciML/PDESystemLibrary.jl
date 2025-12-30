using PDESystemLibrary
PSL = PDESystemLibrary

using ModelingToolkit, MethodOfLines, DomainSets, OrdinaryDiffEq, NonlinearSolve, Test
using DomainSets: supremum, infimum

N = 100

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
            sol = NonlinerSolve.solve(prob, NewtonRaphsom())
            @test sol.retcode == SciMLBase.ReturnCode.Success
        else
            @parameters t
            disc = MOLFiniteDifference(dxs, t)
            prob = discretize(ex, disc)
            sol = solve(prob, FBDF())
            @test sol.retcode == SciMLBase.ReturnCode.Success
        end
    end
end
