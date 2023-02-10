using PDESystemLibrary
PSL = PDESystemLibrary

using MethodOfLines, DomainSets, OrdinaryDiffEq, NonlinearSolve

for ex in PSL.all_systems
    @testset "Example: $(ex.name)" begin
        ivs = filter(x -> !isequal(Symbol(x), :t), ex.ivs)
        dxs = map(x -> x => (supremum(d) - infimum(d))/10, ivs)
        if length(ivs) == 0
            continue
        elseif length(ivs) == length(ex.ivs)
            disc = MOLFiniteDifference(dxs)
            prob = discretize(ex, disc)
            sol = NonlinerSolve.solve(prob, NewtonRaphsom())
        else
            @parameters t
            disc = MOLFiniteDifference(dxs, t)
            prob = discretize(ex, disc)
            sol = solve(prob, FBDF())
        end
    end
end
