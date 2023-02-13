using PDESystemLibrary
PSL = PDESystemLibrary

using ModelingToolkit, MethodOfLines, DomainSets, OrdinaryDiffEq, NonlinearSolve

for ex in PSL.all_systems
    @testset "Example with MethodOfLines.jl: $(ex.name)\n Equations: $(ex.eqs) \nBCs/ICs: $(ex.bcs)" begin
        ivs = filter(x -> !isequal(Symbol(x), :t), ex.ivs)
        dxs = map(x -> x => (supremum(d) - infimum(d))/10, ivs)
        if length(ivs) == 0
            continue
        elseif length(ivs) == length($(ex.ivs))
            disc = MOLFiniteDifference(dxs)
            prob = discretize($ex, disc)
            sol = NonlinerSolve.solve(prob, NewtonRaphsom())
        else
            @parameters t
            disc = MOLFiniteDifference(dxs, t)
            prob = discretize($ex, disc)
            sol = solve(prob, FBDF())
        end
    end
end
