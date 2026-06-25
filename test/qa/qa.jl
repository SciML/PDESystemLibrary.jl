using PDESystemLibrary, Test
using SciMLTesting
using JET

# JET runs in `:typo` mode (the SciMLTesting `run_qa` default). The library's
# problem builders construct symbolic expressions from `@variables`/`@parameters`;
# `Symbolics.wrap` has a union return type that JET widens to include the callable
# wrapper `Symbolics.CallAndWrap` (only produced by `@variables u(..)`), so basic-mode
# analysis emits spurious `no matching method found` union-split reports for `+/-/*/cos/sin`
# on scalar variables that are always `Num` at runtime. Typo mode still catches real
# undefined-name errors without these dependency-inference false positives.
#
# ExplicitImports findings:
#  * no_implicit_imports: tracked-broken. The module `using`s the symbolic/solver
#    umbrellas (ModelingToolkit, DomainSets, OrdinaryDiffEq, Interpolations) that
#    re-export from Symbolics/ModelingToolkitBase/IntervalSets/SciMLBase/CommonSolve,
#    so a `using X: name` rewrite would pin owner packages that reshuffle names across
#    releases. Two entries (`limit`, `domain`) are also local definitions EI flags
#    conservatively, so a blanket explicit-import rewrite would be wrong. Tracked in
#    https://github.com/SciML/PDESystemLibrary.jl/issues/63 ; kept broken
#    (auto-flags Unexpected Pass once resolved).
#  * all_qualified_accesses_are_public: `Random.seed!` is the canonical seeding API
#    but the Random stdlib does not declare it `public`; ignored.
run_qa(
    PDESystemLibrary; explicit_imports = true,
    jet_kwargs = (; target_modules = (PDESystemLibrary,), mode = :typo),
    ei_kwargs = (;
        all_qualified_accesses_are_public = (; ignore = (:seed!,)),  # Random.seed! (Random stdlib, not declared public)
    ),
    ei_broken = (:no_implicit_imports,)
)
