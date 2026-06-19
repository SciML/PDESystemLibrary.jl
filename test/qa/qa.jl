using PDESystemLibrary, Aqua, JET
using SciMLTesting

# JET runs in `:typo` mode (the SciMLTesting `run_qa` default). The library's
# problem builders construct symbolic expressions from `@variables`/`@parameters`;
# `Symbolics.wrap` has a union return type that JET widens to include the callable
# wrapper `Symbolics.CallAndWrap` (only produced by `@variables u(..)`), so basic-mode
# analysis emits spurious `no matching method found` union-split reports for `+/-/*/cos/sin`
# on scalar variables that are always `Num` at runtime. Typo mode still catches real
# undefined-name errors without these dependency-inference false positives.
run_qa(
    PDESystemLibrary; Aqua = Aqua, JET = JET, jet = true,
    jet_kwargs = (; target_modules = (PDESystemLibrary,), mode = :typo)
)
