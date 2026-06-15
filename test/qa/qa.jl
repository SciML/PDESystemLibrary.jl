using PDESystemLibrary, Aqua, JET
using SciMLTesting

run_qa(
    PDESystemLibrary; Aqua = Aqua, JET = JET, jet = true,
    jet_kwargs = (; target_defined_modules = true)
)
