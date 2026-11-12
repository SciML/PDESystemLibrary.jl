using PDESystemLibrary, Test
using SciMLTesting
# ExtFloat (lib/linear_convection.jl) calls MPFR through Base's non-public libmpfr handle.
run_qa(
    PDESystemLibrary;
    ei_kwargs = (; all_qualified_accesses_are_public = (; ignore = (:MPFR, :libmpfr)))
)
