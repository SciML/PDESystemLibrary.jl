using PDESystemLibrary, BenchmarkTools

const SUITE = BenchmarkGroup()

# =============================================================================
# PDE system library queries
# =============================================================================

SUITE["query"] = BenchmarkGroup()

SUITE["query"]["diffusion_all"] = @benchmarkable get_pdesys_with_tags(["Diffusion"])
SUITE["query"]["heat_or_1d"] = @benchmarkable get_pdesys_with_tags(
    ["1D", "Heat"]; f = any
)
SUITE["query"]["smooth_diffusion"] = @benchmarkable get_pdesys_with_tags(
    ["Diffusion"]; without = ["Discontinuous"]
)
SUITE["query"]["all"] = @benchmarkable get_pdesys_with_tags(String[])
