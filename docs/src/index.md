# PDESystemLibrary.jl

Load `PDEBase` to make the unit-square heat equation available through `get_pdesys_with_tags(["MeshDomain"])`. MethodOfLines also loads PDEBase. The example provides a symbolic model and an analytic solution. Discretization requires a backend mesh. To omit mesh-domain examples from a MethodOfLines query, pass `without = ["MeshDomain"]`.

```@docs
get_pdesys_with_tags
```
