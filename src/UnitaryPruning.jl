```@contents
```
module UnitaryPruning

using Printf
using BenchmarkTools
using DataStructures
using InteractiveUtils
using PauliOperators


include("stochastic_evolution.jl")
include("deterministic_evolution.jl")
include("eagle.jl")
include("unitary_sequence.jl")
include("energy_dfs.jl")


end # module
