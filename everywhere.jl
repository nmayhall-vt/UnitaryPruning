@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.instantiate()
@everywhere using UnitaryPruning

@everywhere dir = $dir
@everywhere include("read_in.jl")
