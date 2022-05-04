using UnitaryPruning 
using Test
using Random

Random.seed!(1234567)

@testset "UnitaryPruning" begin
    include("test_PauliStrings.jl")
    include("test_PauliMask.jl")
    include("test_dfs.jl")
end

