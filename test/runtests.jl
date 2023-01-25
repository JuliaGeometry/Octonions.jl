using Test
using Octonions
using Aqua

Aqua.test_all(Octonions)

@testset "Octonions.jl" begin
    include("helpers.jl")
    include("octonion.jl")
end
