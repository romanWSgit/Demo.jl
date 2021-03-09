using Test
using Demo


@testset "SBFEM.jl" begin
    # Test for hamiltonian matrix
    @test greeting() == nothing
    @test ishamiltonian([1 2; 2 -1]) == true
    @test ishamiltonian([22 9.3; 200 -1]) == false
end
