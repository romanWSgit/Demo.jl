using Test
using Demo


@testset "SBFEM.jl" begin
    #@test greeting() === nothing

    # Test for hamiltonian matrix
    
    @test ishamiltonian([1 2; 2 -1]) == true
    @test ishamiltonian([22 9.3; 200 -1]) == false

    # Test
    @test lagrange_sample_points(1) == [-1.0; 1.0]
    @test lagrange_sample_points(2) == [-1.0; 0.0; 1.0]
    # Test lagrange_poly
    #   - test for linear case:
    @test lagrange_poly_or_deriv(-1.0, 1, [-1.0; 1.0]) == 1
    @test lagrange_poly_or_deriv(1.0, 1, [-1.0; 1.0]) == 0
    @test lagrange_poly_or_deriv(-1.0, 2, [-1.0; 1.0]) == 0
    @test lagrange_poly_or_deriv(1.0, 2, [-1.0; 1.0]) == 1
   
end