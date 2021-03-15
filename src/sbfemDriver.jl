#=
Main Sbfem functionality
=#
using LinearAlgebra



@doc raw"""
    shape_mat(η, poly_ord) -> Matrix

returns the shape matrix for standard shape functions

# Arguments
- `η`: the circumferential variable.
- `poly_ord`: shape function order

# Examples
```julia-repl
julia> shape_mat(1.0, 1)
2×4 Array{Float64,2}:
 -0.0  -0.0  1.0  0.0
 -0.0  -0.0  0.0  1.0
```
"""
function shape_mat(η::Float64, poly_ord::Int)::Matrix
    lsp = lagrange_sample_points(poly_ord)
    Nmat = lagrange_poly(η,1,lsp) * Matrix(I,2,2)
    for i in 1:(poly_ord + 1)
        if i > 1
            Nmat = hcat(Nmat, lagrange_poly(η,i,lsp)*Matrix(I,2,2))
        end
    end
    return Nmat
end

@doc raw"""
    isHamiltonian(matrix, [error]) -> Bool

Checks if the `matrix` is a hamiltonian matix.

# Arguments
- `matrix`: the matrix to check.
- `error`: the value an entry of `matrix` can defer from 0

# Examples
```julia-repl
julia> ishamiltonian(matrix::[1 2; 2 -1], error::Float64 = 1e-12)
true
```
"""
function rhat_vec(ξ::Real, η::Real, poly_ord ,coord_vec::Vector, scaling_centre::Vector = [0.0, 0.0])::Vector
    return (ξ * (shape_mat(η, poly_ord)*coord_vec) + scaling_centre)
end
