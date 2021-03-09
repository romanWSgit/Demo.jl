#=
mandatory math stuff
=#
using LinearAlgebra
"""
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
function ishamiltonian(matrix::Matrix, error::Float64 = 1e-12)
    m = matrix
    dim = size(m)[1] ÷ 2
    i = Matrix(1.0*I, dim, dim)
    z = Matrix(0.0*I, dim, dim)
    JN2 = hvcat(2, z, i, -i, z)
    Res = transpose(JN2 * m) - (JN2 * m)
    checkHamiltonian = maximum(Res) < error && minimum(Res) > -error
    return checkHamiltonian
end
