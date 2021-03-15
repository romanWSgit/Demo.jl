#=
mandatory math stuff
=#
using LinearAlgebra
using Plots

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





@doc raw"""
    lagrange_poly(atPoint, atIndex, atDataPoints) -> Float64

Returns the value of the Lagrange

# Arguments
- `atPoint`: the matrix to check.
- `atIndex`: the value an entry of `matrix` can defer from 0
- `atDataPoints`:
# Examples
```julia-repl
julia> 
```
"""
function lagrange_sample_points(polyOrd)::Vector{Float64}
    return collect(-1.0:(2.0/polyOrd):1.0)
end


@doc raw"""
    lagrange_poly(atPoint, atIndex, atDataPoints) -> Float64

Returns the value of the Lagrange

# Arguments
- `atPoint`: the matrix to check.
- `atIndex`: the value an entry of `matrix` can defer from 0
- `atDataPoints`:
# Examples
```julia-repl
julia> 
```
"""
function lagrange_poly(atPoint::Float64, atIndex::Int64, atDataPoints::Vector{Float64})::Float64
	x, i, xm = atPoint, atIndex, atDataPoints
	y = 1
	n::Int64 = length(xm)
	for index = 1:n
		if i != index 
			y *= (x - xm[index]) / (xm[i] - xm[index])
		end
	end
	return y
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
function lagrange_poly_or_deriv(atPoint::Float64, atIndex::Int64, atDataPoints::Vector{Float64},
     calcDeriv::Bool = false)::Float64
	x, i, xm, deriv = atPoint, atIndex, atDataPoints, calcDeriv
    n::Int64  = length(xm)
    if deriv
        k = 0
        y = 0
        for l = 1:n 
            if l != i 
                k = 1 / (xm[i] - xm[l])
                for m = 1:n 
                    if (m != i) && (m != l) 
                        k = k * ((x-xm[m])/(xm[i]-xm[m]))
					end
				end
                y = y + k
			end
		end
    else
        y = 1
        for index = 1:n 
            if i != index 
                y *= (x - xm[index]) / (xm[i] - xm[index])
			end
        end
	end
    return y
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
function plot_lagrange_sf(polyOrd::Int64) 
	x = collect(-1.0:0.01:1.0)
	y = []
	arr = collect(-1.0:(2.0/polyOrd):1.0)
	for i = 1:(polyOrd + 1)
		y1 = map(s->lagrange_poly(s,i,arr), x)
		if i == 1
			y = y1
		else
			y = hcat(y, y1) 
		end
	end
	plot(x, y)
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
function plot_lagrange_sf_or_deriv(polyOrd::Int64, deriv::Bool = false) 
	x = collect(-1.0:0.01:1.0)
	y = []
	arr = collect(-1.0:(2.0/polyOrd):1.0)
	for i = 1:(polyOrd + 1)
		y1 = map(s->lagrange_poly_or_deriv(s,i,arr,deriv), x)
		if i == 1
			y = y1
		else
			y = hcat(y, y1) 
		end
	end
	plot(x, y)
end