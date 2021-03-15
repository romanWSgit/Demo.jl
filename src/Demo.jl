module Demo
include("math.jl")

export ishamiltonian
export lagrange_sample_points
export lagrange_poly
export lagrange_poly_or_deriv
export plot_lagrange_sf
export plot_lagrange_sf_or_deriv

include("sbfemDriver.jl")
export shape_mat
export rhat_vec

# greeting() = print("Hello World!")
# export greeting

end # module
