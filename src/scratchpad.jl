using LinearAlgebra
using Plots
using Demo

poly_ord = 1
Nmat = shape_mat(0.0,poly_ord)
print(Nmat)
@time rvec =rhat_vec(1.0, 1.0, poly_ord, [1.0, 1.0, 0.4, 0.2])
print(rvec)

plot_lagrange_sf_or_deriv(1)
plot_lagrange_sf_or_deriv(1,true)
plot_lagrange_sf_or_deriv(2)
plot_lagrange_sf_or_deriv(2,true)
plot_lagrange_sf_or_deriv(5)
plot_lagrange_sf_or_deriv(5,true)