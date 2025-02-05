module Parameter
    # parameters of sketched Newton
    struct StoNewton
		# stopping parameters
        MaxIter::Int                       # Maximum iteration
        # fixed parameters
        Rep::Int                           # Number of independent runs
        tau::Int64                         # Number of iterations for inexact solver
        c_1::Float64                       # beta_t = c_1/t^{c_2}
        c_2::Float64
        c_3::Float64                       # chi_t = beta_t^{c_3}
        Sigma::Matrix{Float64}             # Sampling covariance matrix
        # data parameters
        D::Int64                           # Problem dimension
		sigma::Float64                     # Noise level
    end
end

# Equi-corr matrix
D = 5
RR = 0.3
Sigma = RR*ones(D,D)+(1-RR)*Matrix(1.0I,D,D)


StoNewtonSet = Parameter.StoNewton(3e5,              # Max_Iter
	                  200,                           # Rep
	                  0,                             # tau
					  1,                             # c_1
	                  0.505,                         # c_2
					  2,                             # c_3
					  Sigma,                         # Sigma
					  D,                             # d
					  1)                             # sigma
