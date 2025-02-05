module Parameter
    # Parameters of sketched Newton method
    struct StoNewton
        # stopping parameters
        MaxIter::Int                       # Maximum iteration
        # fixed parameters
        Rep::Int                           # Number of independent runs
        tau::Array{Int64}                  # Number of iterations for inexact solver
        c_1::Float64                       # beta_t = c_1/t^{c_2}
        c_2::Float64
        c_3::Float64                       # chi_t = beta_t^{c_3}
        SigToe::Array{Float64}             # Toeplitz covariance parameter
        SigEqui::Array{Float64}            # Equi-Corr covariance parameter
        # data parameters
	        D::Array{Int64}                # Problem dimension
		mu::Float64                        # Regularization parameter
    end
end

StoNewtonSet = Parameter.StoNewton(3e5,              # Max_Iter
	                  200,                           # Rep
	                  [0,10,20,40],                  # tau
					  1,                             # c_1
	                  0.505,                         # c_2
					  2,                             # c_3
					  [0,0.4,0.5,0.6],               # RToe
					  [0.1,0.2,0.3],                 # REqu
					  [20,40,60,100],                # d
					  0.1)                           # mu
