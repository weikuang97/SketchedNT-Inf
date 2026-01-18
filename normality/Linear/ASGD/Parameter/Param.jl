module Parameter
    # Parameters of ASGD method
    struct ASGD
		# stopping parameters
        MaxIter::Int                       # Maximum iteration
        # fixed parameters
        Rep::Int                           # number of independent runs
        c_1::Float64                       # beta_t = c_1/t^{c_2}
        c_2::Float64
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


ASGDSet = Parameter.ASGD(3e5,                        # Max_Iter
	                  5,                           # Rep
					  0.5,                           # c_1
	                  0.505,                         # c_2
					  Sigma,                         # Sigma
					  D,                             # d
					  1)                             # sigma
