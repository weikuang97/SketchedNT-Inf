module Parameter
    # Parameters of ASGD method
    struct SGD
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
RR = 0.2
Sigma = RR*ones(D,D)+(1-RR)*Matrix(1.0I,D,D)

# Sigma = Matrix(1.0I,D,D)
# Sigma = [RR^abs(i-j) for i=1:D, j=1:D]


SGDSet = Parameter.SGD(3e5,                        # Max_Iter
	                  200,                           # Rep
					  1,                           # c_1
	                  1,                         # c_2
					  Sigma,                         # Sigma
					  D,                             # d
					  1)                             # sigma
