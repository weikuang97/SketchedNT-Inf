module Parameter
    # Parameters of StoSQP with augmented Lagrangian
    struct StoSQP
		# stopping parameters
        MaxIter::Int64                     # Maximum iteration
        # fixed parameters
        Rep::Int64                         # Number of independent runs
        tau::Array{Int64}                  # Number of iterations for inexact solver
        c_1::Float64                       # beta_t = c_1/t^{c_2}
        c_2::Float64
        c_3::Float64                       # chi_t = beta_t^{c_3}
        sigma2::Array{Float64}             # Sampling variance
		# problem set
		Prob::String					   # Problem set
    end
end



StoSQPSet = Parameter.StoSQP(3e5,                    # Max_Iter
	                  200,                           # Rep
	                  [40],                          # tau
					  1,                             # c_1
	                  0.505,                         # c_2
					  2,                             # c_3
					  [1,1e-1,1e-2,1e-4],		     # sigma2
					  "MARATOS")                     # problem
