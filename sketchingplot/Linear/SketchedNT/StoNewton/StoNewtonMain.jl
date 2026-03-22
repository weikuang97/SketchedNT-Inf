include("trueCov.jl")
include("StoNewton.jl")
## Implement sketched Newton inference framework for whole problem set
# StoNewtonSet: parameters

function StoNewtonMain(StoNewtonSet)
	# load parameters
	Max_Iter,Rep,tau = StoNewtonSet.MaxIter,StoNewtonSet.Rep,StoNewtonSet.tau
	c_1,c_2,c_3 = StoNewtonSet.c_1,StoNewtonSet.c_2,StoNewtonSet.c_3
	Sigma,sigma = StoNewtonSet.Sigma,StoNewtonSet.sigma
	D = StoNewtonSet.D



	q_vec = [1]


	# go over all repetitions
	for q in q_vec
		for ttau in tau
			nx = D
			# true parameters and true covariance matrix
		    X_true = (1/nx)*ones(nx)
			Xistar = trueCov(Sigma, nx, ttau, sigma)

		    # go over all repetitions
			for IdRep = 1:Rep
				println("q-", q, "-tau-", ttau, "-Rep-", IdRep)
				Time,Err_SC,IDCov_SC = StoNewton(c_1,c_2,c_3,Max_Iter,ttau,nx,X_true,Sigma,Xistar,q,sigma)
		        println("Time:", Time)

				path1 = string("../Solution/q",q,"/tau",ttau,"/rep",IdRep)
				if !isdir(path1)
					mkpath(path1)
				end
				pathcov = string(path1,"/Cov.csv")
				df_cov = DataFrame(ErrSC = Err_SC, CovSC = IDCov_SC)
				CSV.write(pathcov, df_cov)
			end
		end
	end
end
