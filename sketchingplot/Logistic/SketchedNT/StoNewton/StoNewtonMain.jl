include("trueCov.jl")
include("StoNewton.jl")
## Implement sketched Newton inference framework for whole problem set
# StoNewtonSet: parameters

function StoNewtonMain(StoNewtonSet)
	# load parameters
	Max_Iter,Rep,tau = StoNewtonSet.MaxIter,StoNewtonSet.Rep,StoNewtonSet.tau
	c_1,c_2,c_3 = StoNewtonSet.c_1,StoNewtonSet.c_2,StoNewtonSet.c_3
	Sigma,mu = StoNewtonSet.Sigma,StoNewtonSet.mu
	D = StoNewtonSet.D


	nx = D
	# true parameters and true covariance matrix
	X_true = (1/nx)*ones(nx)
	Bstar, Omegastar = Bstar_Omegastar(Sigma, X_true, nx, mu)
	Xistar = trueCov(Omegastar, Bstar, nx, tau)

	q_vec = [1,2,3,4,5]


	# go over all repetitions
	for q in q_vec:
		for ttau in tau:
			for IdRep = 1:Rep
				println("Rep:", IdRep)
				Time,IDCov_Xistar,Err_PI,IDCov_PI,Err_SC,IDCov_SC = StoNewton(c_1,c_2,c_3,Max_Iter,ttau,nx,X_true,Sigma,Xistar,q,mu)
		        println("Time:", Time)

				path1 = string("../Solution/q",q,"/ttau",tau,"/rep",IdRep)
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
