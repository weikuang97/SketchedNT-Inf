include("StoNewton.jl")
## Implement sketched Newton inference framework for whole problem set
# StoNewtonSet: parameters

function StoNewtonMain(StoNewtonSet)
	# load parameters
	Max_Iter,Rep,tau = StoNewtonSet.MaxIter,StoNewtonSet.Rep,StoNewtonSet.tau
	c_1,c_2,c_3 = StoNewtonSet.c_1,StoNewtonSet.c_2,StoNewtonSet.c_3
	Sigma,sigma = StoNewtonSet.Sigma,StoNewtonSet.sigma
	D = StoNewtonSet.D


	nx = D
	# true parameters and true covariance matrix
    X_true = (1/nx)*ones(nx)

	# true limiting cov
	if c_2 == 1
		Xistar = inv(Sigma)
	else
		Xistar = 0.5 * inv(Sigma)
	end

    # go over all repetitions
	for IdRep = 1:Rep
		println("Rep:", IdRep)
		Time,diff_vec = StoNewton(c_1,c_2,c_3,Max_Iter,tau,nx,X_true,Sigma,Xistar,sigma)
        println("Time:", Time)

		path1 = string("../Solution/tau",tau,"/rep",IdRep)
		if !isdir(path1)
			mkpath(path1)
		end
		pathcov = string(path1,"/Diff.csv")
		df_diff = DataFrame(diff_vec = diff_vec)
		CSV.write(pathcov, df_diff)
	end
end
