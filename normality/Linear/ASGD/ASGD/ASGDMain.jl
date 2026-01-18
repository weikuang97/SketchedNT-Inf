include("ASGD.jl")


## Implement SGD framework for whole problem set
# ASGDSet: parameters of ASGD inference
function ASGDMain(ASGDSet)
	# load parameters
	Max_Iter,Rep = ASGDSet.MaxIter,ASGDSet.Rep
	c_1,c_2 = ASGDSet.c_1,ASGDSet.c_2
	Sigma,sigma = ASGDSet.Sigma,ASGDSet.sigma
	D = ASGDSet.D


	nx = D
	# true parameters and true covariance matrix
    X_true = (1/nx)*ones(nx)
	Xistar = (sigma^2)*inv(Sigma)

	# go over all repetitions
	for IdRep = 1:Rep
		println("Rep:", IdRep)
		Time,diff_vec = ASGD(c_1,c_2,Max_Iter,nx,X_true,Sigma,Xistar,sigma)
        println("Time:", Time)

		path1 = string("../Solution/rep",IdRep)
		if !isdir(path1)
			mkpath(path1)
		end
		pathcov = string(path1,"/Diff.csv")
		df_diff = DataFrame(diff_vec = diff_vec)
		CSV.write(pathcov, df_diff)
	end
end
