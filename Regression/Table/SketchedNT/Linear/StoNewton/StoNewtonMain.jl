include("trueCov.jl")
include("StoNewton.jl")
struct StoNewtonResult
	Time::Float64					# running time
	trueCov::Matrix{Float64}		# true limiting covariance matrix
	IdCov_PI::Int64					# 1 or 0, true parameter covered by confidence interval or not (plug-in)
	COV_value_PI::Float64			# variance estimator of parameter mean (plug-in)
	IdCov_SC::Int64					# 1 or 0, true parameter covered by confidence interval or not (weighted sample cov)
	COV_value_SC::Float64			# variance estimator of parameter mean (weighted sample cov)
end


## Implement sketched Newton framework for whole problem set
# StoNewtonSet: parameters of Newton inference

function StoNewtonMain(StoNewtonSet)
	# load parameters
	Max_Iter,Rep,tau = StoNewtonSet.MaxIter,StoNewtonSet.Rep,StoNewtonSet.tau
	c_1,c_2,c_3 = StoNewtonSet.c_1,StoNewtonSet.c_2,StoNewtonSet.c_3
	SigToe,SigEqui = StoNewtonSet.SigToe,StoNewtonSet.SigEqui
	D = StoNewtonSet.D
	sigma = StoNewtonSet.sigma
	LenTau,LenDim = length(tau),length(D)
	LenSigToe,LenSigEqui = length(SigToe),length(SigEqui)

	# Go over all cases
	for IdDim=1:LenDim, IdTau=1:LenTau
		nx = D[IdDim]
		X_true = (1/nx)*ones(nx) # true parameter

		# Toeplitz Covariance
		for IdSigToe = 1:LenSigToe
			Sigma = [(SigToe[IdSigToe])^abs(i-j) for i=1:nx, j=1:nx]
			Bstar = Sigma
			Xistar = trueCov(Bstar, nx, tau[IdTau], sigma) # true limiting covariance

			for IdRep = 1:Rep
			    println("D-",IdDim,"-Tau-",IdTau,"-CovT-",IdSigToe,"-Rep-",IdRep)

			    Time,trueCov,IdCov_PI,COV_value_PI,IdCov_SC,COV_value_SC =
			    StoNewton(c_1,c_2,c_3,Max_Iter,tau[IdTau],nx,X_true,Sigma,Xistar,sigma)
			    println("Time,", Time)

				# save results
				path1 = string("../Solution/Dim",IdDim,"Tau",IdTau,"/Toe",IdSigToe)
				if !isdir(path1)
				    mkpath(path1)
				end
				path = string(path1,"/rep",IdRep,".jld2")
				Result = StoNewtonResult(Time,trueCov,Int.(IdCov_PI),COV_value_PI,Int.(IdCov_SC),COV_value_SC)
				@save path Result
			end
		end

		# Equi Covariance
		for IdSigEqui = 1:LenSigEqui
			Sigma = SigEqui[IdSigEqui]*ones(nx,nx)+(1-SigEqui[IdSigEqui])*Matrix(1.0I,nx,nx)
			Bstar = Sigma
			Xistar = trueCov(Bstar, nx, tau[IdTau], sigma) # true limiting covariance

			for IdRep = 1:Rep
			    println("D-",IdDim,"-Tau-",IdTau,"-CovE-",IdSigEqui,"-Rep-",IdRep)

			    Time,trueCov,IdCov_PI,COV_value_PI,IdCov_SC,COV_value_SC =
			    StoNewton(c_1,c_2,c_3,Max_Iter,tau[IdTau],nx,X_true,Sigma,Xistar,sigma)
			    println("Time,", Time)

				# save results
				path1 = string("../Solution/Dim",IdDim,"Tau",IdTau,"/Equ",IdSigEqui)
				if !isdir(path1)
				    mkpath(path1)
			    end
				path = string(path1,"/rep",IdRep,".jld2")
				Result = StoNewtonResult(Time,trueCov,Int.(IdCov_PI),COV_value_PI,Int.(IdCov_SC),COV_value_SC)
				@save path Result
			end
		end
	end

end
