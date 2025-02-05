include("trueCov.jl")
include("StoSQP.jl")
## Implement sketched Newton inference framework for whole problem set
# StoNewtonSet: parameters

struct StoSQPResult
	Time::Float64				       # running time
	Xistar::Matrix{Float64}			   # true limiting covariance matrix
	Index::Vector{Int64}			   # Inference index
	IDCov_SC::Vector{Int64}			   # 1 or 0, true parameter covered by confidence interval or not
	COV_value_SC::Vector{Float64}	   # variance estimation
	IdSing::Int64					   # 1 or 0, iterate converges or not
end


function StoSQPMain(StoNewtonSet)
	# load parameter
	Max_Iter,Rep,tau = StoNewtonSet.MaxIter,StoNewtonSet.Rep,StoNewtonSet.tau
	c_1,c_2,c_3 = StoNewtonSet.c_1,StoNewtonSet.c_2,StoNewtonSet.c_3
	sigma2 = StoNewtonSet.sigma2
	Prob = StoNewtonSet.Prob


	nlp = CUTEstModel(Prob)
	stats = ipopt(nlp,print_level=0)
	X_true, Lam_true = stats.solution, stats.multipliers
	nx, nlam = nlp.meta.nvar, nlp.meta.ncon

	if stats.solver_specific[:internal_msg] == :Solve_Succeeded
		# go over all cases
		for Idsigma2 = 1:length(sigma2)
			# go over all tau
			for Idtau = 1:length(tau)
				# count for Rep convergence runs
				IdRep = 1
			    while(IdRep <= Rep)
				    println("Rep:", IdRep)
				    Time,Xistar,Index,IDCov_SC,COV_value_SC,IdSing = StoSQP(nlp,c_1,c_2,c_3,sigma2[Idsigma2],Max_Iter,tau[Idtau],X_true,Lam_true)
		            println("Time:", Time)
				    if IdSing == 1
						println("Rep:",IdRep," not converge")
					else
					    path1 = string("../Solution/sigma",Idsigma2,"/tau",tau[Idtau],)
					    if !isdir(path1)
						    mkpath(path1)
				     	end

						path = string(path1,"/rep",IdRep,".jld2")
						Result = StoSQPResult(Time,Xistar,Index,IDCov_SC,COV_value_SC,IdSing)
						@save path Result
						IdRep = IdRep + 1
					end
			    end
			end
	    end
    end

	finalize(nlp)

end
