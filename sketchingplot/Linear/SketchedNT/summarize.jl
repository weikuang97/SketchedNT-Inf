using CSV
using Statistics
using DataFrames
using Plots
using LaTeXStrings


# workdir = "/.../sketchingplot/Linear/SketchedNT"
include(string(workdir, "/Parameter/Param.jl"))
include(string(workdir, "/StoNewton/StoNewtonMain.jl"))

# load parameter
Max_Iter,Rep,tau = StoNewtonSet.MaxIter,StoNewtonSet.Rep,StoNewtonSet.tau
c_1,c_2,c_3 = StoNewtonSet.c_1,StoNewtonSet.c_2,StoNewtonSet.c_3
Sigma = StoNewtonSet.Sigma
D = StoNewtonSet.D


Buff= Int(1e5)
# initialize
Err_SC_ave = zeros(Max_Iter-Buff) # relative covariance estimation error (weighted sample cov)
IDCov_SC_ave = zeros(Max_Iter-Buff) #  coverage rate of oracle confidence interval (weighted sample cov)




path1 = string(workdir,"/Solution")

# go over all cases
q_vec = [1]

for q in q_vec
    for ttau in tau
        for IdRep = 1:Rep
            Time = time()
            println("q-", q, "-tau-", ttau, "-Rep-", IdRep)
            path = string(path1,"/q",q,"/tau",ttau,"/rep",IdRep,"/Cov.csv")
            df_cov = CSV.File(path; header=true) |> DataFrame

            global Err_SC_ave = ((IdRep-1)/IdRep).*Err_SC_ave + (1/IdRep).*df_cov.ErrSC
            global IDCov_SC_ave = ((IdRep-1)/IdRep).*IDCov_SC_ave + (1/IdRep).*df_cov.CovSC

            Time = time()-Time
            println("Time:",Time)
        end

        # save the results
        path2 = string(workdir, "/Solution/Figures/q",q,"/tau",ttau)
        if !isdir(path2)
            mkpath(path2)
        end

        pathcov = string(path2,"/Cov_ave.csv")
        Df_cov = DataFrame(ErrSC = Err_SC_ave, CovSC = IDCov_SC_ave)
        CSV.write(pathcov, Df_cov)
    end
end
