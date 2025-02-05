using CSV
using Statistics
using DataFrames
using Plots
using LaTeXStrings


workdir = "/.../Regression/Plot/Linear/SketchedNT"

include(string(workdir, "/Parameter/Param.jl"))
include(string(workdir, "/StoNewton/StoNewtonMain.jl"))

# load parameter
Max_Iter,Rep,tau = StoNewtonSet.MaxIter,StoNewtonSet.Rep,StoNewtonSet.tau
c_1,c_2,c_3 = StoNewtonSet.c_1,StoNewtonSet.c_2,StoNewtonSet.c_3
Sigma = StoNewtonSet.Sigma
D = StoNewtonSet.D


Buff= Int(1e5)
# initialize
IDCov_Xistar_ave = zeros(Max_Iter-Buff) # coverage rate of oracle confidence interval
Err_PI_ave = zeros(Max_Iter-Buff) # relative covariance estimation error (plug-in)
IDCov_PI_ave = zeros(Max_Iter-Buff)# coverage rate of oracle confidence interval (plug-in)
Err_SC_ave = zeros(Max_Iter-Buff) # relative covariance estimation error (weighted sample cov)
IDCov_SC_ave = zeros(Max_Iter-Buff) #  coverage rate of oracle confidence interval (weighted sample cov)



path1 = string(workdir,"/Solution/tau",tau)

# go over all repetitions
for IdRep = 1:Rep
    Time = time()
    println("IDRep:",IdRep)
    path = string(path1,"/rep",IdRep,"/Cov.csv")
    df_cov = CSV.File(path; header=true) |> DataFrame

    global IDCov_Xistar_ave = ((IdRep-1)/IdRep).*IDCov_Xistar_ave + (1/IdRep).*df_cov.CovXistar

    global Err_PI_ave = ((IdRep-1)/IdRep).*Err_PI_ave + (1/IdRep).*df_cov.ErrPI
    global IDCov_PI_ave = ((IdRep-1)/IdRep).*IDCov_PI_ave + (1/IdRep).*df_cov.CovPI

    global Err_SC_ave = ((IdRep-1)/IdRep).*Err_SC_ave + (1/IdRep).*df_cov.ErrSC
    global IDCov_SC_ave = ((IdRep-1)/IdRep).*IDCov_SC_ave + (1/IdRep).*df_cov.CovSC

    Time = time()-Time
    println("Time:",Time)
end

# save the results
path2 = string(workdir, "/Solution/Figures/tau",tau)
if !isdir(path2)
    mkpath(path2)
end

pathcov = string(path2,"/Cov_ave.csv")
Df_cov = DataFrame(CovXistar = IDCov_Xistar_ave, ErrPI = Err_PI_ave, CovPI = IDCov_PI_ave, ErrSC = Err_SC_ave, CovSC = IDCov_SC_ave)
CSV.write(pathcov, Df_cov)
