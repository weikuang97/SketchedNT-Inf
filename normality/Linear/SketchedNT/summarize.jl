using CSV
using Statistics
using DataFrames
using Plots
using LaTeXStrings


# workdir = "/.../Regression/Plot/Linear/SketchedNT"

include(string(workdir, "/Parameter/Param.jl"))
include(string(workdir, "/StoNewton/StoNewtonMain.jl"))

# load parameter
Max_Iter,Rep,tau = StoNewtonSet.MaxIter,StoNewtonSet.Rep,StoNewtonSet.tau
c_1,c_2,c_3 = StoNewtonSet.c_1,StoNewtonSet.c_2,StoNewtonSet.c_3
Sigma = StoNewtonSet.Sigma
D = StoNewtonSet.D


diff_mat = zeros(div(Max_Iter, 50000), Rep)

path1 = string(workdir,"/Solution/tau",tau)

# go over all repetitions
for IdRep = 1:Rep
    Time = time()
    println("IDRep:",IdRep)
    path = string(path1,"/rep",IdRep,"/Diff.csv")
    df_diff = CSV.File(path; header=true) |> DataFrame
    diff_mat[:,IdRep] = df_diff.diff_vec

    Time = time()-Time
    println("Time:",Time)
end

# save the results
path2 = string(workdir, "/Solution/Figures/tau",tau)
if !isdir(path2)
    mkpath(path2)
end

pathcov = string(path2,"/Diff_mat.csv")
Df_diff = DataFrame(diff_t1 = diff_mat[1,:], diff_t2 = diff_mat[2,:], diff_t3 = diff_mat[3,:], diff_t4 = diff_mat[4,:], diff_t5 = diff_mat[5,:], diff_t6 = diff_mat[6,:])
CSV.write(pathcov, Df_diff)
