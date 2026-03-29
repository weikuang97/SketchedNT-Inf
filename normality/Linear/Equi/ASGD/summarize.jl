using CSV
using Statistics
using DataFrames
using Plots
using LaTeXStrings
using LinearAlgebra


workdir = "/.../normality/Linear/Toep/ASGD"

include(string(workdir, "/Parameter/Param.jl"))
include(string(workdir, "/ASGD/ASGDMain.jl"))

# load parameters
Max_Iter,Rep = ASGDSet.MaxIter,ASGDSet.Rep
c_1,c_2 = ASGDSet.c_1,ASGDSet.c_2
Sigma = ASGDSet.Sigma
D = ASGDSet.D

# initialize
diff_vec = zeros(Rep)

path1 = string(workdir,"/Solution")
# go over all repetitions
for IdRep = 1:Rep
    Time = time()
    println("IDRep:",IdRep)
    path = string(path1,"/rep",IdRep,"/Diff.csv")
    df_diff = CSV.File(path; header=true) |> DataFrame
    diff_vec[IdRep] = df_diff.diff[1]

    Time = time()-Time
    println("Time:",Time)
end


# save the results
path2 = string(workdir, "/Solution/Figures")
if !isdir(path2)
    mkpath(path2)
end

pathcov = string(path2,"/Diff_mat.csv")
Df_diff = DataFrame(diff_vec = diff_vec)
CSV.write(pathcov, Df_diff)
