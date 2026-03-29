using CSV
using Statistics
using DataFrames
using Plots
using LaTeXStrings
using LinearAlgebra


workdir = "/.../normality/Linear/Toep/SGD"

include(string(workdir, "/Parameter/Param.jl"))
include(string(workdir, "/SGD/SGDMain.jl"))

# load parameters
Max_Iter,Rep = SGDSet.MaxIter,SGDSet.Rep
c_1,c_2 = SGDSet.c_1,SGDSet.c_2
Sigma = SGDSet.Sigma
D = SGDSet.D

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
