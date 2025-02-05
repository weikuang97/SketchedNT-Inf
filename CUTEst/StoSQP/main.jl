using NLPModels
using JuMP
using LinearOperators
using OptimizationProblems
using MathProgBase
using ForwardDiff
using CUTEst
using NLPModelsJuMP
using LinearAlgebra
using Distributed
using Ipopt
using DataFrames
using Glob
using DelimitedFiles
using Random
using Distributions
using NLPModelsIpopt
using JLD2



cd("/.../CUTEst/StoSQP")



#using Main.Parameter
include("../Parameter/Param.jl")

include("StoSQPMain.jl")
#######################################
#########  run main file    ###########
#######################################
function main()
    StoSQPMain(StoSQPSet)
end

main()
