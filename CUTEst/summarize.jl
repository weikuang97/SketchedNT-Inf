using CSV
using Statistics
using DataFrames
using Plots
using JLD2

workdir = "/.../CUTEst"

include(string(workdir, "/Parameter/Param.jl"))
include(string(workdir, "/StoSQP/StoSQPMain.jl"))

# load parameter
Max_Iter,Rep,tau = StoSQPSet.MaxIter,StoSQPSet.Rep,StoSQPSet.tau
c_1,c_2,c_3 = StoSQPSet.c_1,StoSQPSet.c_2,StoSQPSet.c_3
sigma2 = StoSQPSet.sigma2
Prob = StoSQPSet.Prob


struct Summarize
    IDCov_SC_ave::Vector{Float64} # coverage rate of confidence intervals
    Vardiff_SC::Vector{Float64} # variance estimation error
end



Buff= Int(1e5)

IDCov_SC_ave = zeros(Max_Iter-Buff)
COV_value_SC_ave = zeros(Max_Iter-Buff)
COV_value_Oracle = 0



for Idsigma2 = 1:length(sigma2)
    for Idtau = 1:length(tau)
        path1 = string(workdir, "/Solution/sigma", Idsigma2, "/tau",tau[Idtau])

        IDCov_SC_ave = zeros(Max_Iter-Buff)
        COV_value_SC_ave = zeros(Max_Iter-Buff)
        COV_value_Oracle = 0

        Vardiff_PI = zeros(Max_Iter-Buff)
        Vardiff_SC = zeros(Max_Iter-Buff)

        for IdRep = 1:Rep
            println("Rep:",IdRep)
            path = string(path1,"/rep",IdRep,".jld2")
            @load path Result
            Xistar = Result.Xistar
            Index = Int.(Result.Index)
            COV_value_Oracle = sum(Xistar[Index,Index])/(length(Index))^2

            global IDCov_SC_ave = ((IdRep-1)/IdRep).*IDCov_SC_ave + (1/IdRep).*Result.IDCov_SC

            global COV_value_SC_ave = ((IdRep-1)/IdRep).*COV_value_SC_ave + (1/IdRep).*Result.COV_value_SC
        end
        Vardiff_SC .= (COV_value_SC_ave.-COV_value_Oracle)/COV_value_Oracle

        path2 = string(workdir, "/Solution/Figures/sigma", Idsigma2, "/tau",tau[Idtau])
        if !isdir(path2)
            mkpath(path2)
        end

        pathcov = string(path2, "/Cov.csv")
        Df_COV = DataFrame(IDCov_SC = IDCov_SC_ave, Vardiff_SC = Vardiff_SC)
        CSV.write(pathcov, Df_COV)
    end
end
