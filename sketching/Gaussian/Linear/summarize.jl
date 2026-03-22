using JLD2
using Statistics

workdir = "/.../sketching/Gaussian/Linear"


include(string(workdir, "/Parameter/Param.jl"))
include(string(workdir, "/StoNewton/StoNewtonMain.jl"))

# load parameters
R_Tail=StoNewtonSet.MaxIter; Rep=StoNewtonSet.Rep; tau=StoNewtonSet.tau;
c_1=StoNewtonSet.c_1; c_2=StoNewtonSet.c_2; c_3=StoNewtonSet.c_3;
SigToe=StoNewtonSet.SigToe; SigEqui=StoNewtonSet.SigEqui;
D=StoNewtonSet.D;
LenTau=length(tau); LenDim=length(D);
LenSigToe=length(SigToe);LenSigEqui = length(SigEqui)

# initialize
Err_SC = zeros(Rep) # error in variance estimation of parameter mean (weighted sample cov)
RateCov_SC = zeros(Rep) # coverage rate of confidence interval (weighted sample cov)
VarErr_SC = zeros(Rep)
width_SC = zeros(Rep)


# go over all cases
for IdDim=1:LenDim, IdTau=1:LenTau
    path1 = string(workdir,"/Solution/Dim",IdDim,"Tau",IdTau)
    nx = D[IdDim]
    q_vec = [1,5,10,20]

    # Toeplitz Covariance
    for IdSigToe = 1:LenSigToe
        for q in q_vec
            # reset
            Err_PI = zeros(Rep)
            Err_SC = zeros(Rep)
            RateCov_PI = zeros(Rep)
            RateCov_SC = zeros(Rep)


            # load experiment results
            for IdRep = 1:Rep
                println("IDRep:",IdRep)
                path = string(path1,"/Toe",IdSigToe,"/q",q,"/rep",IdRep,".jld2")

                @load path Result
                COV_value_orc = mean(Result.trueCov)
                VarErr_SC[IdRep] = (Result.COV_value_SC - COV_value_orc)/COV_value_orc
                RateCov_SC[IdRep] = Result.IdCov_SC
                Err_SC[IdRep] = Result.Err_SC
                width_SC[IdRep] = Result.width_SC
            end

            # summarize results
            path2 = string(workdir, "/Figures")
            if !isdir(path2)
                mkpath(path2)
            end

            io = open(string(path2, "/Error.txt"), "a");
            write(io, string("Dim-",D[IdDim],"-Tau-",tau[IdTau],"-Toe-",SigToe[IdSigToe],"-q-",q,"\n"))
            write(io, string("Err_SC:", mean(Err_SC),'\n'))
            write(io, string("VarErr_SC:", mean(VarErr_SC),'\n'))
            close(io)


            io = open(string(path2, "/Rate.txt"), "a");
            write(io, string("Dim-",D[IdDim],"-Tau-",tau[IdTau],"-Toe-",SigToe[IdSigToe],"-q-",q,"\n"))
            write(io, string("CovRate_SC:", mean(RateCov_SC),'\n'))
            write(io, string("width_SC:", mean(width_SC),'\n'))
            close(io)
        end
    end


    # Equi-corr Covariance
    for IdSigEqui = 1:LenSigEqui
        for q in q_vec
            # reset
            Err_PI = zeros(Rep)
            Err_SC = zeros(Rep)
            RateCov_PI = zeros(Rep)
            RateCov_SC = zeros(Rep)


            # load experiment results
            for IdRep = 1:Rep
                println("IDRep:",IdRep)
                path = string(path1,"/Equ",IdSigEqui,"/q",q,"/rep",IdRep,".jld2")

                @load path Result
                COV_value_orc = mean(Result.trueCov)
                VarErr_SC[IdRep] = (Result.COV_value_SC - COV_value_orc)/COV_value_orc
                RateCov_SC[IdRep] = Result.IdCov_SC
                Err_SC[IdRep] = Result.Err_SC
                width_SC[IdRep] = Result.width_SC
            end

            # summarize results
            path2 = string(workdir, "/Figures")
            if !isdir(path2)
                mkpath(path2)
            end

            io = open(string(path2, "/Error.txt"), "a");
            write(io, string("Dim-",D[IdDim],"-Tau-",tau[IdTau],"-Equ-",SigEqui[IdSigEqui],"-q-",q,"\n"))
            write(io, string("Err_SC:", mean(Err_SC),'\n'))
            write(io, string("VarErr_SC:", mean(VarErr_SC),'\n'))
            close(io)


            io = open(string(path2, "/Rate.txt"), "a");
            write(io, string("Dim-",D[IdDim],"-Tau-",tau[IdTau],"-Equ-",SigEqui[IdSigEqui],"-q-",q,"\n"))
            write(io, string("CovRate_SC:", mean(RateCov_SC),'\n'))
            write(io, string("width_SC:", mean(width_SC),'\n'))
            close(io)
        end
    end
end
