# true limiting covariance
function trueCov(Kstar, CovM, nx, nlam, tau)
    d = nx+nlam
    Omegastar = inv(Matrix(Kstar))*hcat(vcat(CovM, zeros(nlam,nx)), zeros(d,nlam))*inv(Matrix(Kstar))
    if (tau==0)
        Xistar = Omegastar/2
    else
        Kstar_norm = mapslices(normalize, Kstar, dims=1)
        Dstar_j = (Kstar_norm*Kstar_norm')./d
        Cstar_j = Matrix(1.0I,d,d)-Dstar_j
        Cstar = (Cstar_j^tau).*(-1)

        eigen_result = eigen(Hermitian(Matrix(1.0I,d,d)+Cstar))
        sigma = eigen_result.values
        U = eigen_result.vectors
        Theta = inv.(sigma*ones(d)'+ones(d)*sigma')

        mid_temp = Omegastar
        for i = 1:tau
            mid_temp = onestep(mid_temp, Kstar_norm, Dstar_j, d)
        end
        exp_term = Omegastar + Cstar*Omegastar + Omegastar*Cstar' + mid_temp

        Xistar = U*(Theta.*(U'*exp_term*U))*U'
    end

    return Xistar
end


function onestep(Omegastar, Kstar_norm, Dstar_j, d)
    wt = mapslices(vec -> weight(vec, Omegastar), Kstar_norm, dims=1)[1,:]
    term = (Kstar_norm*Matrix(Diagonal(wt))*Kstar_norm')./d
    res = Omegastar - Dstar_j*Omegastar - Omegastar*Dstar_j' + term
    return res
end

function weight(vec, Omegastar)
    return vec'*Omegastar*vec
end
