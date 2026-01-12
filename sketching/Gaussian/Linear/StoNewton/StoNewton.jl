function StoNewton(c1,c2,c3,Max_Iter,ttau,nx,X_true,Sigma,Xistar,sigma,q)
    ## Set Buffer Size
    Buff= Int(1e5)

    ## Initialize variables
    t,X_t,NewDir_t = 1,zeros(nx),zeros(nx) # iterates

    cum_bar2x_f_t,B_t = Matrix(1.0I,nx,nx),zeros(nx,nx) # sample cov estimator
    W_t,V_t,u_t,barx_t = zeros(nx,nx),zeros(nx),0,zeros(nx)
    Xi_t_SC = zeros(nx,nx)

    alpha_t = 0 # stepsize

    ## Start the loop
    Time = time()
    while t <= Max_Iter
        # Step 0: compute deterministic quantities
        B_t = cum_bar2x_f_t

        beta_t = c1/t^c2 # stepsize
        chi_t = beta_t^c3
        phi_t = beta_t+chi_t/2

        # Step 1: generate random gradient and Hessian
        a_t, eps_t = rand(MvNormal(zeros(nx),Sigma)), sigma*rand(Normal(0,1))
        barg_t = a_t*(a_t'*(X_t-X_true)) - eps_t*a_t
        bar_nab_x2f_t = a_t*a_t'
        cum_bar2x_f_t = t/(t+1)*cum_bar2x_f_t + 1/(t+1)*bar_nab_x2f_t


        # Step 2: solve Newton system via Randomized Solve
        if ttau == 0
            NewDir_t = lu(B_t)\-barg_t
        else
            NewDir_t = zeros(nx)
            for inner_iter = 1:ttau
                # generate sketching matrix
                S_t = randn(nx, q)
                NewDir_t = NewDir_t - (B_t*S_t) * pinv((B_t*S_t)' * (B_t*S_t)) * S_t' * (B_t * NewDir_t + barg_t)
            end
        end


        # Step 3: update the iterate
        alpha_t = rand(Uniform(beta_t,beta_t+chi_t))
        X_t = X_t + alpha_t*NewDir_t


        # Step 4: covariance estimation
        if t > Buff
            # update weighted sample covariance matrix
            W_t = (t-Buff-1)/(t-Buff)*W_t + 1/(t-Buff)*X_t*X_t'/phi_t
            V_t = (t-Buff-1)/(t-Buff)*V_t + 1/(t-Buff)*X_t/phi_t
            u_t = (t-Buff-1)/(t-Buff)*u_t + 1/(t-Buff)*(1/phi_t)
            barx_t = (t-Buff-1)/(t-Buff)*barx_t + 1/(t-Buff)*X_t
        end

        t = t + 1
    end
    Time = time() - Time

    # compute weighted sample covariance matrix
    Xi_t_SC = W_t-V_t*barx_t'-barx_t*V_t'+u_t*barx_t*barx_t'
    Err_SC = opnorm(Xi_t_SC-Xistar, 2)/opnorm(Xistar,2)


    # statistical inference
    COV_value_SC = sum(Xi_t_SC)/(nx)^2
    if COV_value_SC<1e-6
        COV_value_SC = 0
    end
    radius_SC = 1.96*sqrt(alpha_t)*sqrt(COV_value_SC)
    width_SC = 2*radius_SC


    Xmean = mean(X_t)
    IdCov_SC = (mean(X_true)>=(Xmean-radius_SC)) && (mean(X_true)<=(Xmean+radius_SC))

    return Time,Xistar,IdCov_SC,COV_value_SC,Err_SC,width_SC
end
