function StoNewton(c1,c2,c3,Max_Iter,ttau,nx,X_true,Sigma,Xistar,sigma)
    ## Initialize variables
    t,X_t,NewDir_t = 1,zeros(nx),zeros(nx) # iterates
    cum_bar2x_f_t,B_t = Matrix(1.0I,nx,nx),zeros(nx,nx)
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
                # pick a random index
                j = sample(1:nx)
                NewDir_t = NewDir_t - (B_t[j,:]'*NewDir_t+barg_t[j])/(B_t[j,:]'*B_t[:,j])*B_t[:,j]
            end
        end


        # Step 3: update the iterate
        alpha_t = rand(Uniform(beta_t,beta_t+chi_t))
        X_t = X_t + alpha_t*NewDir_t

        t = t + 1
    end

    diff = sqrt(1/alpha_t) * sum(X_t-X_true)

    Time = time() - Time
    return Time,diff
end
