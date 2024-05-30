#* For now, I'm going to just focus on parameter estimation. This is an easier problem than random effects prediction
#* Later, I can introduce prediction of random effects

using Distributions, Random, Statistics
using LogExpFunctions  
using DataFrames
using MixedModels
using ProgressMeter, DrWatson
using JLD2
using GLM   # For linear and logistic regression with fixed effects

include("../src/Helpers - General.jl")
include("../src/Helpers - Bin-Resp, Bin-Med, Fixed.jl")

# include("src/Helpers - General.jl")
# include("src/Helpers - Bin-Resp, Bin-Med, Fixed.jl")



Random.seed!(11111)

num_reps = 1000    # Number of datasets to generate

n = 100         # Sample size for each group
K = 5           # Number of groups

p_conf = 3      # Number of confounders


# M model
## Fixed effects
a_1 = 1
A_2 = repeat([1], p_conf)
a_0 = - (a_1 + sum(A_2))    # Offset mean effect of other coefficients in linear predictor for M

## Random effects
sigma_a_0 = 0.2 * abs(a_0)
sigma_a_1 = 0.2 * abs(a_1)



# Y model
## Fixed effects
b_1 = 1
b_2 = 1
B_3 = repeat([1], p_conf)
b_0 = - (b_1 + 0.5*b_2 + sum(B_3))    # Approximately offset mean effect of other coefficients in linear predictor for Y

## Random effects
sigma_b_0 = 0.2 * abs(b_0)
# sigma_b_1 = 0.2 * abs(b_1)
sigma_b_2 = 0.2 * abs(b_2)



# Reference values at which we compute the total mediation effect
x_pred = 0
W_pred = repeat([1], p_conf)    # W=[1,1,1]


# X_dist = Bernoulli(0.5)
# W_dist = Bernoulli(0.5)

X_dist = Normal(1, 1)
W_dist = Normal(1, 1)



all_a_hats = []
all_b_hats = []

all_a_SEs = []
all_b_SEs = []

all_med_hats = []
all_med_SEs = []

# @showprogress Threads.@threads for _ in 1:M
@showprogress for _ in 1:num_reps

    
    # Make X and W
    all_Xs = [rand(X_dist, n) for _ in 1:K]
    all_Ws = [rand(W_dist, n, p_conf) for _ in 1:K]
    
    # Make data
    ## M
    all_eta_vecs_fixed = [a_0 .+ a_1 .* X .+ W * A_2 for (X, W) in zip(all_Xs, all_Ws)]

    ### Add random effects
    all_eta_vecs = [eta_vec_fixed .+ rand(Normal(0, sigma_a_0), n) .+ rand(Normal(0, sigma_a_1), n) .* X for (eta_vec_fixed, X) in zip(all_eta_vecs_fixed, all_Xs)]

    ### Get probabilities
    all_p_M_vecs = [expit.(eta_vec) for eta_vec in all_eta_vecs]

    ### Generate M
    all_Ms = [rand.(Bernoulli.(p_M_vec)) for p_M_vec in all_p_M_vecs]


    ## Y
    all_zeta_vecs_fixed = [b_0 .+ b_1 .* M .+ b_2 .* X .+ W * B_3 for (M, X, W) in zip(all_Ms, all_Xs, all_Ws)]

    ### Add random effects
    all_zeta_vecs = [zeta_vec_fixed .+ rand(Normal(0, sigma_b_0), n) .+ rand(Normal(0, sigma_b_2), n) .* X for (zeta_vec_fixed, X) in zip(all_zeta_vecs_fixed, all_Xs)]
    
    ### Get probabilities
    all_p_Y_vecs = [expit.(zeta_vec) for zeta_vec in all_zeta_vecs]

    ### Generate Y
    all_Ys = [rand.(Bernoulli.(p_Y_vec)) for p_Y_vec in all_p_Y_vecs]


    ## Combine data from all groups
    X = vcat(all_Xs...)
    W = vcat(all_Ws...)
    M = vcat(all_Ms...)
    Y = vcat(all_Ys...)
    group = repeat(1:K, inner=n)




    ## Fit M
    M_data = DataFrame(M=M, X=X, group=group)
    for p in 1:p_conf
        M_data[!,Symbol("W$p")] = W[:, p]
    end

    ### This commented code tries to automate construction of the formula. I ran into issues specifying the random effects
    # M_formula_fix = term(:M) ~ sum(term.(names(M_data, Not(:M, :group))))
    # RE_terms = term(1) + term(:X)
    # group_term = term(:group)
    # M_formula = M_formula_fix + (RE_terms | group_term)

    M_formula = @formula(M ~ X + W1 + W2 + W3 + (1 + X | group))

    M_model = fit(MixedModel, M_formula, M_data, Bernoulli(), verbose=false, progress=false)

    # M_model.InterceptTerm

    a_hat = coef(M_model)
    a_SE = stderror(M_model)
    a_cov = vcov(M_model)

    push!(all_a_hats, a_hat)
    push!(all_a_SEs, a_SE)

    push!(all_a_hats, coef(M_model))
    push!(all_a_SEs, stderror(M_model))

    ## Fit Y
    Y_data = DataFrame(Y=Y, M=M, X=X)
    for p in 1:p_conf
        Y_data[!,Symbol("W$p")] = W[:, p]
    end

    Y_formula = term(:Y) ~ sum(term.(names(Y_data, Not(:Y))))

    Y_model = glm(Y_formula, Y_data, Binomial(), LogitLink())
    b_hat = coef(Y_model)
    b_SE = stderror(Y_model)
    b_cov = vcov(Y_model)

    push!(all_b_hats, b_hat)
    push!(all_b_SEs, b_SE)





    


    ## Get coefficients from models for M and Y
    a_hat = coef(M_model)
    b_hat = coef(Y_model)

    a_x = a_hat[2]
    b_m = b_hat[2]
    b_x = b_hat[3]



    ## Compute linear predictor for M and Y at reference levels of X and W    
    eta_hat = a_hat[1] + a_hat[2] * x_pred + sum(a_hat[3:end] .* W_pred)    # Lin pred for M
    zeta_hat = b_hat[1] + b_hat[3] * x_pred + sum(b_hat[4:end] .* W_pred)   # Lin pred for Y (without M contribution)


    ## Compute mediation effect
    med_hat = get_odds_ratio(eta_hat, a_x, zeta_hat, b_x, b_m)
    push!(all_med_hats, med_hat)




    ## Analytical standard error

    ### Regression coefficients' asymptotic covariance matrix
    a_size = size(a_cov, 1)
    b_size = size(b_cov, 1)
    total_size = a_size + b_size
    reg_covs = zeros(total_size, total_size)
    reg_covs[1:a_size, 1:a_size] = a_cov
    reg_covs[a_size+1:end, a_size+1:end] = b_cov

    reg_asymp_cov = n .* reg_covs


    ### Derivatives of gamma wrt each parameter
    dOR_dt = d_OR_d_theta(eta_hat, a_x, zeta_hat, b_x, b_m, x_pred, W_pred)

    ### Compute asymp var of gamma_hat
    asymp_var = dOR_dt' * reg_asymp_cov * dOR_dt
    asymp_SE = sqrt(asymp_var)

    ### Compute SE of gamma_hat
    med_SE = asymp_SE / sqrt(n)



    push!(all_med_SEs, med_SE)
end



par_list_name = @savename num_reps n p_conf
@save datadir("Bin-Resp, Bin-Med, Fixed-$par_list_name.jld2") all_a_hats all_a_SEs all_b_hats all_b_SEs all_med_hats all_med_SEs
println("Work complete!")




