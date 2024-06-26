#* For now, I'm going to just focus on parameter estimation. This is an easier problem than random effects prediction
#* Later, I can introduce prediction of random effects

using Distributions, Random, Statistics
using LogExpFunctions   # For the logit and expit functions
using DataFrames
using MixedModels
using ProgressMeter, DrWatson
using JLD2
using GLM   # For linear and logistic regression with fixed effects

# include("src/Helpers - General.jl")
# include("src/Helpers - Cont-Resp, Bin-Med, Fixed.jl")

include("../src/Helpers - General.jl")
include("../src/Helpers - Cont-Resp, Bin-Med, Fixed.jl")




Random.seed!(1)

num_reps = 1000    # Number of datasets to generate

n = 1000
p_conf = 3      # Number of confounders

a_1 = 1
A_2 = repeat([1], p_conf)
a_0 = - (a_1 + sum(A_2))    # Offset mean effect of other coefficients in linear predictor for M

b_1 = 1
b_2 = 1
B_3 = repeat([1], p_conf)
b_0 = - (b_1 + 0.5*b_2 + sum(B_3))    # Approximately offset mean effect of other coefficients in linear predictor for Y
sigma_e = 0.2   # Residual SD for Y



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

# @showprogress Threads.@threads for _ in 1:num_reps
@showprogress for _ in 1:num_reps

    
    # Make X and W
    X = rand(X_dist, n)
    W = rand(W_dist, n, p_conf)
    
    # Make data
    ## M
    eta_M_vec = a_0 .+ a_1 .* X .+ W * A_2
    p_M_vec = expit.(eta_M_vec)
    M = rand.(Bernoulli.(p_M_vec))


    ## Y
    epsilon = rand(Normal(0, sigma_e), n)
    Y = b_0 .+ b_1 .* M .+ b_2 .* X .+ W * B_3 .+ epsilon




    ## Fit M
    M_data = DataFrame(M=M, X=X)
    for p in 1:p_conf
        M_data[!,Symbol("W$p")] = W[:, p]
    end

    M_formula = term(:M) ~ sum(term.(names(M_data, Not(:M))))

    M_model = glm(M_formula, M_data, Binomial(), LogitLink())
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

    Y_model = lm(Y_formula, Y_data)
    b_hat = coef(Y_model)
    b_SE = stderror(Y_model)
    b_cov = vcov(Y_model)

    push!(all_b_hats, b_hat)
    push!(all_b_SEs, b_SE)





    ## Compute linear predictor for M at reference levels of X and W    
    eta_hat = a_hat[1] + a_hat[2] * x_pred + sum(a_hat[3:end] .* W_pred)



    ## Get relevant coefficients from models for M and Y
    a_x = a_hat[2]
    b_m = b_hat[2]
    b_x = b_hat[3]


    ## Compute mediation effect
    delta = get_delta(eta_hat, a_x)
    med_hat = get_gamma(delta, b_m, b_x)

    push!(all_med_hats, med_hat)

    ## Standard error

    ### Regression coefficients' asymptotic covariance matrix
    a_size = size(a_cov, 1)
    b_size = size(b_cov, 1)
    total_size = a_size + b_size
    reg_covs = zeros(total_size, total_size)
    reg_covs[1:a_size, 1:a_size] = a_cov
    reg_covs[a_size+1:end, a_size+1:end] = b_cov

    reg_asymp_cov = n .* reg_covs


    ### Derivatives of gamma wrt each parameter
    dg_dt = d_gamma_d_theta(eta_hat, x_pred, W_pred, a_hat, b_hat)

    ### Compute asymp var of gamma_hat
    asymp_var = dg_dt' * reg_asymp_cov * dg_dt
    asymp_SE = sqrt(asymp_var)

    ### Compute SE of gamma_hat
    med_SE = asymp_SE / sqrt(n)



    push!(all_med_SEs, med_SE)
end



par_list_name = @savename num_reps n p_conf
@save datadir("Cont-Resp, Bin-Med, Fixed-$par_list_name.jld2") all_a_hats all_a_SEs all_b_hats all_b_SEs all_med_hats all_med_SEs
println("Work complete!")




