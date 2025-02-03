

j = 1

n = N_K_combs[j, "N"]
K = N_K_combs[j, "K"]

# Reference values for X and W, the confounders
## Note: We compute the mediation effect at X_pred vs X_pred + 1, so X_pred=0 is a good choice
X_pred = 0
W_pred = c(1,1,1)



# Use tryCatch to omit any datasets leading to convergence issues
tryCatch({
  ############## Un-quote
  all_Xs = list()
  all_Ws = list()

  for (k in 1:K) {
    X = rnorm(n, mean = 0, sd = 1)
    W = matrix(rnorm(n * p_conf, mean = 0, sd = 1),
               nrow = n,
               ncol = p_conf)

    all_Xs[[k]] = X
    all_Ws[[k]] = W
  }

  # Generate M
  all_Ms = list()
  for (k in 1:K) {
    eta_vec_fixed = a_0 + a_1 * all_Xs[[k]] + all_Ws[[k]] %*% A_2

    ## Add random effects
    a_ran = mvrnorm(1, mu = rep(0, 2), Sigma = Sigma_a)
    eta_vec = eta_vec_fixed + a_ran[1] + a_ran[2] * all_Xs[[k]]

    ## Generate M
    p_M_vec = expit(eta_vec)
    M = rbinom(n, size = 1, prob = p_M_vec)
    all_Ms[[k]] = M
  }

  # Generate Y
  all_Ys = list()
  for (k in 1:K) {
    zeta_vec_fixed = b_0 + b_1 * all_Ms[[k]] + b_2 * all_Xs[[k]] + all_Ws[[k]] %*%
      B_3

    ## Add random effects
    b_ran = mvrnorm(1, mu = rep(0, 2), Sigma = Sigma_b)
    zeta_vec = zeta_vec_fixed + b_ran[1] + b_ran[2] * all_Xs[[k]]

    ## Generate Y
    p_Y_vec = expit(zeta_vec)
    Y = rbinom(n, size = 1, prob = p_Y_vec)
    all_Ys[[k]] = Y
  }


  # Consolidate groups
  X = do.call(c, all_Xs)
  W = do.call(rbind, all_Ws)
  M = do.call(c, all_Ms)
  Y = do.call(c, all_Ys)
  group = rep(1:K, each = n)


  # Fit models

  ## M
  M_data = data.frame(
    M,
    X,
    W1 = W[, 1],
    W2 = W[, 2],
    W3 = W[, 3],
    group = group
  )
  M_model = glmer(M ~ X + W1 + W2 + W3 + (1 + X |
                                            group),
                  data = M_data,
                  family = binomial(link = "logit"))
  M_model_info = attributes(VarCorr(M_model)$group)

  ### Fitted parameters
  a_hat = fixef(M_model)
  a_RE_sds = M_model_info$stddev
  a_RE_cor = M_model_info$correlation[2, 1]
  theta_hat = c(a_RE_sds, a_RE_cor)
  if (any(is.nan(theta_hat)))
    stop("NaNs in theta_hat")  # Skip rest of current analysis if correlation is 0/0

  ### Estimated SE covariance matrix
  M_cov = vcov(M_model, full = TRUE, ranpar = "sd")

  # some_a_hats[[i]] = a_hat
  # some_theta_hats[[i]] = theta_hat
  # some_M_covs[[i]] = M_cov


  ## Y
  Y_data = data.frame(
    Y,
    M,
    X,
    W1 = W[, 1],
    W2 = W[, 2],
    W3 = W[, 3],
    group = group
  )
  Y_model = glmer(Y ~ M + X + W1 + W2 + W3 + (1 + X |
                                                group),
                  data = Y_data,
                  family = binomial(link = "logit"))
  Y_model_info = attributes(VarCorr(Y_model)$group)

  ### Fitted parameters
  b_hat = fixef(Y_model)
  b_RE_sds = Y_model_info$stddev
  b_RE_cor = Y_model_info$correlation[2, 1]
  gamma_hat = c(b_RE_sds, b_RE_cor)
  if (any(is.nan(gamma_hat)))
    stop("NaNs in gamma_hat")  # Skip rest of current analysis if correlation is 0/0

  ### Estimated SE covariance matrix
  Y_cov = vcov(Y_model, full = TRUE, ranpar = "sd")

  # some_b_hats[[i]] = b_hat
  # some_gamma_hats[[i]] = gamma_hat
  # some_Y_covs[[i]] = Y_cov



  # Estimate mediation effect

  ## Fixed-effects

  a_0_hat = a_hat[1]
  a_x_hat = a_hat[2]
  A_2_hat = a_hat[3:5]

  b_0_hat = b_hat[1]
  b_m_hat = b_hat[2]
  b_x_hat = b_hat[3]
  B_3_hat = b_hat[4:6]


  ## Linear predictors
  eta_hat = as.numeric(a_0_hat + a_x_hat * x_pred + W_pred %*% A_2_hat)
  zeta_hat = as.numeric(b_0_hat + b_x_hat * x_pred + W_pred %*% B_3_hat)


  ## Random effects covariances
  s_M_0 = a_RE_sds[1]
  s_M_x = a_RE_sds[2]
  rho_M = a_RE_cor

  s_Y_0 = b_RE_sds[1]
  s_Y_x = b_RE_sds[2]
  rho_Y = b_RE_cor


  ## Sigma functions
  sigma_M1 = sigma_fun(x_pred, s_M_0, s_M_x, rho_M)
  sigma_M2 = sigma_fun(x_pred + 1, s_M_x, s_M_0, rho_M)

  sigma_Y1 = sigma_fun(x_pred, s_Y_0, s_Y_x, rho_Y)
  sigma_Y2 = sigma_fun(x_pred + 1, s_Y_x, s_Y_0, rho_Y)

  ## Mediation effect
  ### See Helpers.R for the function Phi, which computes the mediation effect on odds-ratio scale
  med_hat = Phi(eta_hat,
                zeta_hat,
                a_x_hat,
                b_m_hat,
                b_x_hat,
                sigma_M2,
                sigma_Y2,
                sigma_M1,
                sigma_Y1)
  # some_med_hats[[i]] = med_hat








  # Estimate SE

  ## Gradient of OR wrt regression coefficients
  grad_Phi_obs = grad_Phi(eta_hat,
                          zeta_hat,
                          a_x_hat,
                          b_m_hat,
                          b_x_hat,
                          sigma_M2,
                          sigma_Y2,
                          sigma_M1,
                          sigma_Y1)
  grad_xi_obs = grad_xi(a_hat, theta_hat, b_hat, gamma_hat, x_pred, W_pred)
  d_Phi_d_GLMM_pars = grad_Phi_obs %*% grad_xi_obs

  ## Get asymptotic SE using delta method

  ### Build joint covariance matrix of regression coefficients
  M_length = nrow(M_cov)
  Y_length = nrow(Y_cov)
  joint_cov = matrix(0, nrow = M_length + Y_length, ncol = M_length +
                       Y_length)
  joint_cov[1:M_length, 1:M_length] = M_cov
  joint_cov[(M_length + 1):(M_length + Y_length), (M_length + 1):(M_length +
                                                                    Y_length)] = Y_cov

  ### Convert to asymptotic covariance matrix
  asymp_reg_cov = n * joint_cov

  ### Pre- and post-multiply asymptotic covariance by gradient of Phi wrt GLMM parameters
  med_asymp_var = d_Phi_d_GLMM_pars %*% asymp_reg_cov %*% t(d_Phi_d_GLMM_pars)

  ### Get small-sample standard error
  med_asymp_SE = sqrt(med_asymp_var)
  med_SE = med_asymp_SE / sqrt(n)

  # some_med_SEs[[i]] = med_SE

  this_output = list(
    a_hat = a_hat,
    theta_hat = theta_hat,
    b_hat = b_hat,
    gamma_hat = gamma_hat,
    M_cov = M_cov,
    Y_cov = Y_cov,
    med_hat = med_hat,
    med_SE = med_SE
  )
  # }, error = function(e){grad_Phi_obs <<- e}) # For troubleshooting parallel errors
}, error = function(e) {

})                                                       ############ Un-quote





library(devtools)
load_all("D:/William/Research/MultiMedUQ/")


scale = "OR"
(this_MEs = all_MEs_models(scale, W_pred, Y_model, M_model))
(this_ME_covs = all_cov_MEs(scale, W_pred, Y_model, M_model))
