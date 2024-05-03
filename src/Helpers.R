# General ####

logit = function(x) log(x/(1-x))

expit = function(x) exp(x)/(1+exp(x))



# Cont Resp, Bin Med, Fixed-Effs ####
## Some notation:
  ## eta: a_0 + a_1 x + A_2 W = a_0 + a_x * x + A_2' w. Linear predictor for M
  ## zeta: b_0 + b_1 * m + b_2 * x + B_3' w = b_0 + b_m * m + b_x * x + B_3' w. Linear predictor for Y

  ## delta: E(M | X=x+1) - E(M | X=x)
  ## gamma: E(Y | X=x+1) - E(Y | X=x). Total effect of X on Y


get_delta = function(eta, a_x){
  Q1 = exp(-eta)
  Q2 = (1 - exp(-a_x))
  Q3 = (1 + exp(-eta - a_x))
  Q4 = (1 + exp(-eta))

  return(Q1 * Q2 / (Q3 * Q4))
}


get_gamma = function(delta, b_m, b_x){
  return(delta * b_m + b_x)
}


## Gradient of total effect wrt the a's and b's

### First, some intermediate quantities
d_delta_d_eta = function(eta, a_x){
  delta = get_delta(eta, a_x)

  Q1 = 1 - exp(-2*eta - a_x)
  Q2 = 1 + exp(-eta)
  Q3 = 1 + exp(-eta - a_x)

  return( (-1) * delta * Q1 / (Q2 * Q3))
}

d_delta_d_a_x = function(eta, a_x){
  Q1 = exp(-2*eta - a_x)
  Q2 = exp(-eta - a_x)
  Q3 = 1 + exp(-eta)
  Q4 = 1 + exp(-eta - a_x)

  return((Q1 + Q2)/(Q3 * Q4^2))
}


### Now, the gradient
d_gamma_d_a_0 <- function(dd_de, b_m) {
  return(dd_de * b_m)
}

d_gamma_d_a_1 <- function(dd_de, dd_da_x, x, b_m) {
  return(b_m * x * dd_de + b_m * dd_da_x)
}

d_gamma_d_A_2 <- function(dd_de, W, b_m) {
  return(c(dd_de * b_m) * W)
}




d_gamma_d_b_0 <- function() {
  return(0)
}

d_gamma_d_b_1 <- function(eta, a_x) {
  return(get_delta(eta, a_x))
}

d_gamma_d_b_2 <- function() {
  return(1)
}

d_gamma_d_B_3 <- function(W) {
  return(rep(0, times = length(W)))
}

d_gamma_d_theta <- function(eta, x, W, a, b) {
  a_0 <- a[1]
  a_x <- a[2]
  A_2 <- a[3:length(a)]

  b_0 <- b[1]
  b_m <- b[2]
  b_x <- b[3]
  B_3 <- b[4:length(b)]

  dd_de <- d_delta_d_eta(eta, a_x)
  dd_da_x <- d_delta_d_a_x(eta, a_x)

  dg_da_0 <- d_gamma_d_a_0(dd_de, b_m)
  dg_da_1 <- d_gamma_d_a_1(dd_de, dd_da_x, x, b_m)
  dg_dA_2 <- d_gamma_d_A_2(dd_de, W, b_m)

  dg_db_0 <- d_gamma_d_b_0()
  dg_db_1 <- d_gamma_d_b_1(eta, a_x)
  dg_db_2 <- d_gamma_d_b_2()
  dg_dB_3 <- d_gamma_d_B_3(W)

  return(c(dg_da_0, dg_da_1, dg_dA_2, dg_db_0, dg_db_1, dg_db_2, dg_dB_3))
}
