
# Compute mediation effect

## Increment in conditional expectation of M for a unit increase in X
function get_delta(eta, a_x)
    Q1 = exp(-eta)
    Q2 = (1 - exp(-a_x))
    Q3 = (1 + exp(-eta - a_x))
    Q4 = (1 + exp(-eta))

    return Q1 * Q2 / (Q3 * Q4)
end

## Increment in conditional expectation of Y for a unit increase in X
## I.e. Mean over M of increment in conditional expectation for Y given M and X
function get_gamma(delta, b_m, b_x)
    return b_m * delta + b_x
end


# Compute SE of estimated mediation effect
# We use the delta method, so many derivatives are required


## Helpful intermedeate quantities
function d_delta_d_eta(eta, a_x)
    delta = get_delta(eta, a_x)

    Q1 = 1 - exp(-2*eta - a_x)
    Q2 = 1 + exp(-eta)
    Q3 = 1 + exp(-eta - a_x)

    return - delta * Q1 / (Q2 * Q3)
end

function d_delta_d_a_x(eta, a_x)
    Q1 = exp(-2*eta - a_x)
    Q2 = 1 + exp(-eta)
    Q3 = 1 + exp(-eta - a_x)

    return Q1 / (Q2 * Q3^2)
end


## Derivative of gamma wrt each parameter
function d_gamma_d_a_0(dd_de, b_x)
    return dd_de * b_x
end

function d_gamma_d_a_1(dd_de, dd_da_x, x, b_x)
    return b_x * x * dd_de + b_x * dd_da_x
end

function d_gamma_d_A_2(dd_de, W, b_x)
    return dd_de * b_x .* W
end


function d_gamma_d_b_0()
    return 0
end

function d_gamma_d_b_1(eta)
    return eta
end

function d_gamma_d_b_2()
    return 1
end

function d_gamma_d_B_3(W)
    return 0 .* W
end



## Gradient of gamma wrt all parameters
function d_gamma_d_theta(eta, x, W, a, b)
    a_0, a_x = a[1:2]
    A_2 = a[3:end]

    b_0, b_m, b_x = b[1:3]
    B_3 = b[4:end]

    dd_de = d_delta_d_eta(eta, a_x)
    dd_da_x = d_delta_d_a_x(eta, a_x)


    dg_da_0 = d_gamma_d_a_0(dd_de, b_x)
    dg_da_1 = d_gamma_d_a_1(dd_de, dd_da_x, x, b_x)
    dg_dA_2 = d_gamma_d_A_2(dd_de, W, b_x)

    dg_db_0 = d_gamma_d_b_0()
    dg_db_1 = d_gamma_d_b_1(eta)
    dg_db_2 = d_gamma_d_b_2()
    dg_dB_3 = d_gamma_d_B_3(W)

    return [dg_da_0; dg_da_1; dg_dA_2; dg_db_0; dg_db_1; dg_db_2; dg_dB_3]
end












function PY1(eta, zeta, beta)
    Q1 = 1 + exp(-zeta)
    Q2 = 1 + exp(-eta)
    Q3 = 1 + exp(-zeta - beta)

    num = Q1 + (Q2 - 1) * Q3
    den = Q1 * Q2 * Q3

    return num / den
end

function PY0(eta, zeta, beta)
    Q1 = exp(-zeta)
    Q2 = exp(-eta)
    Q3 = exp(-beta)

    num = Q1 * (Q2 + Q3 + Q1 * Q3 + Q1 * Q2 * Q3)
    den = (1 + Q1) * (1 + Q2) * (1 + Q1*Q3)

    return num / den
end


function get_odds_ref(eta, zeta, beta)
    num = PY1(eta, zeta, beta)
    den = PY0(eta, zeta, beta)

    return num / den

end

function get_odds(eta, zeta, beta)
    Q1 = exp(-zeta)
    Q1_inv = 1/Q1
    Q2 = exp(-eta)
    Q3 = exp(-beta)

    num = 1 + Q1_inv + Q2 * (Q3 + Q1_inv)
    den = Q2 + Q3 + Q1 * Q3 + Q1 * Q2 * Q3

    return num / den
end

function get_odds_ratio(eta, a_x, zeta, b_x, b_m)
    odds_1 = get_odds(eta + a_x, zeta + b_x, b_m)
    odds_2 = get_odds(eta, zeta, b_m)
    OR_hat = odds_1 / odds_2
end



# Gradients of OR, obtained from Maple
## I verified all of these against a simple finite difference calculation

function d_OR_d_eta(eta, a_x, zeta, b_x, b_m)
    -exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m)) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) - (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) ^ 2 / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-eta - a_x) - exp(-zeta - b_x - eta - a_x - b_m)) + (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) ^ 2 * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * exp(-eta) * (1 + exp(-zeta - b_m)) + (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (-exp(-eta) - exp(-zeta - eta - b_m))
end



function d_OR_d_a_x(eta, a_x, zeta, b_x, b_m)
    -exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m)) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) - (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) ^ 2 / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-eta - a_x) - exp(-zeta - b_x - eta - a_x - b_m))
end


function d_OR_d_zeta(eta, a_x, zeta, b_x, b_m)
    (-exp(-zeta - b_x) - exp(-eta - a_x) * exp(-zeta - b_x - b_m)) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) - (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) ^ 2 / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-zeta - b_x - b_m) - exp(-zeta - b_x - eta - a_x - b_m)) - (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) ^ 2 * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-zeta) - exp(-eta) * exp(-zeta - b_m)) + (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (-exp(-zeta - b_m) - exp(-zeta - eta - b_m))
end


function d_OR_d_b_x(eta, a_x, zeta, b_x, b_m)
    (-exp(-zeta - b_x) - exp(-eta - a_x) * exp(-zeta - b_x - b_m)) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) + (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) - (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) ^ 2 / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-zeta - b_x - b_m) - exp(-zeta - b_x - eta - a_x - b_m))
end


function d_OR_d_b_m(eta, a_x, zeta, b_x, b_m)
    -exp(-eta - a_x) * exp(-zeta - b_x - b_m) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) - (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) ^ 2 / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-b_m) - exp(-zeta - b_x - b_m) - exp(-zeta - b_x - eta - a_x - b_m)) + (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) ^ 2 * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * exp(-eta) * exp(-zeta - b_m) + (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (-exp(-b_m) - exp(-zeta - b_m) - exp(-zeta - eta - b_m))
end





# Build gradient wrt reg pars
function d_OR_d_a0(eta, a_x, zeta, b_x, b_m)
    return d_OR_d_eta(eta, a_x, zeta, b_x, b_m)
end

function d_OR_d_a1(eta, a_x, zeta, b_x, b_m, x_ref)
    return x_ref * d_OR_d_eta(eta, a_x, zeta, b_x, b_m) + d_OR_d_a_x(eta, a_x, zeta, b_x, b_m)
end

function d_OR_d_A2(eta, a_x, zeta, b_x, b_m, W_ref)
    return W_ref .* d_OR_d_eta(eta, a_x, zeta, b_x, b_m)
end


function d_OR_d_b0(eta, a_x, zeta, b_x, b_m)
    return d_OR_d_zeta(eta, a_x, zeta, b_x, b_m)
end

function d_OR_d_b1(eta, a_x, zeta, b_x, b_m)
    return d_OR_d_b_m(eta, a_x, zeta, b_x, b_m)
end

function d_OR_d_b2(eta, a_x, zeta, b_x, b_m, x_ref)
    return x_ref * d_OR_d_zeta(eta, a_x, zeta, b_x, b_m) + d_OR_d_b_x(eta, a_x, zeta, b_x, b_m)
end

function d_OR_d_B3(eta, a_x, zeta, b_x, b_m, W_ref)
    return W_ref .* d_OR_d_zeta(eta, a_x, zeta, b_x, b_m)
end


function d_OR_d_theta(eta, a_x, zeta, b_x, b_m, x_ref, W_ref)
    return [d_OR_d_a0(eta, a_x, zeta, b_x, b_m); d_OR_d_a1(eta, a_x, zeta, b_x, b_m, x_ref); d_OR_d_A2(eta, a_x, zeta, b_x, b_m, W_ref); d_OR_d_b0(eta, a_x, zeta, b_x, b_m); d_OR_d_b1(eta, a_x, zeta, b_x, b_m); d_OR_d_b2(eta, a_x, zeta, b_x, b_m, x_ref); d_OR_d_B3(eta, a_x, zeta, b_x, b_m, W_ref)]
end

