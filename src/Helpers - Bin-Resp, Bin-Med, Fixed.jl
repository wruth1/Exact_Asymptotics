
# Compute mediation effect

## Probability of Y=1 given inputs
function PY1(eta, zeta, beta)
    Q1 = 1 + exp(-zeta)
    Q2 = 1 + exp(-eta)
    Q3 = 1 + exp(-zeta - beta)

    num = Q1 + (Q2 - 1) * Q3
    den = Q1 * Q2 * Q3

    return num / den
end

## Probability of Y=0 given inputs
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

