

sigma_I <- function (s, S, sigma, a) {
  sqrt(
    a^2 * s^2 + (1-a)^2 * S^2 + 2*a*(1-a)*sigma
  )
}

A <- function (s, S, sigma, a) (a*s^2 + (1-a)*sigma)/sigma_I(s, S, sigma, a)

C <- function (s, S, sigma, a) (a*sigma + (1-a)*S^2)/sigma_I(s, S, sigma, a)


cov_x1_x2 <- function (s, S, sigma, a, gamma, tau, theta) {
  A <- A(s, S, sigma, a)
  C <- C(s, S, sigma, a)
  
  gamma * tau^2/2 * (s^2 + A^2) + 
  theta * tau/2 *(sigma + A*C) + 
  gamma
}

corr_x1_x2 <- function (s, S, sigma, a, gamma, tau, theta) {
  A <- A(s, S, sigma, a)
  C <- C(s, S, sigma, a)
  
  cov_x1_x2 <- cov_x1_x2(s, S, sigma, a, gamma, tau, theta)
  
  var_x1 <- tau^2/2 * (s^2+A^2) + 1
  
  var_x2 <- gamma^2 * tau^2/2 * (s^2 + A^2) +
            gamma * tau * theta * (sigma + A*C) +
            theta^2/2 * (S^2 + C^2) +
            1 + gamma^2
  
  cov_x1_x2/sqrt(var_x1 * var_x2)
}

corr_x1_x2_a_b <- function (s, S, sigma, a, b, tau, theta) {
  A <- A(s, S, sigma, a)
  B <- A(s, S, sigma, b)
  C <- C(s, S, sigma, a)
  D <- C(s, S, sigma, b)
  
  cov <- (tau * theta / 4) * (2*sigma + A*D + B*C)
  var_s <- (tau^2 / 2) * (s^2 + A*B + 1)
  var_S <- (theta^2 / 2) * (S^2 + C*D + 1)
  
  cov/sqrt(var_s*var_S)
}


long_run_corr <- function (a, gamma, tau, theta) {
  s2 <- 1
  S2 <- 1
  sigma <- 0
  
  diff <- Inf
  epsilon <- 1e-5
  while (diff > epsilon) {
    old_s2 <- s2
    old_S2 <- S2
    old_sigma <- sigma
    old_corr <- sigma/sqrt(s2*S2)

    A <- A(sqrt(s2), sqrt(S2), sigma, a)
    C <- C(sqrt(s2), sqrt(S2), sigma, a)
    s2 <- tau^2/2 * (s2+A^2) + 1
    S2 <- gamma^2 * tau^2/2 * (s2 + A^2) +
            gamma * tau * theta * (sigma + A*C) +
            theta^2/2 * (S2 + C^2) +
            1 + gamma^2
    sigma <- cov_x1_x2(sqrt(s2), sqrt(S2), sigma, a, gamma, tau, theta)
    
    corr <- sigma/sqrt(s2*S2)
    diff <- abs(corr - old_corr)
  }
  
  return(corr)
}

# example of where theta decreases correlation for positive gamma
curve(corr_x1_x2(1, 1, sigma = 0, a = .1, gamma = 0.25, tau = 0.95, theta = x), ylim = c(0, 0.33))
curve(corr_x1_x2(1, 1, sigma = 0, a = .1, gamma = 0, tau = 0.95, theta = x), add=T,col='red')

# new stuff
# 
#  
lambda <- function (mu, a, theta, tau) {
  X <- (1/2 * theta * tau - 1)
  numer <- mu * X * ((1-a)^2 + 2 * a * (1-a) * mu) + 
    1/2 * theta * tau * mu * (1-a) * (1 - a + a*mu)
  
  denom <- (a^2 * mu * -X) - (a/2 * theta * tau * (1 - a + a*mu))
  
  numer/denom
}

F <- function (mu, a, theta, tau) {
  lambda <- lambda(mu = mu, a = a, theta = theta, tau = tau)
  
  sigma_I_ish <- a^2 * lambda + (1-a)^2 + 2*a*(1-a)*mu
  
  A_ish <- (a * lambda + (1-a) * mu)^2 / sigma_I_ish
  
  C_ish <- (a * mu + (1-a))^2 / sigma_I_ish
  
  (lambda - 1/2 * tau^2 * (A_ish + lambda)) -
    (1 - 1/2 * theta^2 * (C_ish + 1))
}

phi <- function (a, theta, tau) {
  F_to_solve <- \(mu) F(mu, a, theta, tau)
  roots <- fstats::uniroot(F_to_solve, lower = 0, upper = 10)
}

plot_F <- function (a, theta, tau, from = 0, to = 10, ...) {
  F_to_solve <- \(mu) F(mu, a, theta, tau)
  plot(F_to_solve, from = from, to = to, ...)
}
