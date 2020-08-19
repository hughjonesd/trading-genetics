

q <- function (k, rho) {
  2 * rho*k*(1-k)  
}

z <- function (k) k^2 + (1-k)^2

rho2_AB <- function (k, rho) {
  z <- z(k)
  (2*k - 1)^2 * rho^2 / (z + q(k, rho)) / (z - q(k, rho))
}

var_B <- function (k, rho) {
  z(k) - q(k, rho)  
}

dcov_GS <- function (k, rho) {
  k * (1-k) / z(k)^2 * (1 - rho2_AB(k, rho)) * var_B(k, rho)
}

corr_GS <- function (k, rho) {
  cov_GS <- rho + dcov_GS(k, rho)
  
  rho2_AB <- rho2_AB(k, rho)
  var_A <- z(k) + q(k, rho)
  var_Bp <- z(k) - q(k, rho)
  var_Bc <- (1/2 + rho2_AB/2) * var_Bp
  cov_ABc <- sqrt(rho2_AB) * sqrt(1/2 + rho2_AB/2) * sqrt(var_A*var_Bc)
    
  c <- k/z(k)
  d <- (1-k)/z(k)
  var_G <- c^2 * var_A + d^2 * var_Bc + 2 * c * d * cov_ABc
  var_S <- d^2 * var_A + c^2 * var_Bc - 2 * c * d * cov_ABc
  
  cov_GS / sqrt(var_G * var_S)
}

rho <- -0.9
curve(rho + dcov_GS(k, rho), xname = "k")
abline(h = rho, col = "darkgreen", lty = 2)

