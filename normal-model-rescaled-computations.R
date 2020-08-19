
z <- function (k) sqrt(k^2 + (1-k)^2)
C <- function (k) k/z(k)
d <- function (k)(1-k)/z(k)

rho2_AB <- function (k, rho) {
  C <- C(k)
  d <- d(k)
  (C - d)^2 * rho^2 / z(k)^2 / (1 + 2*C*d*rho) / (1 - 2*C*d*rho)
}

var_Gc <- function (k, rho) {
  d <- d(k)
  C <- C(k)
  
  d^2 *  (1 + 2 * C * d * rho) + 
  C^2 * (1/2 + rho2_AB(k, rho)/2) * (1 - 2 * C * d * rho) +
  2 * C * d * (C - d) * rho / z(k)
}


var_Sc <- function (k, rho) {
  d <- d(k)
  C <- C(k)
  
  C^2 *  (1 + 2 * C * d * rho) + 
  d^2 * (1/2 + rho2_AB(k, rho)/2) * (1 - 2 * C * d * rho) -
  2 * C * d * (C - d) * rho / z(k)
}

denom <- function (k, rho) sqrt(var_Sc(k, rho) * var_Gc(k, rho))

cov_Gc_Sc <- function (k, rho) {
  d <- d(k)
  C <- C(k)
  
  nasty_bit <- (C - d)^2 * rho^2 / z(k)^2 / (1 + 2 * C * d * rho)
  rho + (C * d / 2)*(1 - 2 * C * d * rho - nasty_bit)
  
}

corr_Gc_Sc <- function (k, rho) {
  cov_Gc_Sc(k, rho)/denom(k, rho)
}

plot_both <- function (rho, ylim = c(-1, 1)) {
  curve(cov_Gc_Sc(k, rho), xname = 'k', col = 'darkgreen', 0, 0.5, ylim = ylim)
  curve(denom(k, rho), xname = 'k', add = TRUE, col = 'red')
}
