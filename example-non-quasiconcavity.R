
cov_Gc_Sc <- function (Gp, Sp, k) {
  stopifnot(length(Gp) == length(Sp))
  stopifnot(isTRUE(all.equal(mean(Gp), 0)))
  stopifnot(isTRUE(all.equal(mean(Sp), 0)))
  stopifnot(isTRUE(all.equal(var(Gp ), 1)))
  stopifnot(isTRUE(all.equal(var(Sp), 1)))
    
  Ap <- (1 - k) * Sp + k * Gp
  Bp <- k * Sp - (1 - k) * Gp
  
  Gp <- Gp[order(Ap)]
  Sp <- Sp[order(Ap)]
  
  Gp <- matrix(Gp, nrow = 2)
  Sp <- matrix(Sp, nrow = 2)
  
  Gc <- colMeans(Gp)
  Sc <- colMeans(Sp)
  
  # for covariance
  # return(cov(Gc, Sc))
  return(cor(Gc, Sc))
}
cov_Gc_Sc_v <- Vectorize(cov_Gc_Sc, vectorize.args = "k")

N <- 1e6

# normal
Gp <- rnorm(N)
Sp <- rnorm(N)
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))

# correlated
rho <- .5
zeta <- rmvnorm(N, mean = c(0, 0), sigma = matrix(c(1, rho^2, rho^2, 1), 2, 2))
Gp <- zeta[, 1]
Sp <- zeta[, 2]
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))

# L-shaped
N <- 3e5
Gp <- c(runif(N/3, 0, 1/2), runif(N/3, 0, 1/2), runif(N/3, -1/2, 0))
Sp <- c(runif(N/3, -1/2, 0), runif(N/3, 0, 1/2), runif(N/3, 0, 1/2))
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))
# other way round
Gp <- -Gp
# tests
cov_Gc_Sc(Gp, Sp, 0) - cov(Gp, Sp)
cov_Gc_Sc(Gp, Sp, 1) - cov(Gp, Sp)
cov_Gc_Sc(Gp, Sp, .5) - cov(Gp, Sp)
# here, covariance is NOT quasiconcave!!
# though it is symmetric
curve(cov_Gc_Sc_v(Gp, Sp, k), xname = "k")
curve(cov_Gc_Sc_v(Gp, Sp, k), xname = "k", from = 0.4, to = 0.6)
