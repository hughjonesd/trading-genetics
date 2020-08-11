
cor_Gc_Sc <- function (Gp, Sp, k) {
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
cor_Gc_Sc_v <- Vectorize(cor_Gc_Sc, vectorize.args = "k")
dcor_Gc_Sc <- function (Gp, Sp, k) {
  cor_Gc_Sc(Gp, Sp, k) - cor(Gp, Sp)
}
dcor_Gc_Sc_v <- Vectorize(dcor_Gc_Sc, vectorize.args = "k")


# EXAMPLES ====

# normal
N <- 1e6
Gp <- rnorm(N)
Sp <- rnorm(N)
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))

# correlated
# symmetric and concave
N <- 1e6
rho <- 0.99
zeta <- mvtnorm::rmvnorm(N, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), 2, 2))
Gp <- zeta[, 1]
Sp <- zeta[, 2]
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))

# negatively correlated
# Result: symmetric, quasiconcave, negative near k = 0 and k = 1.
N <- 1e6
rho <- - 0.9
zeta <- mvtnorm::rmvnorm(N, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), 2, 2))
Gp <- zeta[, 1]
Sp <- zeta[, 2]
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))

# L-shaped
# Result: symmetric but not quasiconcave
N <- 3e5
Gp <- c(runif(N/3, 0, 1/2), runif(N/3, 0, 1/2), runif(N/3, -1/2, 0))
Sp <- c(runif(N/3, -1/2, 0), runif(N/3, 0, 1/2), runif(N/3, 0, 1/2))
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))
Gp <- -Gp

# L-shaped and "thin" in one section, positively correlated
# Result: neither symmetric nor q-c
N <- 3e5
Gp <- c(runif(N/3, 0, 1/20), runif(N/3, 0, 1/20), runif(N/3, -1/2, 0))
Sp <- c(runif(N/3, -1/2, 0), runif(N/3, 0, 1/2), runif(N/3, 0, 1/2))
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))
Gp <- -Gp

# uniform, independent
# Result: concave and symmetric
N <- 1e6
Gp <- runif(N)
Sp <- runif(N)
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))

# non q-c distribution, independent (four "squares")
# Result: symmetric, not q-c
# hee hee 
N <- 1e6
Gp <- c(runif(N/2, 0, 1/3), runif(N/2, 2/3, 1))
Sp <- c(runif(N/2, 0, 1/3), runif(N/2, 2/3, 1))
Sp <- sample(Sp)
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))


# non q-c distribution, independent (two oblongs) but
# not invariant to rotations
# Result: not symmetric and not q-c (random?)
N <- 1e6
Gp <- c(runif(N/2, 0, 1/10), runif(N/2, 9/10, 1))
Sp <- runif(N)
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))

# q-c, symmetric in G and S, and independent (thicker in the middle) but
# not invariant to rotations
# Result: symmetric, concave
N <- 1e6
Gp <- c(runif(N/4, 0, 1/3), runif(N/2, 1/3, 2/3), runif(N/4, 2/3, 1))
Sp <- runif(N)
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))


# q-c, symmetric in G and S, independent (thick bar down the middle)
# but not invariant to rotations
# Result: not symmetric, quasiconcave
N <- 1e6
Gp <- c(
         runif(N/10, 0, 0.48), 
         runif(8 * N/10, 0.48, 0.52), 
         runif(N/10, 0.52, 1)
       )
Sp <- runif(N)
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))

# q-c, asymmetric, not independent, A triangle.
# Result: still quasi-concave
N <- 1e6
Gp <- runif(N)
Sp <- runif(N)
upper <- Gp < Sp
Gp <- Gp[upper]
Sp <- Sp[upper]
Gp <- Gp[1:4e5]
Sp <- Sp[1:4e5]
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))

# q-c, lump at top
N <- 1e6
Gp <- c(runif(N * 0.2), runif(N * 0.8, 0.9, 1))
Sp <- c(runif(N * 0.2), runif(N * 0.8, 0.9, 1))
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))


# q-c but should be thinner at rotation of 20 deg or so
N <- 1e6
Gp <- runif(N)
Sp <- runif(N)
rot <- -pi/12
m <- matrix(c(cos(rot), sin(rot), -sin(rot), cos(rot)), 2, 2)
rotated <- m %*% rbind(Gp, Sp)
Gp <- rotated[1,]
Sp <- rotated[2,]
Gp <- c(scale(Gp))
Sp <- c(scale(Sp))

# conclusions so far
# independence doesn't guarantee you Q-C or symmetry
# in general, NONE of symmetry, quasi-concavity, positive change in correlation
# are guaranteed
# maybe independence and QC guarantees you Q-C?
# (or maybe QC alone?)

# EVALUATION ====

# picture
spl <- sample(length(Gp), 5000)
plot(Gp[spl], Sp[spl])
# tests
dcor_Gc_Sc(Gp, Sp, 0)
dcor_Gc_Sc(Gp, Sp, 1)
dcor_Gc_Sc(Gp, Sp, .5)

curve(dcor_Gc_Sc_v(Gp, Sp, k), xname = "k")
