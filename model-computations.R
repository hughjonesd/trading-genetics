

type_names <- c("b", "g", "s", "n")

type_to_GS <- matrix(c(1, 1, 1, 0, 0, 1, 0, 0), 2, 4,
                     dimnames = list(c("G", "S"), type = type_names))


GS_to_type <- matrix(c("n", "g", "s", "b"), 2, 2, 
                     dimnames = list(G = c("0", "1"), S = c("0", "1")))

G_types <- type_to_GS["G", type_names]
G_parents <- outer(G_types, G_types, `+`)
S_types <- type_to_GS["S", type_names]
S_parents <- outer(S_types, S_types, `+`)

# types of children from types of parents
calc_children_types <- function (types, alpha) {
  stopifnot(sum(types) == 1, all(types >= 0))
  stopifnot(names(types) == type_names)
  

  # copy names and dimensions:
  G_children <- S_children <- G_parents
  # proportion of children with G is average of parents G: 1, 1/2 or 0
  # proportion of children with S is 1 if both parents have S = 1, alpha if only one does, 0 otherwise
  G_children[] <- c(0, 0.5, 1)[G_parents + 1]   # + 1 for 1-indexing
  S_children[] <- c(0, alpha, 1)[S_parents + 1]
  
  # within each set of parents, G and S of children are independent,
  # so we can do this:
  type_children <- list(
                      b = G_children * S_children, 
                      g = G_children * (1 - S_children),
                      s = (1 - G_children) * S_children,
                      n = (1 - G_children) * (1 - S_children)
                    )
  
  # now we can multiply by the proportions of each parent
  parents_props <- types %o% types
  type_children <- lapply(type_children, `*`, parents_props)
  type_children <- vapply(type_children, sum, numeric(1))
  
  type_children
}


calc_children_from_params <- function (gamma, sigma, theta, alpha, k) {
  stopifnot(k >= 0, k <= 1)
  stopifnot(alpha >= 0, alpha <= 1)
  stopifnot(gamma > 0, gamma < 1, sigma > 0, sigma < 1)
  stopifnot(theta > gamma * sigma, theta <= sigma, theta <= gamma)
  
  
  # just for checking:
  types <- c(
              b = gamma * sigma, 
              g = gamma * (1 - sigma), 
              s = (1 - gamma) * sigma, 
              n = (1 - gamma) * (1 - sigma)
            )
  
  phi <- theta - sigma * gamma
  types_H <- c(
                b = gamma * sigma / theta, 
                g = (1 - k) * phi/theta, 
                s = k * phi/theta, 
                n = 0
              )
  types_L <- c(
                b = 0, 
                g = (gamma * (1 - sigma) - (1 - k) * phi)/(1 - theta),
                s = ((1 - gamma) * sigma - k * phi)/(1 - theta),
                n = (1 - gamma) * (1 - sigma)/(1 - theta)
              )
            
  types_H_chn <- calc_children_types(types_H, alpha = alpha)
  types_L_chn <- calc_children_types(types_L, alpha = alpha)
  types_chn <- theta * types_H_chn + (1 - theta) * types_L_chn
  
  # only for alpha = .5:
  # stopifnot(isTRUE(all.equal(types_chn["b"], types["b"] + types["g"] * types["s"]/2)))
            
  return(types_chn)
}

calc_GS_cor <- function(types) {
  with(as.list(types), 
    b - (b + g)*(b + s)
  )
}

calc_cor <- Vectorize(function (gamma, sigma, theta, alpha, k) {
  calc_GS_cor(calc_children_from_params(gamma, sigma, theta, alpha, k))
})


calc_children_algebra <- function (gamma, sigma, theta, alpha, k) {
  
  b <- gamma * sigma
  g <- gamma * (1 - sigma)
  s <- (1 - gamma) * sigma
  n <- (1 - gamma) * (1 - sigma)
  
  b_H <- b/theta
  g_H <- (1 - k)*(theta - gamma*sigma)/theta
  types_chn <- numeric(4)
  names(types_chn) <- type_names

  types_chn["b"] <- b - theta * (1 - 2 * alpha) * b_H * g_H + alpha * g * s
  types_chn["g"] <- g + theta * (1 - 2 * alpha) * b_H * g_H - alpha * g * s
  types_chn["s"] <- s - (1 - alpha) * g * s
  types_chn["n"] <- n + (1 - alpha) * g * s

  types_chn
}

gamma <- .2
sigma <- 0.6
theta <- .2
k <- .35
curve(calc_cor(gamma, sigma, theta, alpha, k), xname = "alpha", ylim = c(0, .2))
