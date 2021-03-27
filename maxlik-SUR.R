library(maxLik)
library(mvtnorm)

# PSEA.y ~ a1 + PSEA1 * PSEA.x + BO1 * birth_order.x + BM1 * factor(birth_mon.x)
#            + NS1 * factor(n_sibs.x) + e1
# birth_order.y ~ a2 + b3 * PSEA.x + b4 * birth_order.x + BM2 * factor(birth_mon.x)
#            + NS2 * factor(n_sibs.x) + e2
#           
# e1, e2 ~ bivariate normal(sigma1, sigma2, rho) 
loglik_sur <- function (params, data, restricted = FALSE) {
  
  if (restricted) {
    params["BO2"] <- params["PSEA2"] * params["BO1"] / params["PSEA1"]
    if (params["PSEA1"] == 0) params["BO2"] <- 0
  }
  
  f_psea <- EA3.y ~ EA3.x + birth_order.x + moth_age_birth.x + 
              factor(birth_mon.x) + factor(n_sibs.x)
  f_bo   <- update(f_psea, birth_order.y ~ .)
  
  mm <- model.matrix(f_psea, data)
  
  yhat_psea <- params["a1"] + 
                  params["PSEA1"] * mm[, "EA3.x"] + 
                  params["BO1"] * mm[, "birth_order.x"] + 
                  params["MAB1"] * mm[, "moth_age_birth.x"] + 
                  mm[, grepl("birth_mon", colnames(mm))] %*% 
                    matrix(params[grepl("BM1", names(params))]) +
                  mm[, grepl("n_sibs", colnames(mm))] %*% 
                    matrix(params[grepl("NS1", names(params))])
  
  yhat_bo <- params["a2"] + 
               params["PSEA2"] * mm[, "EA3.x"] + 
               params["BO2"] * mm[, "birth_order.x"] + 
               params["MAB2"] * mm[, "moth_age_birth.x"] + 
               mm[, grepl("birth_mon", colnames(mm))] %*% 
                 matrix(params[grepl("BM2", names(params))]) +
               mm[, grepl("n_sibs", colnames(mm))] %*% 
                 matrix(params[grepl("NS2", names(params))])
  
  ehat_psea <- yhat_psea - data$EA3.y
  ehat_bo   <- yhat_bo   - data$birth_order.y
  
  cov_matrix <- matrix(
                  c(params["sigma1"], params["rho"], params["rho"], 
                        params["sigma2"]), 
                  2, 2
                )
  ll <- mvtnorm::dmvnorm(cbind(ehat_psea, ehat_bo), sigma = cov_matrix, log = TRUE)
  
  sum(ll)
}

# that's the log likelihood. What's the gradient?
# also a question of whether standard errors are correct when we impose
# the constraint
  
mf_pairs_sf <- mf_pairs_reg %>% 
                 filter(
                   ! is.na(birth_order.y),
                   ! is.na(EA3.y),
                   ! is.na(birth_mon.x),
                   ! is.na(EA3.x),
                   ! is.na(moth_age_birth.x)
                 )

data <- mf_pairs_sf %>% filter(! female.x)
start <- c(rep(0, 38), 1, 1, 0)
names(start) <- c("a1", "PSEA1", "BO1", "MAB1",
                  paste0("BM1.", 2:12),  paste0("NS1.", 3:6),
                  "a2", "PSEA2", "BO2", "MAB2",
                  paste0("BM2.", 2:12),  paste0("NS2.", 3:6),
                  "sigma1", "sigma2", "rho")
m <- maxLik::maxLik(loglik_sur, start = start, data = data)

omit <- which(names(start) == "BO2")
m_rest <- maxLik::maxLik(loglik_sur, start = start[-omit], data = data, restricted = TRUE)

BO2 <- coef(m_rest)["PSEA2"] * coef(m_rest)["BO1"] / coef(m_rest)["PSEA1"]

lrtest(m, m_rest)
  
