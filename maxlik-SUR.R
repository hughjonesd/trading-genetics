
library(maxLik)
library(mvtnorm)
library(dplyr)
library(rumpel) # just to play with it!


loglik_sur <- function (params, mm, y_psea, y_bo, restricted = NULL) {
  
  if (! is.null(restricted)) {
    to_restrict <- paste("bo", restricted, sep = "_")
    from_restrict <- paste("psea", restricted, sep = "_")
    # params["bo_birth_order.x"] will be in 'fixed',
    # so now we calculate it
    params[to_restrict] <- params[from_restrict] * 
                             params["bo_EA3.x"] / 
                             params["psea_EA3.x"]
    if (params["psea_EA3.x"] == 0) params[to_restrict] <- 0
  }
  
  psea_params <- params %>% named_starting("psea")
  yhat_psea <- mm %*% matrix(psea_params)
  bo_params <- params %>% named_starting("bo")
  yhat_bo   <- mm %*% matrix(bo_params)

  ehat_psea <- yhat_psea - y_psea
  ehat_bo   <- yhat_bo   - y_bo
  
  cov_matrix <- matrix(
                  c(params["sigma1"], params["rho"], params["rho"], 
                    params["sigma2"]), 
                  2, 2
                )
  
  ll <- mvtnorm::dmvnorm(
          cbind(ehat_psea, ehat_bo), 
          sigma = cov_matrix, 
          log   = TRUE
        )
  
  ll
}


#' Estimate full and restricted SUR regressions on spouse PSEA and birth order,
#' using maximum likelihood
#'
#' @param fml One-sided formula of independent variables
#' @param data Data frame
#' @param restricted Character: regex. Matching independent variables will be
#'   restricted.
#' 
#' In the restricted model, the coefficient of restricted variables 
#' on birth order will not be estimated, rather calculated as "effect of
#' X on PSEA * effect of PSEA on birth order / effect of PSEA on PSEA".
#' 
#' Estimation is done by [maxLik::maxLik()].
#' 
#' @return 
#' A list of two models, "full" and "restricted". The "restricted"
#' attribute of the "restricted" model contains calculated values for restricted
#' coefficients.
#'
estimate_surs <- function (fml, data, restricted) {
  
  mm <- model.matrix(fml, data)
  
  # Ensure EA3.y and birth_order.y stay in `data`:
  fml_ext <- update(fml, ~ . + EA3.y + birth_order.y) 
  # Get rid of NAs:
  data <- model.frame(fml_ext, data)
  
  restricted <- grep(restricted, colnames(mm), value = TRUE)
  
  start <- rep(0, ncol(mm) * 2)
  names(start) <- paste(
                    rep(c("psea", "bo"), each = ncol(mm)),
                    rep(colnames(mm), 2),
                    sep = "_"
                  )
  start <- c(start, sigma1 = 1, sigma2 = 1, rho = 0)
  
  # constraints: sigma1 > 0, sigma2 > 0
  constraintsA <- matrix(0, 2, length(start))
  constraintsA[1, length(start) - 2] <- 1
  constraintsA[2, length(start) - 1] <- 1
  constraintsB <- matrix(0, 2, 1)
  
  parscale <- rep(1, length(start))
  parscale[grepl("birth_order.x", names(start))] <- 10
  parscale[grepl("EA3.x", names(start))] <- 10
  
  mod_full <- maxLik::maxLik(
                loglik_sur, 
                start       = start, 
                mm          = mm, 
                y_psea      = data$EA3.y,
                y_bo        = data$birth_order.y,
                constraints = list(
                                ineqA = constraintsA, 
                                ineqB = constraintsB
                              ),
                method      = "BFGS",
                parscale    = parscale
              )
  
  bo_restricted <- paste("bo", restricted, sep = "_")
  fixed_index <- which(names(start) %in% bo_restricted)
  mod_restricted <- maxLik::maxLik(
                      loglik_sur, 
                      start       = start, 
                      mm          = mm, 
                      y_psea      = data$EA3.y,
                      y_bo        = data$birth_order.y,
                      restricted  = restricted,
                      constraints = list(
                                      ineqA = constraintsA, 
                                      ineqB = constraintsB
                                    ),
                      method      = "BFGS",
                      parscale    = parscale,
                      fixed       = fixed_index
                    )
  
  coefs <- coef(mod_restricted)
  psea_restricted <- paste("psea", restricted, sep = "_")
  r_coefs <- coefs[psea_restricted] * coefs["bo_EA3.x"] / coefs["psea_EA3.x"]
  attr(mod_restricted, "restricted") <- setNames(r_coefs, bo_restricted)
  
  list_mod <- list(full = mod_full, restricted = mod_restricted)
  
  list_mod
}

mf_pairs_sf <- mf_pairs_reg %>% 
                 filter(
                   ! is.na(birth_order.y),
                   ! is.na(EA3.y),
                   ! is.na(birth_mon.x),
                   ! is.na(EA3.x),
                   ! is.na(moth_age_birth.x)
                 )


fml <- ~ EA3.x + factor(birth_order.x) + moth_age_birth.x + factor(birth_mon.x) + 
           factor(n_sibs.x) + university.x + fluid_iq.x + height.x

list_mod_male <- estimate_surs(fml, mf_pairs_sf %>% filter(! female.x), 
                                 restricted = "birth_order.x")

summary(list_mod_male$full)
summary(list_mod_male$restricted)
attr(list_mod_male$restricted, "restricted")
lmtest::lrtest(list_mod_male$full, list_mod_male$restricted)


list_mod_female <- estimate_surs(fml, mf_pairs_sf %>% filter(female.x), 
                                 restricted = "birth_order.x")

summary(list_mod_female$full)
summary(list_mod_female$restricted)
attr(list_mod_female$restricted, "restricted")
lmtest::lrtest(list_mod_female$full, list_mod_female$restricted)
