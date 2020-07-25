
# find the "assortative mating score"

library(drake)
library(magrittr)

drake::loadd(score_names)
drake::loadd(mf_pairs)

mf_pairs$score.m <- mf_pairs$EA3_excl_23andMe_UK_resid.m
mf_pairs$score.f <- mf_pairs$EA3_excl_23andMe_UK_resid.f


score_names <- paste0(score_names, "_resid")
fml_m <- as.formula(paste("score.m ~ ", paste0(score_names, ".f", collapse = " +")))
fml_f <- as.formula(paste("score.f ~ ", paste0(score_names, ".m", collapse = " +")))

coef_m <- rep(-100, 34)
coef_f <- rep(-100, 34)
epsilon <- 1e-5
r <- 0
while (r <= 100) {
  mod_m <- lm(fml_m, mf_pairs, na.action = na.exclude)
  mod_f <- lm(fml_f, mf_pairs, na.action = na.exclude)
  delta <- max(c(abs(coef(mod_m) - coef_m), abs(coef(mod_f) - coef_f)))
  if (delta < epsilon) break
  mf_pairs$score.m <- predict(mod_f)
  mf_pairs$score.f <- predict(mod_m)
  mf_pairs[, c("score.m", "score.f")] %<>% scale()
  coef_m <- coef(mod_m)
  coef_f <- coef(mod_f)
  cat(Sys.time(), " : ", r, " ", delta, "\n")
  r <- r + 1
}
