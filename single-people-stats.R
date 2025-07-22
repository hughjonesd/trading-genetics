
# Questions asked by Oana

# What about single people?
# We'll look at their birth order and maybe PGS

drake::loadd(famhist)
library(dplyr)
library(fixest)

# Analogous to our original sample
famhist <- famhist |> filter(birth_order <= 6)


## Stats which don't control for number of siblings

tb <- table(famhist$with_partner, famhist$birth_order)
tb <- proportions(tb, 1)
round(tb, 3)
barplot(tb, beside = TRUE, legend.text = c("No partner", "Has partner"))

t.test(birth_order ~ with_partner, data = famhist)
wilcox.test(birth_order ~ with_partner, data = famhist)

t.test(EA3 ~ with_partner, data = famhist)


## Controlling for number of siblings, as in Table 2 of the paper

fml_with_pnr <- list()
fml_with_pnr[[1]] <- with_partner ~ birth_order | factor(n_sibs)
fml_with_pnr[[2]] <- with_partner ~ birth_order + EA3 + factor(birth_mon) |
  factor(n_sibs) + factor(age_at_recruitment)
fml_with_pnr[[3]] <- with_partner ~ birth_order + EA3 + factor(birth_mon) +
  par_age_birth | factor(n_sibs) + factor(age_at_recruitment)

fml_with_pnr <- lapply(fml_with_pnr, Formula::as.Formula)

mod_with_pnr <- lapply(fml_with_pnr, fixest::feols, 
                           data = famhist,
                           notes = FALSE)
huxtable::huxreg(mod_with_pnr, coefs = c("birth_order", "EA3", "par_age_birth"), 
                 note = "Other controls: dummies for birth month, number of siblings, and own age")


