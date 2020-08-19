
library(drake)
library(dplyr)
library(fixest)
library(huxtable)

loadd(famhist)
loadd(mf_pairs)
loadd(mf_pairs_twice)

famhist$birth_year_m <- famhist$birth_year * 12 + famhist$birth_mon

famhist$dbym <- famhist$birth_year_m - (1957 * 12 + 9)
famhist$rosla <- famhist$dbym >= 0


mod_rosla <- lm(university ~ rosla + rosla*dbym + rosla*I(dbym^2) + dbym + 
                  I(dbym^2), famhist)
summary(mod_rosla)




fml_bo_psea <- list()

fml_bo_psea[[1]] <- EA3.y ~ n_older_sibs.x + EA3.x | factor(n_sibs.x)
fml_bo_psea[[2]] <- EA3.y ~ age_fulltime_edu.x + university.x + n_older_sibs.x + EA3.x | factor(n_sibs.x)
# fml_bo_psea[[3]] <- EA3.y ~ income_cat.x + n_older_sibs.x + EA3.x | factor(n_sibs.x)
# fml_bo_psea[[4]] <- EA3.y ~ income_cat.x + age_fulltime_edu.x + n_older_sibs.x + EA3.x |
#       factor(n_sibs.x)

mf_pairs_reg <- mf_pairs_twice %>% 
                  filter(
                    n_sibs.x >= 2, 
                    n_sibs.x <= 8, 
                    ! is.na(n_older_sibs.x),
                    ! is.na(university.x),
                    ! is.na(age_fulltime_edu.x),
                    ! is.na(income_cat.x)
                  )

mod_bo_psea <- lapply(fml_bo_psea, fixest::feols, 
                 data = mf_pairs_reg,
                 notes = FALSE
               )


huxreg(mod_bo_psea, 
        coefs = c(
          "Birth order" = "n_older_sibs.x", 
          "University"  = "university.xTRUE",
          "Age FTE"     = "age_fulltime_edu.x",
          "Own EA3"     = "EA3.x"
        ),
         note = "{stars}. Standard errors clustered by spouse pair.",
         tidy_args = list(conf.int = FALSE, cluster = list(mf_pairs_reg$couple_id))
       ) 
