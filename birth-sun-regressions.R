
# regressions with birth sun years. 
# Seemingly a null result: once you control tightly enough for geography,
# the effect of sunshine during pregnancy goes away.
# Furthermore, there is a correlation with EA3 which goes away at the same time.

library(drake)
library(santoku)
library(fixest)
library(dplyr)
library(huxtable)

drake::loadd(famhist)
famhist$birth_sun_s <- c(scale(famhist$birth_sun))

famhist$birth_lat_10k <- santoku::chop_width(famhist$birth_lat.0.0, 10000)
famhist$birth_lon_10k <- santoku::chop_width(famhist$birth_lon.0.0, 10000)
famhist$birth_lat_50k <- santoku::chop_width(famhist$birth_lat.0.0, 50000)
famhist$birth_lon_50k <- santoku::chop_width(famhist$birth_lon.0.0, 50000)

 

mod_uni_sun <- mod_income_sun <- mod_ea3_sun <- list()

# subset arg to lm doesn't work for poly(), which can't deal with NA values:
mod_uni_sun$poly <- famhist %>% 
                      filter(! is.na(birth_lat.0.0), ! is.na(birth_lon.0.0)) %>% 
                      {fixest::feglm(university ~ 
                                    birth_sun_s + 
                                    poly(birth_lat.0.0, birth_lon.0.0, degree = 5) |
                                    factor(birth_year) + 
                                    factor(birth_mon), 
                                  data   = .,
                                  family = "binomial"
                                )}

mod_income_sun$poly <- famhist %>% 
                        filter(! is.na(birth_lat.0.0), ! is.na(birth_lon.0.0)) %>% 
                        {lm(income_cat ~ 
                                      birth_sun_s + 
                                      poly(birth_lat.0.0, birth_lon.0.0, degree = 5) +
                                      factor(birth_year) + 
                                      factor(birth_mon), 
                                    data = .
                        )}

mod_ea3_sun$poly <- famhist %>% 
                      filter(! is.na(birth_lat.0.0), ! is.na(birth_lon.0.0)) %>% 
                      {lm(EA3 ~ 
                                    birth_sun_s + 
                                    poly(birth_lat.0.0, birth_lon.0.0, degree = 5) +
                                    factor(birth_year) + 
                                    factor(birth_mon), 
                                  data = .
                      )}

mod_uni_sun$fe50 <- fixest::feglm(
                      university ~ 
                      birth_sun_s | 
                      factor(birth_year) + 
                      factor(birth_mon) +
                      birth_lat_50k^birth_lon_50k,
                      data   = famhist,
                      family = "binomial"
                   )

mod_income_sun$fe50 <- fixest::feols(
                          income_cat ~ 
                          birth_sun_s | 
                          factor(birth_year) + 
                          factor(birth_mon) +
                          birth_lat_50k^birth_lon_50k,
                          data = famhist
                       )

mod_ea3_sun$fe50 <- fixest::feols(
                          EA3 ~ 
                          birth_sun_s | 
                          factor(birth_year) + 
                          factor(birth_mon) +
                          birth_lat_50k^birth_lon_50k,
                          data = famhist
                       )

mod_uni_sun$fe10 <- fixest::feglm(
                      university ~ 
                      birth_sun_s | 
                      factor(birth_year) +
                      factor(birth_mon) +
                      birth_lat_10k^birth_lon_10k,
                      data   = famhist,
                      family = "binomial"
                   )

mod_income_sun$fe10 <- fixest::feols(
                          income_cat ~ 
                          birth_sun_s | 
                          factor(birth_year) + 
                          factor(birth_mon) +
                          birth_lat_10k^birth_lon_10k,
                          data = famhist
                       )

mod_ea3_sun$fe10 <- fixest::feols(
                          EA3 ~ 
                          birth_sun_s | 
                          factor(birth_year) + 
                          factor(birth_mon) +
                          birth_lat_10k^birth_lon_10k,
                          data = famhist
                       )

huxreg(mod_uni_sun, 
        coefs = c("Sunshine" = "birth_sun_s"),
        tidy_args = list(se = "white")
      ) %>% 
      insert_row(after = 3, "Geography", "Polynomial", "50k dummies", "10k dummies") %>% 
      set_number_format(4, everywhere, NA) %>% 
      set_latex_float("ht")

huxreg(mod_income_sun, 
        coefs = c("Sunshine" = "birth_sun_s"),
        tidy_args = list(se = "white")
      ) %>% 
      insert_row(after = 3, "Geography", "Polynomial", "50k dummies", "10k dummies") %>% 
      set_number_format(4, everywhere, NA) %>% 
      set_latex_float("ht")

huxreg(mod_ea3_sun, 
        coefs = c("Sunshine" = "birth_sun_s"),
        tidy_args = list(se = "white")
      ) %>% 
      insert_row(after = 3, "Geography", "Polynomial", "50k dummies", "10k dummies") %>% 
      set_number_format(4, everywhere, NA) %>% 
      set_latex_float("ht")
```
