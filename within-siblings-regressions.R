
# run setup chunk from trading-genetics.Rmd


fml_bo_psea_fixef <- list()

fml_bo_psea_fixef[[1]] <- EA3.y ~ birth_order.x | n_sibs.x + sib_group.x
fml_bo_psea_fixef[[2]] <- EA3.y ~ birth_order.x + EA3.x |
                            n_sibs.x + sib_group.x
  
fml_bo_psea_fixef[[3]] <- EA3.y ~ birth_order.x + EA3.x + 
                            moth_age_birth.x | n_sibs.x + sib_group.x

# not a fixed effects model
# fml_bo_psea_fixef[[4]] <- EA3.y ~ birth_order.x + EA3.x + 
#                             moth_age_birth.x + moth_afb.x | n_sibs.x

fml_bo_psea_fixef %<>% map(Formula::as.Formula)


mod_bo_psea_fixef <- lapply(fml_bo_psea_fixef, fixest::feols, 
                 data = mf_pairs_twice
               )


huxreg(mod_bo_psea_fixef, 
         coefs = c(
           "Birth order"           = "birth_order.x", 
           "Own PSEA"              = "EA3.x",
           "Mother's age at birth" = "moth_age_birth.x"
         ),
         statistics = c(
           "N"  = "nobs", 
           "R2" = "r.squared"
          ),
         stars = my_stars,
         note = "{stars}. Robust SEs in brackets.",
         tidy_args = list(
           conf.int = FALSE, 
           se = "hetero"
           #cluster  = list(mf_pairs_reg$sib_group.x, mf_pairs_reg$n_sibs.x)
         )
       ) %>% 
       insert_row(after = 7, "Family size dummies", rep("Yes", 3)) %>% 
       insert_row(after = 8, "Sibling group dummies", "Yes", "Yes", "Yes") %>% 
       set_bottom_border(8, everywhere, 0) %>% 
       set_number_format(2:7, -1, 4) %>% 
       set_width(0.5) %>% 
       set_position("left") %>% 
       set_caption("Regressions of spouse PSEA on birth order: within-family")

