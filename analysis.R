
# setup ====
# 
library(drake)
library(dplyr)
library(AER)
library(ggplot2)

drake::loadd(mf_pairs)
drake::loadd(famhist)


# own education and spouse's PGS ====
famhist %>% 
      group_by(f.6138.0.0) %>% 
      summarize(
        across(c(EA3_excl_23andMe_UK, height_combined), mean, na.rm = TRUE)
      )


# own birth month and education ====
tbl1 <- table(famhist$f.52.0.0, famhist$f.6138.0.0)
chisq.test(tbl1)
ptbl1 <- proportions(tbl1, 1)
round(ptbl1 * 100, 2)


# narrow down to "university"
proportions(table(famhist$university, famhist$f.52.0.0), 2)
# it's not much to work with... a 2 percentage point difference!


# controlling for age and gender ====
# 
mod_uni_bm <- lm(university ~ factor(f.52.0.0) + age_at_recruitment + 
             I(age_at_recruitment^2) + female, famhist)
# the biggest effect is about 1/3 the effect of sex.


# own birth month and own genetics ====
# 
famhist %>% 
      group_by(f.52.0.0) %>% 
      summarize(
        across(c(university, EA3_excl_23andMe_UK), mean, na.rm = TRUE)
      ) %>% 
      ggplot(aes(f.52.0.0)) + 
        geom_line(aes(y = university), color = "red") +
        geom_line(aes(y = EA3_excl_23andMe_UK), color = "blue") +
        labs(x = "birth month", y = "EA3 (blue)/Prob. university (red)")

# yes, there are some significant differences!
mod_ea3_bm <- lm(EA3_excl_23andMe_UK ~ factor(f.52.0.0), famhist)
# note that these peak in july, whereas university peaks in sept/oct,
# suggesting that genetics and environment push in *opposite* directions


# the below is a bit pointless because of unmeasured genetics:
mod_uni_bm2 <- lm(university ~ factor(f.52.0.0) + EA3_excl_23andMe_UK + 
                  poly(age_at_recruitment, 5) + female, 
                  famhist)

# own birth order and university ====

famhist %>% 
      group_by(n_sibs, n_older_sibs) %>% 
      filter(! is.na(n_older_sibs), n_sibs <= 8) %>% 
      summarize(mean_uni = mean(university, na.rm = TRUE)) %>% 
      ggplot(aes(n_older_sibs, mean_uni, colour = factor(n_sibs), group = n_sibs)) + 
      geom_line() + geom_point()
# note that only children do very badly! Maybe poverty?

# we only want to look at effect of n_older_sibs within each sibling group size
mod_uni_bo <- lm(university ~ factor(n_sibs) + n_older_sibs:factor(n_sibs) +
                 age_at_recruitment + I(age_at_recruitment^2) + female, 
                 famhist, n_sibs >= 2 & n_sibs <= 8)
# this is highly significant and effect sizes are rather chunky :-)

# unfortunately, there are also "effects" on EA3 for birth order, especially
# among larger families. 
famhist %>% 
      group_by(n_sibs, n_older_sibs) %>% 
      filter(! is.na(n_older_sibs), n_sibs <= 8) %>% 
      summarize(mean_EA3 = mean(EA3_excl_23andMe_UK, na.rm = TRUE)) %>% 
      ggplot(aes(n_older_sibs, mean_EA3, colour = factor(n_sibs), group = n_sibs)) + 
      geom_line() + geom_point()
# Note that again, EA3 and university go the opposite way!
# This isn't a miscoding:
famhist %>% 
      group_by(university) %>%
      summarize(mean_EA3 = mean(EA3_excl_23andMe_UK, na.rm = TRUE))

mod_ea3_bo <- lm(EA3_excl_23andMe_UK ~ factor(n_sibs) + n_older_sibs:factor(n_sibs),
                 famhist, n_sibs > 1 & n_sibs <= 8 & ! is.nan(n_older_sibs))
# Maybe if your first kid seems dumb, you decide to have more. So
# then, among e.g. 3-child families, the first kid will be dumber than
# average (conditioning on the parents), while the second and subsequent 
# children will be average
# NB also: weird values for n_older_sibs, needs investigating
# NB that there are also very big effect sizes for large families - as you
# might expect given your other paper


# couple regressions ====



mod_ea3f_iv <- AER::ivreg(EA3_excl_23andMe_UK.f ~ university.m + factor(n_sibs.m) 
               + EA3_excl_23andMe_UK.m | 
               factor(f.52.0.0.m) + factor(n_sibs.m) + 
               n_older_sibs.m:factor(n_sibs.m) + EA3_excl_23andMe_UK.m, 
               data = mf_pairs, subset = n_sibs.m > 1 & n_sibs.m <= 8 & ! is.nan(n_older_sibs.m))


mod_ea3m_iv <- AER::ivreg(EA3_excl_23andMe_UK.m ~ university.f + factor(n_sibs.f) 
               + EA3_excl_23andMe_UK.f | 
               factor(f.52.0.0.f) + factor(n_sibs.f) + 
               n_older_sibs.f:factor(n_sibs.f) + EA3_excl_23andMe_UK.f, 
               data = mf_pairs, subset = n_sibs.f > 1 & n_sibs.f <= 8)

