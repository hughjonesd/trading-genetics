
suppressPackageStartupMessages({
  requireNamespace("tibble")
  requireNamespace("tidyr")
  library(drake)
  requireNamespace("matrixStats")
  library(dplyr)
  library(readr)
  library(readxl)
  library(magrittr)
  library(purrr)
  library(forcats)
  library(santoku)
  library(tidync)
  library(abind)
})

# the project import-ukbb-data can also be downloaded from git with:
# git clone https://github.com/hughjonesd-private/import-ukbb-data.git
# then change the file path below
source("~/import-ukbb-data/import-ukbb-data.R")


not.na <- Negate(is.na)


make_parent_child <- function (relatedness, famhist_raw) {
  
  parent_child <- relatedness %>% 
                    filter(relation == "parents") %>% 
                    select(ID1, ID2)

  fh_age <- famhist_raw %>% 
              select(f.eid, age_at_reqruitment) %>% 
              rename(age = age_at_reqruitment)
  
  parent_child %<>% 
                  left_join(fh_age, by = c(ID1 = "f.eid")) %>% 
                  left_join(fh_age, by = c(ID2 = "f.eid"), suffix = c("1", "2")) %>% 
                  mutate(
                    parent_id = ifelse(age1 > age2, ID1, ID2),
                    child_id  = ifelse(age1 > age2, ID2, ID1),
                    eid.x     = parent_id, # duplicate columns so we can easily 
                    eid.y     = child_id   # add pairs data
                  ) %>% 
                  filter(
                    abs(age1 - age2) > 10    # also removes 924 "parent-children"...
                  )

  # 2 people with 3 "parents" (maybe siblings?)
  parent_child %<>% 
                  group_by(child_id) %>% 
                  add_count() %>% 
                  filter(n <= 2) %>% 
                  ungroup()

  parent_child %<>% 
                  mutate(
                    age_parent = ifelse(parent_id == ID1, age1, age2)
                  ) %>% 
                  select(-ID1, -ID2)
  
  parent_child
}


make_parent_pairs <- function (parent_child, famhist, resid_scores) {
  two_parents <- parent_child %>% 
                    select(parent_id, child_id) %>%  # remove parent-child data
                    left_join(parent_child, by = "child_id") %>% 
                    filter(
                      parent_id.x != parent_id.y,
                      parent_id.x < parent_id.y   # get rid of 1 of each 2 dupes
                    ) %>% 
                    select(
                      eid.x = parent_id.x,  
                      eid.y = parent_id.y
                    )
  
  parent_pairs <- add_data_to_pairs(two_parents, famhist, resid_scores)  
  parent_pairs %<>% 
                  filter(female.x != female.y) # probably a false positive
  
  parent_pairs$EA3.x <- parent_pairs$EA3_excl_23andMe_UK_resid.x
  parent_pairs$EA3.y <- parent_pairs$EA3_excl_23andMe_UK_resid.y
  
  parent_pairs
}


make_howe_pairs <- function (famhist, resid_scores) {
  fh_w_spouse <- famhist %>% 
                   filter(f.6141.0.0 == 1) %>% 
                   tidyr::drop_na(f.699.0.0, f.709.0.0, f.728.0.0, f.20074.0.0, 
                     f.20075.0.0, f.670.0.0, f.680.0.0, f.54.0.0)
  
  fh_w_spouse %<>% left_join(resid_scores, by = "f.eid")                
                 
  howe_pairs <- inner_join(fh_w_spouse, fh_w_spouse, 
                             by = c("f.699.0.0", "f.709.0.0",  "f.728.0.0", 
                                    "f.670.0.0", "f.680.0.0",
                                    "f.20074.0.0", "f.20075.0.0", "f.54.0.0"),
                             suffix = c(".m", ".f")
                           )
  # the above will give two copies of every match,
  # also some people will be matched with themselves

  howe_pairs %<>% filter(
                    # no individuals matched with themselves:
                    f.eid.m != f.eid.f, 
                    # the below removes (a) all-female pairs; (b) one copy out of
                    # two for each heterosexual pair:
                    ! female.m,
                    # the below removes all-male pairs:
                    female.f
                  )

  # no multiple matches
  howe_pairs %<>% 
                add_count(f.eid.m, name = "n.m") %>% 
                filter(n.m == 1) %>% 
                add_count(f.eid.f, name = "n.f") %>% 
                filter(n.f == 1)

  howe_pairs %<>% filter(
                    # no pairs with same age death of both parents:
                    ! (
                        # first three lines avoid matching people with NAs:
                        not.na(f.1807.0.0.m) & not.na(f.1807.0.0.f) &
                        not.na(f.3526.0.0.m) & not.na(f.3526.0.0.f) &
                        f.1807.0.0.m >= 0 & f.3526.0.0.m >= 0 & 
                        f.1807.0.0.m == f.1807.0.0.f & f.3526.0.0.m == f.3526.0.0.f
                      )
                  )

  # should also remove pairs with IBD > 0.1
  howe_pairs
}


filter_fake_pairs <- function (mf_pairs) {
  mf_pairs %<>% 
                mutate(
                  sameness_score = same_rent + same_time_hh + w_spouse.m + 
                                     w_spouse.f + same_n_kids + same_centre +
                                     same_day
                ) %>% 
                filter(
                  sameness_score == 5,
                  heterosexual, # heterosexual couples only
                ) 
  orig_n <- nrow(mf_pairs)
  # remove pairs with duplications
  mf_pairs %<>% 
              add_count(ID.m, name = "n.m") %>% 
              add_count(ID.m, name = "n.f") %>% 
              filter(n.m == 1, n.f == 1) %>% 
              select(-n.m, -n.f)
              
  stopifnot(all(mf_pairs$female.f == TRUE))
  
  mf_pairs
}


weight_van_alten <- function (weights_file) {
  weight_df <- readr::read_table(weights_file)
  names(weight_df) <- c("f.eid", "weights")
  
  weight_df
}

van_alten_weights_file <- "UKBSelectionWeights.tab"

plan <- drake_plan(
  score_names  = target(import_score_names(file_in(!! pgs_dir)), format = "rds"),
  
  famhist_raw = import_famhist(file_in(!! famhist_files), file_in(!! pgs_dir)), 
  
  relatedness = make_relatedness(file_in(!! relatedness_file)),
  
  pcs = import_pcs(file_in(!! pcs_file)),
  
  ashe_income = import_ashe_income(file_in(!! ashe_income_file)),
  
  van_alten_weights = weight_van_alten(file_in(!! van_alten_weights_file)),
  
  famhist = {
    famhist <- clean_famhist(famhist_raw, score_names, sib_groups)
    famhist <- add_ashe_income(famhist, ashe_income)
    famhist <- left_join(famhist, van_alten_weights, by = "f.eid")
    famhist
  },

  resid_scores_raw = compute_resid_scores(famhist_raw, pcs, score_names),
  
  resid_scores = subset_resid_scores(resid_scores_raw, famhist, score_names),
  
  parent_child = {
    parent_child <- make_parent_child(relatedness, famhist_raw)
    parent_child <- add_data_to_pairs(parent_child, famhist, resid_scores)
    parent_child
  },
  
  sib_groups = make_sib_groups(relatedness),
  
  mf_pairs_raw = make_mf_pairs(file_in(!! mf_pairs_file), famhist, resid_scores, 
                                 ashe_income),
  
  mf_pairs = filter_mf_pairs(mf_pairs_raw),
  
  mf_pairs_fake = filter_fake_pairs(mf_pairs_raw),
  
  mf_pairs_twice = make_pairs_twice(mf_pairs, suffix = c(".m", ".f")),
  
  mf_fake_twice = make_pairs_twice(mf_pairs_fake, suffix = c(".m", ".f")),
  
  parent_pairs = make_parent_pairs(parent_child, famhist, resid_scores),
  
  parent_pairs_twice = make_pairs_twice(parent_pairs),
  
  howe_pairs = make_howe_pairs(famhist, resid_scores)
)


drake_config(plan, history = FALSE, log_build_times = FALSE, format = "fst_tbl")
