
suppressPackageStartupMessages({
  requireNamespace("tibble")
  requireNamespace("tidyr")
  library(drake)
  requireNamespace("matrixStats")
  library(raster)
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
                    child_id  = ifelse(age1 > age2, ID2, ID1)
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

  parent_child
}


make_parent_pairs <- function (parent_child, famhist, resid_scores) {
  two_parents <- parent_child %>% 
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


make_pairs_twice <- function (pairs_df) {
  pairs_df_rebadged <- pairs_df
  names(pairs_df_rebadged) <- sub("\\.x$", "\\.tmp", names(pairs_df_rebadged))
  names(pairs_df_rebadged) <- sub("\\.y$", "\\.x", names(pairs_df_rebadged))
  names(pairs_df_rebadged) <- sub("\\.tmp$", "\\.y", names(pairs_df_rebadged))
  
  pairs_twice <- bind_rows(pairs_df, pairs_df_rebadged)
  pairs_twice$x <- ifelse(pairs_twice$female.x, "Female", "Male")
  
  pairs_twice
}

compute_sunshine <- function (famhist, sun_dir) {
  years <- 1941:1970
  sun_files <- list.files(sun_dir, pattern = ".nc$", full.names = TRUE)
  famhist$birth_sun <- famhist$trim_2_sun <- famhist$trim_3_sun <- NA_real_
  
  for (year in years) {
    prev_ras <- if (exists("sun_ras")) {
                    sun_ras
                  } else {
                    prev_file <- grep(paste0("mon_", year - 1), sun_files, value = TRUE)
                    raster::stack(prev_file)
                  }
    sun_file <- grep(paste0("mon_", year), sun_files, value = TRUE)
    sun_ras <- raster::stack(sun_file)
    
    sun_array <- as.array(sun_ras)
    prev_array <- as.array(prev_ras)
    # a very high value used for NA
    sun_array[sun_array == max(sun_array)] <- NA
    prev_array[prev_array == max(prev_array)] <- NA
    sun_array <- abind(prev_array, sun_array, along = 3)
    
    # the second dimension (columns) of sun_array is latitude (y)
    # the first dimension (rows) is latitude (x)
    fh_yr <- famhist %>% filter(birth_year == year)
    fh_yr$col <- raster::colFromX(sun_ras, fh_yr$birth_lon.0.0)
    fh_yr$row <- raster::rowFromY(sun_ras, fh_yr$birth_lat.0.0)
    fh_yr %<>% filter(! is.na(row), ! is.na(col), ! is.na(birth_mon))
    fh_yr$trim_2_sun <- fh_yr$trim_3_sun <- 0 
    for (mon in 1:6) {
      # index into data. Current year starts at 13.
      fh_yr$month <- fh_yr$birth_mon - mon + 12
      indices <- fh_yr %>% 
                    dplyr::select(row, col, month) %>% 
                    as.matrix()  
      month_sun <- sun_array[indices]
      var_name <- if (mon <= 3) "trim_2_sun" else "trim_3_sun"
      fh_yr[[var_name]] <- fh_yr[[var_name]] + month_sun
    }
    fh_yr$birth_sun <- fh_yr$trim_2_sun + fh_yr$trim_3_sun
    famhist %<>% rows_update(fh_yr %>% 
                               dplyr::select(f.eid, birth_sun, trim_2_sun, trim_3_sun), 
                               by = "f.eid"
                             )
  }
  
  return(famhist)
}


plan <- drake_plan(
  score_names  = import_score_names(file_in(!! pgs_dir)),
  
  famhist_raw  = target(
    import_famhist(file_in(!! famhist_files), file_in(!! pgs_dir)), 
    format = "fst_tbl"
  ), 
  
  relatedness = target(
    make_relatedness(file_in(!! relatedness_file)),
    format = "fst_tbl"
  ),
  
  pcs = target(import_pcs(file_in(!! pcs_file)), format = "fst_tbl"),
  
  ashe_income = target(import_ashe_income(file_in(!! ashe_income_file))),
  
  famhist = target({
    famhist <- clean_famhist(famhist_raw, score_names, sib_groups)
    famhist <- add_ashe_income(famhist, ashe_income)
    famhist <- compute_sunshine(famhist, file_in(!! sun_dir))
    famhist
    },
    format = "fst_tbl"
  ),
  
  resid_scores_raw = target(
    compute_resid_scores(famhist_raw, pcs, score_names),
    format = "fst_tbl"
  ),
  
  resid_scores = target(
    subset_resid_scores(resid_scores_raw, famhist, score_names),
    format = "fst_tbl"
  ),
  
  parent_child = target(
    make_parent_child(relatedness, famhist_raw),
    format = "fst_tbl"
  ),
  
  sib_groups = target(make_sib_groups(relatedness), format = "fst_tbl"),
  
  mf_pairs = target(
    make_mf_pairs(file_in(!! mf_pairs_file), famhist, resid_scores),
    format = "fst_tbl"
  ),
  
  mf_pairs_twice = target({
      names(mf_pairs) <- sub("\\.m$", ".x", names(mf_pairs))
      names(mf_pairs) <- sub("\\.f$", ".y", names(mf_pairs))
      make_pairs_twice(mf_pairs)
    },
    format = "fst_tbl"
  ),
  
  parent_pairs = target(
    make_parent_pairs(parent_child, famhist, resid_scores),
    format = "fst_tbl"
  ),
  
  parent_pairs_twice = target(
    make_pairs_twice(parent_pairs),
    format = "fst_tbl"
  )
)


drake_config(plan, history = FALSE, log_build_times = FALSE)
