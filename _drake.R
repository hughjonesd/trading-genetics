
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


import_pcs <- function (pcs_file) {
  pcs <- read_table2(pcs_file)
  pc_names <- grep("PC", names(pcs), value = TRUE)
  pcs[pc_names] <- scale(pcs[pc_names])
  pcs
}


create_relatedness <- function (relatedness_file) {
  relatedness <- read_table2(relatedness_file, col_names = TRUE)

  relatedness %<>% mutate(
                     relation = santoku::chop(Kinship,
                                   breaks = c(
                                     1/2^(9/2), 
                                     1/2^(7/2),
                                     1/2^(5/2), 
                                     1/2^(3/2)
                                   ),
                                   labels = c(
                                     "unrelated", 
                                     "deg3",
                                     "deg2", 
                                     "parentsib",
                                      "mztwins"
                                   )
                                 ),
                     relation = case_when(
                       relation == "parentsib" & IBS0 < 0.0012  ~ "parents",
                       relation == "parentsib" & IBS0 >= 0.0012 ~ "fullsibs",
                       TRUE ~ as.character(relation)
                     )
                   )
  
  relatedness %<>% filter(relation != "unrelated") # remove the two rows with issues from the relatedness file 

  relatedness
}


make_sib_groups <- function (relatedness) {
  
  # this beast takes the pairs of siblings found in relatedness, and
  # converts them into a data frame of sibling groups, with a column for
  # the group and a column for the id
  sib_groups <- relatedness %>% 
               filter(relation %in% c("fullsibs", "mztwins")) %>% 
               select(ID1, ID2) %>% 
               igraph::graph_from_data_frame(directed = FALSE) %>% 
               igraph::max_cliques() %>%  # groups of full siblings
               purrr::map(igraph::as_ids) %>% 
               tibble::tibble() %>% 
               setNames("id") %>% 
               tibble::rowid_to_column("sib_group") %>% 
               tidyr::unnest_longer(id) %>% 
               mutate(
                 sib_group = paste0("sg", sib_group),
                 id = as.numeric(id)  # compatibility when joining
               )

  # we include mztwins because why not, but also because otherwise they
  # create non-overlapping cliques
  # a few people are in overlapping maximal cliques - typically because
  # a-b and b-c are "fullsibs" but a-c is something less, e.g. "deg2".
  # We delete these and then remove any people who have become "singletons"
  sib_groups %<>% 
                distinct(id, .keep_all = TRUE) %>% 
                group_by(sib_group) %>% 
                add_count() %>% 
                filter(n > 1) %>% 
                select(-n) %>% 
                ungroup()
  
  sib_groups
}


import_ashe_income <- function (ashe_income_file) {
  ashe_income <- readxl::read_xls(ashe_income_file, range = "A5:F475")
  
  ashe_income %<>% 
        dplyr::select(Description, Code, Median, Mean) %>% 
        mutate(across(c(Median, Mean), as.numeric)) %>% 
        rename(median_pay = Median, mean_pay = Mean)
  
  ashe_income %<>% filter(! is.na(Code))
  
  ashe_income
} 


compute_resid_scores <- function (famhist, pcs, score_names) {
  famhist <- left_join(famhist, pcs, by = c("f.eid" = "IID"))
  resid_scores <- data.frame(f.eid = famhist$f.eid)
  
  for (score_name in score_names) {
    resid_fml <- paste(score_name, "~", paste0("PC", 1:100, collapse = " + "))
    resid_score <- resid(lm(as.formula(resid_fml), famhist, 
                            na.action = na.exclude))
    resid_scores[[paste0(score_name, "_resid")]] <- resid_score
  }
  
  resid_scores
}


subset_resid_scores <- function (resid_scores_raw, famhist, score_names) {
  resid_scores <- dplyr::semi_join(resid_scores_raw, famhist, by = "f.eid")
  resid_scores[paste0(score_names, "_resid")] %<>% scale()
  
  resid_scores
}


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


add_data_to_pairs <- function (pairs_df, famhist, resid_scores, 
                                 suffix = c(".x", ".y")) {

    fhs <- famhist %>% 
                    left_join(resid_scores, by = "f.eid") %>% 
                    dplyr::select(
                      f.eid, f.6138.0.0, matches("f.6141"), female,
                      matches("_resid$"), nbro, nsis, sib_group,
                      n_sibs, birth_order, university, age_at_recruitment, YOB,
                      age_fulltime_edu, age_fte_cat, income_cat, birth_sun,
                      birth_mon, n_children, fath_age_birth, moth_age_birth,
                      first_job_pay, sr_health, illness, fluid_iq, height, 
                      f.20074.0.0, f.20075.0.0, f.699.0.0,
                      f.709.0.0, f.670.0.0, f.680.0.0, f.52.0.0, f.53.0.0, 
                      f.54.0.0, f.6139.0.0, f.6140.0.0, f.728.0.0
                    )

  pairs_df %<>% 
            left_join(fhs, by = c(eid.x = "f.eid")) %>% 
            left_join(fhs, by = c(eid.y = "f.eid"), suffix = suffix)

  pairs_df
}


make_mf_pairs <- function (mf_pairs_file, famhist, resid_scores) {
  
  mf_pairs <- readr::read_csv(mf_pairs_file)

  mf_pairs$eid.x <- mf_pairs$ID.m
  mf_pairs$eid.y <- mf_pairs$ID.f
  mf_pairs <- add_data_to_pairs(mf_pairs, famhist, resid_scores, 
                                  suffix = c(".m", ".f"))

  # removed f.670 and f.728 because they didn't help predict "having the same kid".
  mf_pairs %<>% filter(
                  f.680.0.0.f  == f.680.0.0.m,  # same rent/own status
                  f.699.0.0.f  == f.699.0.0.m,  # same length of time in hh
                  f.709.0.0.f  == f.709.0.0.m,  # same n occupants of hh
                  n_children.f == n_children.m, # same n kids
                  f.54.0.0.f   == f.54.0.0.m,   # same assessment centre 
                  f.53.0.0.f   == f.53.0.0.m,   # attended assessment same day
                  f.6141.0.0.f == 1,            # both living with spouse
                  f.6141.0.0.m == 1,            #   "     "     "    "
                  female.m != female.f,         # heterosexual couples only
                )
  orig_n <- nrow(mf_pairs)
  # remove pairs with duplications
  mf_pairs %<>% 
              group_by(ID.m) %>% 
              filter(n() == 1) %>% 
              group_by(ID.f) %>% 
              filter(n() == 1) %>% 
              ungroup()

  warning(sprintf("Removed %s pairs with multiple IDs out of %s", orig_n - nrow(mf_pairs), orig_n))
  stopifnot(all(mf_pairs$female.f == TRUE))
  
  mf_pairs$EA3.m <- mf_pairs$EA3_excl_23andMe_UK_resid.m
  mf_pairs$EA3.f <- mf_pairs$EA3_excl_23andMe_UK_resid.f

  mf_pairs$couple_id <- paste0(mf_pairs$ID.m, "_", mf_pairs$ID.f)
  
  mf_pairs
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
    create_relatedness(file_in(!! relatedness_file)),
    format = "fst_tbl"
  ),
  
  pcs = target(import_pcs(file_in(!! pcs_file)), format = "fst_tbl"),
  
  ashe_income = target(import_ashe_income(file_in(!! ashe_income_file))),
  
  famhist = target({
    famhist <- clean_famhist(famhist_raw, score_names, sib_groups, ashe_income)
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


drake_config(plan, history = FALSE)
