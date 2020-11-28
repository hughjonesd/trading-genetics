
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

data_dir       <- "../negative-selection-data"
pgs_dir        <- file.path(data_dir, "polygenic_scores")
sun_dir        <- file.path(data_dir, "sunshine-records")

pcs_file       <- file.path(data_dir, "UKB.HM3.100PCs.40310.txt")
famhist_files <- file.path(data_dir, c(
                    "UKB.EA_pheno.coordinates.QC.csv",
                    "david.family_history.traits.out.csv",
                    "david.family_history.traits.20042020.out.csv",
                    "david.family_history.traits.05052020.out.csv",
                    "david.family_history.traits.16052020.out.csv",
                    "david.family_history.traits.18052020.out.csv",
                    "david.family_history.traits.17062020.out.csv",
                    "david.birthinfo.traits.14072020.out.csv",
                    "david.traits.03112020.out.csv",
                    "david_traits.18112020.I.csv.out.csv",
                    "david_traits.18112020.II.csv.out.csv",
                    "david_traits.18112020.III.out.csv"
                  ))
rgs_file       <- file.path(data_dir, "EA3_rgs.10052019.rgs.csv")
relatedness_file <- file.path(data_dir, "relatedness_file.txt")
mf_pairs_file  <- file.path(data_dir, "spouse_pair_info", 
                            "UKB_out.mf_pairs_rebadged.csv")
alto_file <- file.path(data_dir, "spouse_pair_info", 
                       "UKB_out.all_living_together_rebadged.csv")

ashe_income_file <- file.path(data_dir, 
                      "SOC-income", "Occupation (4) Table 14.7a   Annual pay - Gross 2007.xls") 


# utility function:
negative_to_na <- function (x) {
  x[x < 0] <- NA
  x
}


import_famhist <- function (famhist_files, pgs_dir) {

  col_types <- cols(
            .default               = col_double(), 
            geno_measurement_plate = col_skip(),
            geno_measurement_well  = col_skip(),
            `53-0.0`               = col_date(),
            `53-1.0`               = col_date(),
            `53-2.0`               = col_date()
          )
  
  # this generates a "character" spec for all 22662-x.y columns
  fake_csv_string <- paste0("22662-0.", 0:18, collapse = ",")
  fake_csv_string <- paste0(fake_csv_string, "\n")
  col_types_22662 <- spec_csv(fake_csv_string, 
                                col_types = cols(
                                  .default = col_character()
                                )
                              )
  col_types$cols <- c(col_types$cols, col_types_22662$cols)
  
  fhl <-  lapply(famhist_files, read_csv, col_types = col_types)
  
  fhl[[1]] %<>% rename(f.eid = eid)
  fhl[-1] <- purrr::map(fhl[-1], ~ {
    names(.x) <- paste0("f.", names(.x))
    names(.x) <- gsub("\\-", ".", names(.x)) 
    .x
  })
  
  famhist <- purrr::reduce(fhl, left_join, by = "f.eid")

  # remove and rename some duplicated columns
  dupes <- grep("\\.x$", names(famhist), value = TRUE)
  for (dupe in dupes) {
    dupe_y <- sub("\\.x$", ".y", dupe)
    famhist[[dupe_y]] <- NULL
  }
  famhist %<>% rename_with(
                 .fn = ~ sub("\\.x$", "", .x), 
                 .cols = all_of(dupes)
               )
  
  # only self-identified, and genetically identified, white people
  famhist %<>% filter(f.21000.0.0 == 1001, ! is.na(genetic_ethnic_grouping))
  
  for (pgs_file in list.files(pgs_dir, pattern = "csv$", full.names = TRUE)) {
    score_name <- sub(".*UKB\\.AMC\\.(.*?)\\..*", "\\1", pgs_file, perl = TRUE)
    pgs <- read_delim(pgs_file, delim = " ", col_types = "dd")
    pgs %<>% filter(FID > 0) 
    names(pgs)[2] <- score_name # instead of "SCORE"
    famhist %<>% left_join(pgs, by = c("f.eid" = "FID"))
  }
  
  return(famhist)
}


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


clean_famhist <- function (famhist, score_names, sib_groups, ashe_income) {
  # we get very few extra cases from adding f.2946.1.0 etc, and it makes calculating
  # father's year of birth more complex
  
  # -10 means less than a year, we call that 0 full years
  famhist$f.699.0.0[famhist$f.699.0.0 == -10] <- 0
  famhist$f.699.1.0[famhist$f.699.1.0 == -10] <- 0
  famhist$f.699.2.0[famhist$f.699.2.0 == -10] <- 0
  
  # remove negatives
  famhist %<>% mutate(across(
      c(age_fulltime_edu, starts_with(c(
        "f.2946", "f.1845", "f.2754", "f.738",  "f.2764", "f.2405", "f.2734",
        "f.2149", "f.1873", "f.1883", "f.2784", "f.2794", "f.709", 
        "f.699", "f.3872", "f.728", "f.670", "f.680",
        "f.5057", "birth_lon", "birth_lat"
      ))), 
      negative_to_na
    )
  )
  # -7 means "never went to school"
  famhist$f.6138.0.0[famhist$f.6138.0.0 == -3] <- NA
  
  
  famhist$age_at_recruitment <- famhist$f.21022.0.0
  # questionnaire sex:
  famhist$female <- famhist$f.31.0.0 == 0 
  famhist$birth_year <- famhist$f.34.0.0
  famhist$birth_mon <- famhist$f.52.0.0
  
  # "Field 845 was collected from all participants except those who indicated 
  # they have a College or University degree, as defined by their answers to 
  # Field 6138". So, we impute this to be 21.
  famhist$age_fulltime_edu[is.na(famhist$age_fulltime_edu) & famhist$edu_qual.0.0 == 1] <- 21
  
  # TODO: ask Abdel how edu_qual was made up. Seems to be a "max" of all the
  # individual edu_qual.x.y's; these are just f.6138.x.y.
  # TODO: that doesn't work!!! Warn Abdel.
  # Currently using edu_qual.0.0
  famhist$university <- famhist$edu_qual.0.0 == 1
  famhist$income_cat <- famhist$f.738.0.0
  
  famhist$n_children <- pmax(famhist$f.2405.0.0, famhist$f.2405.1.0,
                             famhist$f.2405.2.0, famhist$f.2734.0.0, famhist$f.2734.1.0, 
                             famhist$f.2734.2.0, 
                             na.rm = TRUE
  )
  
  famhist$n_in_household <- famhist$f.709.0.0
  
  famhist$with_partner   <- famhist$f.6141.0.0 == 1
  # Many NAs, almost all from people living alone i.e. f.709 == 1
  famhist$with_partner[famhist$n_in_household == 1] <- FALSE
  famhist$with_partner[famhist$f.6141.0.0 == -3] <- NA
  

  famhist$age_fte_cat <- santoku::chop(famhist$age_fulltime_edu,
                                       c(16, 18),
                                       c("< 16", "16-18", "> 18"))

  # -7 means never went to school. We recode to 0 for simpliciy
  famhist$edu_qual[famhist$edu_qual == -7] <- 0
  famhist$edu_qual[famhist$edu_qual == -3] <- NA
  
  # we use pmax, assuming that people *can* have given birth for the first
  # time in between surveys.
  famhist$age_flb <-  pmax(
                        famhist$f.3872.0.0, famhist$f.3872.1.0, famhist$f.3872.2.0,
                        famhist$f.2754.0.0, famhist$f.2754.1.0, famhist$f.2754.2.0,
                        na.rm = TRUE
                      )
  famhist$age_flb_cat <- santoku::chop_equally(famhist$age_flb, 3, 
                                               labels = lbl_discrete("-"))
  
  # full brothers and sisters
  famhist$nbro <- pmax(famhist$f.1873.0.0, famhist$f.1873.1.0, 
                       famhist$f.1873.2.0, na.rm = TRUE)
  famhist$nsis <- pmax(famhist$f.1883.0.0, famhist$f.1883.1.0, 
                       famhist$f.1883.2.0, na.rm = TRUE)
  famhist$n_sibs <- famhist$nbro + famhist$nsis + 1
  # a few people give varying answers, we assume median is fine.
  # including later answers picks up c. 10K extra people.
  # some people have a non-integer median; we round down.
  famhist$birth_order <- floor(matrixStats::rowMedians(
    as.matrix(famhist[c("f.5057.0.0", "f.5057.1.0", "f.5057.2.0")]),
    na.rm = TRUE
  ))
  famhist$birth_order[famhist$n_sibs == 1] <- 0
  # number of older siblings, plus one:
  famhist$birth_order <- famhist$birth_order + 1
  # TODO: how does this come about??? Stupid answers?
  famhist$birth_order[famhist$birth_order > famhist$n_sibs] <- NA
  # TODO: why is birth_order often NaN?
  # TODO: why is f.5057 often NA when n_sibs == 1? And why often NA in general
  
  # TODO: minimum fath_age_birth is 3, moth_age_birth is 0... Why?
  famhist$fath_age <- famhist$f.2946.0.0
  famhist$moth_age <- famhist$f.1845.0.0
  famhist$fath_age_birth <- famhist$fath_age - famhist$age_at_recruitment
  famhist$moth_age_birth <- famhist$moth_age - famhist$age_at_recruitment

  # TODO: get f.20191, it is the touchscreen equivalent and will have more data
  famhist$fluid_iq <- rowMeans(famhist %>% select(starts_with("f.20016")), na.rm = TRUE)
  
  # first measurement has almost everyone
  # 1 = excellent, 2 = good, 3 = fair, 4 = poor
  # I reverse-code
  famhist$sr_health <- -1 * negative_to_na(famhist$f.2178.0.0)
  # "longstanding illness, disability or infirmity". 1 = TRUE
  famhist$illness   <- negative_to_na(famhist$f.2188.0.0)
  famhist$height <- famhist$f.50.0.0
  famhist$bmi <- famhist$f.21001.0.0
  
  famhist$num_jobs <- famhist$f.22599.0.0
  # job codes are f.22601.0.x
  # SOC2000 job codes are f.22617.0.x
  # start years   f.22602.0.x
  # end years     f.22603.0.x, which we don't have
  

  famhist[score_names] <- scale(famhist[score_names])

  # TODO: ask Abdel for job code f.20277
  # famhist %<>% left_join(ashe_income, by = c("f.20277" = "Code"))
  
  famhist %<>% 
        mutate(f.22617.0.0 = as.character(f.22617.0.0)) %>% 
        left_join(ashe_income, by = c("f.22617.0.0" = "Code")) %>% 
        select(-Description, -mean_pay) %>% 
        mutate(first_job_pay = median_pay/1000) %>% 
        select(-median_pay)
  
  famhist %<>% left_join(sib_groups, by = c("f.eid" = "id"))
  
  return(famhist)
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
  score_names  = {
    score_names <- sub(
      ".*UKB\\.AMC\\.(.*?)\\..*", 
      "\\1", 
      list.files(file_in(!! pgs_dir), pattern = "csv$"), 
      perl = TRUE
    )
    setNames(score_names, score_names)
  },
  
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
