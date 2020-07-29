
suppressPackageStartupMessages({
  library(drake)
  library(matrixStats)
  loadNamespace("raster")
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


# utility function:
negative_to_na <- function (x) {
  x[x < 0] <- NA
  x
}


import_famhist <- function (famhist_files, pgs_dir) {

  fhl <- lapply(famhist_files, read_csv, col_types = cols(
            .default = "d", 
            geno_measurement_plate = col_skip(),
            geno_measurement_well  = col_skip()
          ))
  
  fhl[[1]] %<>% rename(f.eid = eid)
  fhl[-1] <- purrr::map(fhl[-1], ~ {
    names(.x) <- paste0("f.", names(.x))
    names(.x) <- gsub("\\-", ".", names(.x)) 
    .x
  })
  
  famhist <- purrr::reduce(fhl, left_join, by = "f.eid")

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


import_ashe_income <- function (ashe_income_file) {
  ashe_income <- readxl::read_xls(ashe_income_file, range = "A5:F475")
  ashe_income %<>% 
        select(Description, Code, Median, Mean) %>% 
        mutate(across(c(Median, Mean), as.numeric)) %>% 
        rename(median_pay = Median, mean_pay = Mean)
  
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


clean_famhist <- function (famhist, score_names, ashe_income) {
  # we get very few extra cases from adding f.2946.1.0 etc, and it makes calculating
  # father's year of birth more complex
  
  # remove negatives
  famhist %<>% mutate(across(
      c(age_fulltime_edu, starts_with(c(
        "f.2946", "f.1845", "f.2754", "f.738",  "f.2764", "f.2405", "f.2734",
        "f.2149", "f.1873", "f.1883", "f.2784", "f.2794", "f.709",  "f.3872",
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
  famhist$age_fulltime_edu[is.na(famhist$age_fulltime_edu) & famhist$edu_qual == 1] <- 21
  
  # TODO: ask Abdel how edu_qual was made up. Seems to be a "max" of all the
  # individual edu_qual.x.y's; these are just f.6138.x.y
  famhist$university <- famhist$edu_qual == 1
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
  famhist$n_older_sibs <- floor(matrixStats::rowMedians(
    as.matrix(famhist[c("f.5057.0.0", "f.5057.1.0", "f.5057.2.0")]),
    na.rm = TRUE
  ))
  famhist$n_older_sibs[famhist$n_sibs == 1] <- 0
  # TODO: how does this come about??? Stupid answers?
  famhist$n_older_sibs[famhist$n_older_sibs >= famhist$n_sibs] <- NA
  # TODO: why is n_older_sibs often NaN?
  # TODO: why is f.5057 often NA when n_sibs == 1? And why often NA in general
  
  famhist[score_names] <- scale(famhist[score_names])

  # TODO: ask Abdel for job code f.20277
  # famhist %<>% left_join(ashe_income, by = c("f.20277" = "Code"))
  
  return(famhist)
}


subset_resid_scores <- function (resid_scores_raw, famhist, score_names) {
  resid_scores <- dplyr::semi_join(resid_scores_raw, famhist, by = "f.eid")
  resid_scores[paste0(score_names, "_resid")] %<>% scale()
  
  resid_scores
}


compute_sunshine <- function (famhist, sun_dir) {
  years <- 1941:1970
  sun_files <- list.files(sun_dir, pattern = ".nc$", full.names = TRUE)
  famhist$birth_sun <- NA_real_
  
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
    fh_yr$birth_sun <- 0 
    for (mon in 1:6) {
      # index into data. Current year starts at 13.
      fh_yr$month <- fh_yr$birth_mon - mon + 12
      indices <- fh_yr %>% 
                    dplyr::select(row, col, month) %>% 
                    as.matrix()  
      month_sun <- sun_array[indices]
      fh_yr$birth_sun <- fh_yr$birth_sun + month_sun
    }
    famhist %<>% rows_update(fh_yr %>% dplyr::select(f.eid, birth_sun), by = "f.eid")
  }
  
  return(famhist)
}


make_mf_pairs <- function (mf_pairs_file, famhist, resid_scores) {
  mf_pairs <- read_csv(mf_pairs_file, col_types = "dddccccc")
  
  famhist_tmp <- famhist %>% 
                  left_join(resid_scores, by = "f.eid") %>% 
                  dplyr::select(
                    f.eid, f.6138.0.0, f.52.0.0, matches("_resid$"),
                    n_sibs, n_older_sibs, university, age_at_recruitment,
                    age_fulltime_edu, age_fte_cat, income_cat, birth_sun
                  )
  mf_pairs %<>% 
    left_join(famhist_tmp, by = c("ID.m" = "f.eid")) %>% 
    left_join(famhist_tmp, by = c("ID.f" = "f.eid"), suffix = c(".m", ".f"))
  
  mf_pairs$couple_id <- paste(mf_pairs$ID.m, mf_pairs$ID.f, sep = "_")
  mf_pairs$EA3.m <- mf_pairs$EA3_excl_23andMe_UK_resid.m
  mf_pairs$EA3.f <- mf_pairs$EA3_excl_23andMe_UK_resid.f

  mf_pairs
}


make_mf_pairs_twice <- function (mf_pairs) {
  mf_pairs$x <- "Male"
  
  mf_pairs_rebadged <- mf_pairs
  mf_pairs_rebadged$x <- "Female"
  
  names(mf_pairs_rebadged) <- sub("\\.f$", ".mxxx", names(mf_pairs_rebadged))
  names(mf_pairs_rebadged) <- sub("\\.m$", ".f", names(mf_pairs_rebadged))
  names(mf_pairs_rebadged) <- sub("\\.mxxx$", ".m", names(mf_pairs_rebadged))
  
  mf_pairs_twice <- bind_rows(mf_pairs, mf_pairs_rebadged)
  names(mf_pairs_twice) <- sub("\\.m", ".x", names(mf_pairs_twice))
  names(mf_pairs_twice) <- sub("\\.f", ".y", names(mf_pairs_twice))
  
  mf_pairs_twice
}


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
                    "david.birthinfo.traits.14072020.out.csv"
                  ))
rgs_file       <- file.path(data_dir, "EA3_rgs.10052019.rgs.csv")
mf_pairs_file  <- file.path(data_dir, "spouse_pair_info", 
                            "UKB_out.mf_pairs_rebadged.csv")
ashe_income_file <- file.path(data_dir, 
                      "SOC-income", "Occupation (4) Table 14.7b   Annual pay - Gross 2007 CV.xls") 


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
    format = "fst"
  ), 
  
  pcs = target(import_pcs(file_in(!! pcs_file)), format = "fst"),
  
  ashe_income = target(import_ashe_income(file_in(!! ashe_income_file))),
  
  famhist = target({
    famhist <- clean_famhist(famhist_raw, score_names, ashe_income)
    famhist <- compute_sunshine(famhist, file_in(!! sun_dir))
    famhist
    },
    format = "fst"
  ),
  
  resid_scores_raw = target(
    compute_resid_scores(famhist_raw, pcs, score_names),
    format = "fst"
  ),
  
  resid_scores = target(
    subset_resid_scores(resid_scores_raw, famhist, score_names),
    format = "fst"
  ),
  
  mf_pairs = target(
    make_mf_pairs(file_in(!! mf_pairs_file), famhist, resid_scores),
    format = "fst"
  ),
  
  mf_pairs_twice = target(
    make_mf_pairs_twice(mf_pairs),
    format = "fst"
  )
)


drake_config(plan, history = FALSE)
