
library(drake)
library(dplyr)
library(readr)
library(magrittr)
library(purrr)
library(forcats)
library(santoku)

data_dir       <- "../negative-selection-data"
pgs_dir        <- file.path(data_dir, "polygenic_scores")
pcs_file       <- file.path(data_dir, "UKB.HM3.100PCs.40310.txt")
famhist_file   <- file.path(data_dir, "david.family_history.traits.out.csv")
famhist2_file  <- file.path(data_dir, "david.family_history.traits.20042020.out.csv")
famhist3_file  <- file.path(data_dir, "david.family_history.traits.05052020.out.csv")
famhist4_file  <- file.path(data_dir, "david.family_history.traits.16052020.out.csv")
famhist5_file  <- file.path(data_dir, "david.family_history.traits.18052020.out.csv")
famhist6_file  <- file.path(data_dir, "david.family_history.traits.17062020.out.csv")
famhist7_file  <- file.path(data_dir, "david.birthinfo.traits.14072020.out.csv")
rgs_file       <- file.path(data_dir, "EA3_rgs.10052019.rgs.csv")
mf_pairs_file  <- file.path(data_dir, "spouse_pair_info", 
                            "UKB_out.mf_pairs_rebadged.csv")


# utility function:
negative_to_na <- function (x) {
  x[x < 0] <- NA
  x
}


make_famhist <- function (
  famhist_file,
  famhist2_file,
  famhist3_file,
  famhist4_file,
  famhist5_file,
  famhist6_file,
  famhist7_file,
  pgs_dir
) {

  fhl <- list()
  fhl[[1]] <- read_csv(famhist_file, col_types = strrep("d", 40))
  fhl[[2]] <- read_csv(famhist2_file, col_types = strrep("d", 17))
  fhl[[3]] <- read_csv(famhist3_file, col_types = strrep("d", 4))
  fhl[[4]] <- read_csv(famhist4_file, col_types = strrep("d", 33))
  fhl[[5]] <- read_csv(famhist5_file, col_types = strrep("d", 4))
  fhl[[6]] <- read_csv(famhist6_file, col_types = strrep("d", 22))
  fhl[[7]] <- read_csv(famhist7_file, col_types = strrep("d", 3))
  
  fhl <- purrr::map(fhl, ~ {
    names(.x) <- paste0("f.", names(.x))
    names(.x) <- gsub("\\-", ".", names(.x)) 
    .x
  })
  
  famhist <- purrr::reduce(fhl, left_join, by = "f.eid")

  # only "genetic" whites
  # TODO: uncomment
  # famhist %<>% filter(! is.na(genetic_ethnic_grouping))
  
  for (pgs_file in list.files(pgs_dir, pattern = "csv$", full.names = TRUE)) {
    score_name <- sub(".*UKB\\.AMC\\.(.*?)\\..*", "\\1", pgs_file, perl = TRUE)
    pgs <- read_delim(pgs_file, delim = " ", col_types = "dd")
    pgs %<>% filter(FID > 0) 
    names(pgs)[2] <- score_name # instead of "SCORE"
    famhist %<>% left_join(pgs, by = c("f.eid" = "FID"))
  }
  
  return(famhist)
}


clean_famhist <- function (famhist, score_names) {
  # we get very few extra cases from adding f.2946.1.0 etc, and it makes calculating
  # father's year of birth more complex
  
  names(famhist) <- sub("age_at_reqruitment", "age_at_recruitment", 
                        names(famhist))
  
# TODO: temporary fix
  famhist$age_fulltime_edu <- NA_real_
  # remove negatives
  famhist %<>% mutate(across(
    c(age_fulltime_edu, starts_with(c(
      "f.2946", "f.1845", "f.2754", "f.738",  "f.2764", "f.2405", "f.2734",
      "f.2149", "f.1873", "f.1883", "f.2784", "f.2794", "f.709",  "f.3872",
      "f.5057", "f.6138"
    ))), 
    negative_to_na
  )
  )
  
  # "Field 845 was collected from all participants except those who indicated 
  # they have a College or University degree, as defined by their answers to 
  # Field 6138". So, we impute this to be 21.
  famhist$age_fulltime_edu[is.na(famhist$age_fulltime_edu) & famhist$edu_qual == 1] <- 21
  
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
  
# TODO: uncomment
  # famhist$age_fte_cat <- santoku::chop(famhist$age_fulltime_edu, 
  #                                      c(16, 18), 
  #                                      c("< 16", "16-18", "> 18"))
  # 
  # # -7 means never went to school. We recode to 0 for simpliciy
  # famhist$edu_qual[famhist$edu_qual == -7] <- 0
  # famhist$edu_qual[famhist$edu_qual == -3] <- NA
  
  # we use pmax, assuming that people *can* have given birth for the first
  # time in between surveys.
  famhist$age_flb <- pmax(
    famhist$f.3872.0.0, famhist$f.3872.1.0, famhist$f.3872.2.0,
    famhist$f.2754.0.0, famhist$f.2754.1.0, famhist$f.2754.2.0,
    na.rm = TRUE
  )
  famhist$age_flb_cat <- santoku::chop_equally(famhist$age_flb, 3, 
                                               labels = lbl_discrete("-"))
  
  famhist[score_names] <- scale(famhist[score_names])
  
  return(famhist)
}


make_mf_pairs <- function (mf_pairs_file, famhist) {
  mf_pairs <- read_csv(mf_pairs_file, col_types = "dddccccc")
  famhist_tmp <- famhist %>% dplyr::select(! (starts_with("f.") | 
                                                starts_with("PC") | starts_with("home") |
                                                starts_with("assessment") | starts_with("birth_")), f.eid)
  mf_pairs %<>% 
    left_join(famhist_tmp, by = c("ID.m" = "f.eid")) %>% 
    left_join(famhist_tmp, by = c("ID.f" = "f.eid"), suffix = c(".m", ".f"))
  
  mf_pairs
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
    make_famhist( 
      file_in(!! famhist_file),
      file_in(!! famhist2_file),
      file_in(!! famhist3_file),
      file_in(!! famhist4_file),
      file_in(!! famhist5_file),
      file_in(!! famhist6_file),
      file_in(!! famhist7_file),
      file_in(!! pgs_dir)
    ), 
    format = "fst"
  ), 
  
  famhist = target(
    clean_famhist(famhist_raw, score_names),
    format = "fst"
  ),
  
  pairs = target(
    make_mf_pairs(file_in(!! mf_pairs_file), famhist),
    format = "fst"
  )
)


drake_config(plan, history = FALSE)
