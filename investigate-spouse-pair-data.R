
library(readr)
library(tidyr)
library(dplyr)


  
drake::loadd(famhist)
famhist_tmp <- famhist %>% 
                left_join(resid_scores, by = "f.eid") %>% 
                dplyr::select(
                  f.eid, f.6138.0.0, matches("f.6141"), f.52.0.0, female,
                  matches("_resid$"),
                  n_sibs, n_older_sibs, university, age_at_recruitment, YOB,
                  age_fulltime_edu, age_fte_cat, income_cat, birth_sun,
                  birth_mon, n_children, fath_age_birth, moth_age_birth,
                  first_job_pay, sr_health, illness, fluid_iq
                )

mf_pairs %<>% 
  left_join(famhist_tmp, by = c("ID.m" = "f.eid")) %>% 
  left_join(famhist_tmp, by = c("ID.f" = "f.eid"), suffix = c(".m", ".f"))


mf_pairs %<>% 
  mutate(
          ever_lived_with_spouse = (
                                   (f.6141.0.0.m == 1 & f.6141.0.0.f == 1) 
                                 )
        ) %>% 
  filter(ever_lived_with_spouse) %>% 
  select(-ever_lived_with_spouse)

mf_pairs$couple_id <- paste(mf_pairs$ID.m, mf_pairs$ID.f, sep = "_")
mf_pairs %<>% distinct(couple_id, .keep_all = TRUE)

mf_pairs %<>%
  group_by(ID.m) %>% 
  mutate(n.m = n()) %>% 
  group_by(ID.f, .drop = TRUE) %>% 
  mutate(n.f = n()) %>% 
  ungroup()


table(mf_pairs$n.m)
table(mf_pairs$n.f)
# the majority of people (90% ish) live with more than just one person at time 0
# in the dataset, AND say they are living with their spouse!
# this is weird, because most households surely don't have multiple spouse
# pairs.

dupes_m <- mf_pairs %>% filter(n.m > 1) %>% arrange(ID.m)
dupes_f<- mf_pairs %>% filter(n.f > 1) %>% arrange(ID.f)



# individual-level list of people and coordinate-ids
# there are many duplicate IDs, since people move.

alto_file <- "~/negative-selection-data/spouse_pair_info/UKB_out.all_living_together_rebadged.csv"
alto <- read_csv(alto_file)
alto %<>% 
  group_by(ID) %>% 
  mutate(mpsc = max(people_at_same_coordinate)) 

# we have indeed many people who live with many people at one coordinate,
# at some point. Mean is 6.5
table(alto$mpsc)


# IDEA
# pick out every set of people who are living at the same coordinate at UKB
# assessment 1.
# count the number within each set who are living with their spouse at that time.
# if that number is exactly 2, include both people who are living with their spouse.
# otherwise don't (not a couple, or multiple couples and can't identify)


alto %<>% 
  filter(grepl("1", UKB_assessment.nrs)) %>% 
  left_join(famhist_tmp %>% select(f.eid, f.6141.0.0), by = c("ID" = "f.eid")) %>% 
  group_by(coordinate_ID) %>% 
  mutate(n_with_spouse = sum(f.6141.0.0 == 1)) 

# seems weird that so many people are living in hhs with multiple spouse pairs:
table(alto$n_with_spouse)

spouse_pairs <- alto %>% 
                  filter(n_with_spouse == 2, f.6141.0.0 == 1) %>% 
                  select(ID, coordinate_ID) %>% 
                  arrange(coordinate_ID)

# check we only have couples
ns <- spouse_pairs %>% group_by(coordinate_ID) %>% count() %>% pull(n)
stopifnot(all(ns == 2))

spouse_pairs <- matrix(spouse_pairs$ID, ncol = 2, byrow = TRUE)
spouse_pairs %<>% tibble::as.tibble() %>% rename(ID1 = V1, ID2 = V2)         

real_mf_pairs <- spouse_pairs %>% 
                   left_join(famhist_tmp, by = c("ID1" = "f.eid")) %>% 
                   left_join(famhist_tmp, by = c("ID2" = "f.eid"))


                     
