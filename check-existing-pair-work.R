

# mf_pairs plan

# for each period, find all pairs of ids who are living together
# select all pairs who:
# - both living with spouse 6140
# - have been living at the address for the same length of time 699
# - same n house occupants 709
# - same type of accommodation 670 same rental status 680
# - attended centre at same day (for this period) 53 date, 54 centre
# - 6139,6140 (gas cooking, heating type) same
# - same number vehicles 728 (risky because it says "by *you or* members
#    of your household")

# for comparison, do the same with rounded coordinates
# try to get to 58K
# see how many are in the living together data

# ==== match pairs using the mf_pairs data (unrounded coords)


fhs <- famhist %>% select(f.eid, f.20074.0.0, f.20075.0.0, f.6141.0.0, f.699.0.0,
                          f.709.0.0, f.670.0.0, f.680.0.0, f.53.0.0, f.54.0.0,
                          f.6139.0.0, f.6140.0.0, f.728.0.0, sex)
mf_pairs <- readr::read_csv(mf_pairs_file)
mf_pairs %<>% 
          left_join(fhs, by = c(ID.m = "f.eid")) %>% 
          left_join(fhs, by = c(ID.f = "f.eid"), suffix = c(".m", ".f"))

mf_pairs %<>% filter(
                f.670.0.0.f == f.670.0.0.m, # same kind of house
                f.680.0.0.f == f.680.0.0.m, # same rent/own status
                f.728.0.0.f == f.728.0.0.m, # same n vehicles
                f.699.0.0.f == f.699.0.0.m, # same length of time in hh
                f.709.0.0.f == f.709.0.0.m, # same n occupants of hh
                f.6141.0.0.f == 1,          # both living with spouse
                f.6141.0.0.m == 1
              )


# ==== validate using children ====

# find everyone who has a child in the data

parent_child <- relatedness %>% filter(relation == "parents") %>% select(ID1, ID2)

fh_age <- famhist %>% 
            select(f.eid, age_at_recruitment) %>% 
            rename(age = age_at_recruitment)
parent_child %<>% 
                left_join(fh_age, by = c(ID1 = "f.eid")) %>% 
                left_join(fh_age, by = c(ID2 = "f.eid"), suffix = c("1", "2")) %>% 
                mutate(
                  parent_id = ifelse(age1 > age2, ID1, ID2),
                  child_id  = ifelse(age1 > age2, ID2, ID1)
                ) %>% 
                filter(
                  abs(age1 - age2) > 10    # this also removes 924 "parent-children"...
                )

# 2 people with 3 "parents" (maybe siblings?)
parent_child %<>% 
                group_by(child_id) %>% 
                add_count() %>% 
                filter(n <= 2) %>% 
                ungroup()

parent_child %<>% select(parent_id, child_id)

# merge our dataset of children into our dataset of possible spouses

mf_w_parent <- mf_pairs %>% 
  left_join(parent_child, by = c("ID.m" = "parent_id")) %>% 
  left_join(parent_child, by = c("ID.f" = "parent_id"), suffix = c(".m", ".f")) %>% 
  filter(! is.na(child_id.m) | ! is.na(child_id.f))

# now, of the above:
# some have a child_id.x, a genetic child of x; 
# some have a child_id.y, a genetic child of y.
# we'd like to know how many have *the same* id as a child of both parents

# number of pairs where either have a child in the dataset
n_x_or_y_child <- nrow(mf_w_parent)
# number of pairs where both have the same child
n_xy_same <- mf_w_parent %>% filter(child_id.m == child_id.f) %>% nrow()
# number of pairs where both have a different child
# relying on filter() removing NAs here:
n_xy_different <- mf_w_parent %>% filter(child_id.m != child_id.f) %>% nrow()



# ==== match pairs using rounded coordinates ====

fhs <- famhist %>% select(f.eid, f.20074.0.0, f.20075.0.0, f.6141.0.0, f.699.0.0,
                          f.709.0.0, f.670.0.0, f.680.0.0, f.53.0.0, f.54.0.0,
                          f.6139.0.0, f.6140.0.0, f.728.0.0, sex)

mf_round <- inner_join(fhs, fhs, na_matches = "never",
                       by = c(
                         "f.20074.0.0",  # rounded easting 
                         "f.20075.0.0",  # rounded northing
                         "f.670.0.0",    # same kind of house
                         "f.680.0.0",    # same rent/own status
                         "f.728.0.0",    # same n vehicles
                         "f.699.0.0",    # same length of time in hh
                         "f.709.0.0",    # same n occups in hh
                         "f.54.0.0"      # same recruitment centre
                      ))
mf_round %<>% filter(
                f.eid.x != f.eid.y,   # not everyone who's the same person!
                f.eid.x < f.eid.y,    # delete one of each duplicated row
                sex.x != sex.y,       # no same-sex couples or "couples"
                f.6141.0.0.x == 1,    # both living with spouse 
                f.6141.0.0.y == 1     # both living with spouse 
              )
# they don't use f.6139/6140 re heating and cooking (multiple answers possible) 
# Nor f.53 (same date at recruitment centre,
# perhaps that would be a good one)

# remove multiple pairs
mf_round %<>% 
      group_by(f.eid.x) %>% 
      add_count() %>% 
      filter(n == 1) %>% 
      select(-n)
mf_round %<>% 
      group_by(f.eid.y) %>% 
      add_count() %>% 
      filter(n == 1) %>%
      select(-n)

# ==== validate against parent_child ====

mf_w_parent <- mf_round %>% 
  left_join(parent_child, by = c("f.eid.x" = "parent_id")) %>% 
  left_join(parent_child, by = c("f.eid.y" = "parent_id")) %>% 
  filter(! is.na(child_id.x) | ! is.na(child_id.y))

n_x_or_y_child <- nrow(mf_w_parent)
# number of pairs where both have the same child
n_xy_same <- mf_w_parent %>% filter(child_id.x == child_id.y) %>% nrow()

cat("Prop. genetic kids shared: ", n_xy_same/n_x_or_y_child)

# ==== compare to unrounded data ====

# match against the "unrounded coordinate living together" data
mf_round %<>% mutate(
                ID.m = ifelse(sex.x == 1, f.eid.x, f.eid.y),
                ID.f = ifelse(sex.x == 1, f.eid.y, f.eid.x)
              )
unrounded_coords <- readr::read_csv(
  "~/negative-selection-data/spouse_pair_info/UKB_out.mf_pairs_rebadged.csv")

mf_round_matched <- semi_join(mf_round, unrounded_coords,
                              by = c("ID.m", "ID.f"))

# about 23% of the sample weren't living in the same unrounded coordinate
nrow(mf_round)
nrow(mf_round_matched)

# a simple comparison; of all those pairs living in the same rounded coordinates,
# how many were living in the same unrounded coordinates?

mf_round_all <- inner_join(fhs, fhs, na_matches = "never",
                       by = c(
                         "f.20074.0.0",  # rounded easting 
                         "f.20075.0.0"
                       ))
mf_round_all %<>% filter(
                f.eid.x != f.eid.y,   # not everyone who's the same person!
                f.eid.x < f.eid.y     # delete one of each duplicated row
              )
mf_round_all %<>% mutate(
                ID.m = ifelse(sex.x == 1, f.eid.x, f.eid.y),
                ID.f = ifelse(sex.x == 1, f.eid.y, f.eid.x)
              )
mf_ra_matched <- semi_join(mf_round_all, unrounded_coords, 
                           by = c("ID.m", "ID.f"))

# Answer: less than 2% 
nrow(mf_ra_matched)
nrow(mf_round_all)
