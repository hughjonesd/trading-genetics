# Trading genetics in moba
# Fartein Ask Torvik, 2023

library(tidyverse)
library(vroom)
library(haven)

#-------------------------------------------------------------------------------
# loading
pop = vroom("N:/durable/data/registers/original/csv/w19_0634_faste_oppl_ut.csv")
edu = vroom("N:/durable/data/registers/original/csv/w19_0634_utd_1970_2018_ut.csv")
innt = vroom("N:/durable/data/registers/original/csv/w19_0634_innt_1993_2017_ut.csv")


#-------------------------------------------------------------------------------
# CREATE SAMPLE
#-------------------------------------------------------------------------------
# sample = moba parents
x <- vroom('N:/durable/data/moba/Original files/csv/PDB2601_SV_INFO_v12.csv')
x %>% count(is.na(M_ID_2601))
x %>% count(is.na(F_ID_2601))

# we can remove couples with missing informaiton because we need data from both spouses
x2 <- x %>% filter(!is.na(M_ID_2601)) %>% filter(!is.na(F_ID_2601)) 
x2 %>% count(duplicated(M_ID_2601))
x2 %>% count(duplicated(F_ID_2601))

# we need each man to be linked to one and only one woman and vice versa. 
# first, we can remove couples where both partners are duplicated
x2 %>% mutate(MF_ID = paste0(M_ID_2601,F_ID_2601)) %>% count(duplicated(MF_ID))
x3 <- x2 %>% mutate(MF_ID = paste0(M_ID_2601,F_ID_2601)) %>% filter(!duplicated(MF_ID))
x3
x3 %>% count(duplicated(M_ID_2601)) # 154
x3 %>% count(duplicated(F_ID_2601)) # 120

# there are 154 women with children with muliple moba-men and 120 men with children with multiple moba-men
# remove first duplicated women, then duplicted men. 
# this could be done in different ways, but the low number of cases makes this decision unimportant.
set.seed(1)
x4 <- x3 %>% 
  arrange(sample(nrow(.))) %>% 
  filter(!duplicated(M_ID_2601)) %>% 
  filter(!duplicated(F_ID_2601)) %>% 
  arrange(PREG_ID_2601) %>% 
  select(-MF_ID)
x4

# make a long version of the data, individual-based
x5 <- x4 %>% 
#  slice(1:10) %>% 
  rename(PREGID = PREG_ID_2601, MOBAID_M=M_ID_2601, MOBAID_F=F_ID_2601) %>% 
  mutate(male_M = 0) %>%  # m for mother
  mutate(male_F = 1) %>%  # f for father
  identity() %>% 
  pivot_longer(
    cols = c(-PREGID, -FAAR), 
    names_to = c('.value','role'),
    names_pattern = '(.*)_(.)'
  )

#-------------------------------------------------------------------------------
# ADD BASIC VARIABLES FROM STATISTICS NORWAY
#-------------------------------------------------------------------------------
# add SSB id number
mobassblink <- vroom('N:/durable/data/moba/linkage/PDB2601_kobling_SSB_v12.csv') %>% 
  filter(!rolle == 'SU2PT_CHILD') %>% 
  select(-barn_nr) %>% 
  mutate(
    role = case_when(
      rolle == "SU2PT_MOTHER" ~ "F",
      rolle == "SU2PT_FATHER" ~ "M",
      TRUE ~ rolle  # Keep the original value if it's not "SU2PT_MOTHER" or "SU2PT_FATHER"
    )
  ) %>% 
  select(-rolle) %>% 
  rename(PREGID = PREG_ID_2601)
mobassblink
x6 <- x5 %>% left_join(mobassblink)

# permanent information

# year and month of birth
pop2 <- pop %>% 
  mutate(
    byear = as.integer(substr(foedsels_aar_mnd, 1, 4)),
    bmonth = as.integer(substr(foedsels_aar_mnd, 5, 6))
  )
pop2

# calculate parity (by mother)
pop3 <- pop2 %>% filter(!is.na(mor_lnr)) # remove missing mothers
pop4 <- pop3 %>% filter((foedsels_aar_mnd  != doeds_aar_mnd) %in% c(FALSE, NA)) # remove stillbirths
pop5 <- pop4 %>% arrange(mor_lnr) # help while testing, not needed
pop6 <- pop5 %>% 
#  slice(1:1000) %>% 
  group_by(mor_lnr) %>% 
  arrange(foedsels_aar_mnd) %>% 
  mutate(parity = row_number()) %>% 
  add_tally(name='famsize') %>% 
  ungroup()
count(pop6, parity)
count(pop6, famsize)
pop6

# parental age at birth
minipop <- pop2 %>% select(w19_0634_lnr, byear, bmonth)
motherpop <- minipop %>% rename(byearmother = byear, bmonthmother = bmonth)
fatherpop <- minipop %>% rename(byearfather = byear, bmonthfather = bmonth)
pop7 <- pop6 %>% left_join(motherpop, by=c('mor_lnr'='w19_0634_lnr'))
pop7
pop8 <- pop7 %>% left_join(fatherpop, by=c('far_lnr'='w19_0634_lnr'))
pop8
pop9 <- pop8 %>% 
  mutate(fatherage = byear - byearfather) %>% 
  mutate(fatherageexact = byear+(bmonth-1)/12 - byearfather - (bmonthfather-1)/12) %>% 
  mutate(motherage = byear - byearmother) %>% 
  mutate(motherageexact = byear+(bmonth-1)/12 - byearmother - (bmonthmother-1)/12) 
count(pop9, fatherage)  %>% print(n=100)
count(pop9, motherage)  %>% print(n=100)
cor(pop9$fatherage, pop9$fatherageexact, use='c')  
cor(pop9$motherage, pop9$motherageexact, use='c')  
cor(pop9$motherage, pop9$fatherage, use='c')  

# link basic SSB info to MOBA
poplink <- pop9 %>% 
  rename(sex = kjoenn) %>% 
  select(w19_0634_lnr, sex, parity, famsize, byear, bmonth, 
         fatherage,fatherageexact,motherage,motherageexact)

x7 <- x6 %>% left_join(poplink)

x7 %>% count(is.na(motherage))
x7 %>% count(is.na(fatherage))

#-------------------------------------------------------------------------------
# ADD EDUCATIONAL ATTAINMENT FROM SSB
#-------------------------------------------------------------------------------

#university edu
x7 %>% count(byear) %>% print(n=100)
edu2 <- edu %>% filter(w19_0634_lnr %in% x7$w19_0634_lnr)
edu3 <- edu2 %>% select(w19_0634_lnr, all_of(paste0('BU_',1980:2018)))

edu4 <- edu3 %>% 
  pivot_longer(
    cols=-w19_0634_lnr,
    names_to=c('.value','year'),
    names_sep='_'
  )
edu5 <- edu4 %>% filter(!is.na(BU))
edu5
minipop <- pop2 %>% select(w19_0634_lnr, byear)
edu6 <- edu5 %>% left_join(minipop)
edu6
edu6$year = as.numeric(edu6$year)
edu6
edu7 <- edu6 %>% mutate(eduyear = year-byear)
edu7
# chose closest valid to 30 (30.1 for preference for 31 over 29)
edu8 <- edu7 %>% 
  mutate(diff_from_30 = abs(eduyear-30.1)) %>% 
  group_by(w19_0634_lnr) %>% 
  arrange(diff_from_30) %>% 
  slice(1) %>% 
  ungroup()
edu8 %>% count(eduyear) %>% print(n=50)

edu9 <- edu8 %>% 
  mutate(edulevel = as.integer(substr(BU, 1, 1)))

eduA <- edu9 %>% select(w19_0634_lnr,edulevel,eduyear)
eduA$edulevel[eduA$edulevel == 9] = NA
eduA %>% count(edulevel)

x8 <- x7 %>% left_join(eduA)
x8

#-------------------------------------------------------------------------------
# ADD INCOME FROM SSB
#-------------------------------------------------------------------------------

#income
innt2 <- innt %>% select(w19_0634_lnr, aargang , wies)
minipop <- pop2 %>% select(w19_0634_lnr, byear, kjoenn)
innt3 <- innt2 %>% left_join(minipop)

# make percentile rank in year, by sex and age (e.g. among 30 year old men in 1999)
innt4 <- innt3 %>% group_by(aargang, kjoenn, byear) %>% 
  mutate(wies = ifelse(wies < 0,0,wies)) %>% 
  mutate(incomedec = ntile(wies, 10)) %>% 
  mutate(incomez = as.numeric(scale(wies))) %>% 
  ungroup()
innt4
# filer moba participants
innt5 <- innt4 %>% filter(w19_0634_lnr %in% x8$w19_0634_lnr)
innt5
innt6 <- innt5 %>% mutate(incomeyear = aargang-byear)
innt6
# chose closest valid to 30 (30.1 for preference for 31 over 29)
# YES, THIS MAY HAVE SOME ISSUES FOR INCOME/AGE 30
# THER APPROACHES HAVE ISSUES TOO
# DISCUSS LATER
innt7 <- innt6 %>% 
  mutate(diff_from_25 = abs(incomeyear-25.1)) %>% 
  group_by(w19_0634_lnr) %>% 
  arrange(diff_from_25) %>% 
  slice(1) %>% 
  ungroup()
innt7
innt8 <- innt7 %>% select(w19_0634_lnr,incomedec,incomez,incomeyear)
innt8 %>% count(incomeyear) %>% print(n=100)
innt8$income25dec_strict = innt8$incomedec
innt8$income25dec_strict[innt8$incomeyear!=25] = NA
innt8$income25z_strict = innt8$incomez
innt8$income25z_strict[innt8$incomeyear!=25] = NA




innt7b <- innt6 %>% 
  mutate(diff_from_35 = abs(incomeyear-35.1)) %>% 
  group_by(w19_0634_lnr) %>% 
  arrange(diff_from_35) %>% 
  slice(1) %>% 
  ungroup()
innt7b

innt8b <- innt7b %>% select(w19_0634_lnr,incomedec,incomez,incomeyear) %>% 
  rename(income35dec=incomedec,income35z=incomez,incomeyear35=incomeyear)
innt8b %>% count(incomeyear35) %>% print(n=100)
innt8b$income35dec_strict = innt8b$income35dec
innt8b$income35dec_strict[innt8b$incomeyear35!=35] = NA
innt8b$income35z_strict = innt8b$income35z
innt8b$income35z_strict[innt8b$incomeyear35!=35] = NA




innt7c <- innt6 %>% 
  mutate(diff_from_30 = abs(incomeyear-30.1)) %>% 
  group_by(w19_0634_lnr) %>% 
  arrange(diff_from_30) %>% 
  slice(1) %>% 
  ungroup()
innt7c

innt8c <- innt7c %>% select(w19_0634_lnr,incomedec,incomez,incomeyear) %>% 
  rename(income30dec=incomedec,income30z=incomez,incomeyear30=incomeyear)
innt8c %>% count(incomeyear30) %>% print(n=100)
innt8c$income30dec_strict = innt8c$income30dec
innt8c$income30dec_strict[innt8c$incomeyear30!=30] = NA
innt8c$income30z_strict = innt8c$income30z
innt8c$income30z_strict[innt8c$incomeyear30!=30] = NA


x9 <- x8 %>% left_join(innt8) %>% left_join(innt8b) %>% left_join(innt8c)
x9
x9 %>% count(incomedec) %>% print(n=101)




#-------------------------------------------------------------------------------
# ADD PHENOTYPES FROM MOBA
#-------------------------------------------------------------------------------
# height
# bmi
mothers = vroom('N:/durable/data/moba/Original files/csv/PDB2601_Q1_v12.csv')
mothers2 <- mothers %>% select(PREG_ID_2601, AA85, AA87, AA88, AA89) %>% 
  rename(mother_weight = AA85) %>% 
  rename(mother_height = AA87) %>% 
  rename(father_weight_mreport = AA89) %>% 
  rename(father_height_mreport = AA88) %>% 
  mutate(mother_bmi = mother_weight/(mother_height/100)^2) %>% 
  mutate(father_bmi_mreport = father_weight_mreport/(father_height_mreport/100)^2)
mothers2
# recode implausible values
mothers3 = mothers2
mothers3$mother_height[mothers3$mother_height<50]=NA
mothers3$mother_height[which(mothers3$mother_height<100)]=100+mothers3$mother_height[which(mothers3$mother_height<100)]
mothers3 %>% count(mother_height) %>% print(n=100)
mothers3$father_height_mreport[mothers3$father_height_mreport<50]=NA
mothers3$father_height_mreport[which(mothers3$father_height_mreport<100)]=100+mothers3$father_height_mreport[which(mothers3$father_height_mreport<100)]
mothers3 %>% count(father_height_mreport) %>% print(n=100)
mothers3$mother_bmi[mothers3$mother_bmi<14]=NA
mothers3$mother_bmi[mothers3$mother_bmi>50]=NA
mothers3 %>% count(round(mother_bmi) ) %>% print(n=200)
mothers3$father_bmi_mreport[mothers3$father_bmi_mreport<14]=NA
mothers3$father_bmi_mreport[mothers3$father_bmi_mreport>50]=NA
mothers3 %>% count(round(father_bmi_mreport) ) %>% print(n=200)
mothers_pheno <- mothers3 %>% 
  select(PREG_ID_2601, mother_height, mother_bmi) %>% 
  mutate(role='M')
fathers_pheno_mreport <- mothers3 %>% 
  select(PREG_ID_2601, father_height_mreport, father_bmi_mreport) %>% 
  mutate(role='FM') %>% 
  rename(father_height=father_height_mreport, father_bmi=father_bmi_mreport)%>% 
  filter(!is.na(father_height))
fathers = vroom('N:/durable/data/moba/Original files/csv/PDB2601_QF_v12.csv')
fathers2 <- fathers %>% select(PREG_ID_2601, FF333, FF334) %>% 
  rename(father_weight = FF334) %>% 
  rename(father_height = FF333) %>% 
  mutate(father_bmi = father_weight/(father_height/100)^2) %>% 
  mutate(role='F') 
fathers2
fathers2$father_height[fathers2$father_height<50]=NA
fathers2$father_height[which(fathers2$father_height<100)]=100+fathers2$father_height[which(fathers2$father_height<100)]
fathers2 %>% count(father_height) %>% print(n=100)
fathers2$father_bmi[fathers2$father_bmi<14]=NA
fathers2$father_bmi[fathers2$father_bmi>50]=NA
fathers2 %>% count(round(father_bmi) ) %>% print(n=200)
fathers2 <- fathers2 %>% select(-father_weight)
fathers2

# join fathers, keep own report if available
fathers_pheno <- fathers2 %>% filter(!is.na(father_height))
fathers_pheno_mreport %>% count(is.na(father_height ))
fathers_pheno_mf = bind_rows(fathers_pheno, fathers_pheno_mreport)
fathers_pheno_mf$priority = 0
fathers_pheno_mf$priority[fathers_pheno_mf$role=='MF'] = 1
fathers_pheno_mf = fathers_pheno_mf %>% group_by(PREG_ID_2601) %>% arrange(priority) %>% slice(1) %>% ungroup()
table(fathers_pheno_mf$role)
fathers_pheno_mf$role[fathers_pheno_mf$role=='FM'] = 'F'
table(fathers_pheno_mf$role)
table(duplicated(fathers_pheno_mf$PREG_ID_2601))

# join mothers and fathers
fa_link = fathers_pheno_mf %>% select(-priority) %>% 
  rename(height=father_height, bmi=father_bmi )
fa_link
mo_link=mothers_pheno%>% 
  rename(height=mother_height, bmi=mother_bmi )
mo_link
heightbmi <- bind_rows(fa_link, mo_link) %>% 
  rename(PREGID=PREG_ID_2601)
heightbmi %>% count(role)

#link to dataset
xA <- x9 %>% left_join(heightbmi)
xA

#write_csv(xA, 'xA-tempsave.csv')

#-------------------------------------------------------------------------------
# ADD GENETIC DATA
#-------------------------------------------------------------------------------

# linked IDs from Hans Fredrik (check script)
IDs_AUG22_long <- vroom("N:/durable/projects/HansFredrik/AM/misc_files/IDs_AUG22_long.csv")
# clean up
ids2 <- IDs_AUG22_long %>% 
  filter(Role != 'child') %>% 
  select(-barn_nr, -PREGBARN) %>% 
  rename(role=Role) %>% 
  mutate(role = ifelse(role == "mother", "F", "M"))
ids2

source('hansfredrik/PreparingPGSs_new_pipeline_fato.R')
PGS_use

ids3 <- ids2 %>% left_join(PGS_use)
ids3
ids4 <- ids3 %>% 
  rename(PREGID=PREG_ID_2601) %>% 
  rename(w19_0634_lnr_copy =w19_0634_lnr ) 
ids4  

# merge
xB = left_join(xA, ids4)
table(xB$w19_0634_lnr == xB$w19_0634_lnr_copy)

# 
cor (xB$eapgsresid, xB$edulevel, use='c')
cor (xB$heightpgsresid, xB$height, use='c')




write_csv(xB, path='trading_dataset_20230705.csv')

#-------------------------------------------------------------------------------

# variables currently included:
# PREGID: implicitly spouse ID, use pivot_wider
# FAAR: birth year of moba child (3rd gen)
# role: role for index person: M(opther) or F(ather) (2nd gen; index generation)
# MOBAID: id for index person
# male: sex, coded from "role"
# w19_0634_lnr: SSB id
# sex: sex from SSB (identical to "male", but reveals that 17 fathers are co-mothers)
# parity: parity by mother of index person
# famsize: number of children born to mother index person
# byear: birth year
# bmonth: birth month
# fatherage: father's age when index person born
# fatherageexact: father's age when index person born, including month
# motherage: mother's age when index person born
# motherageexact: mother's age when index person born, including month
# edulevel: NUS edu level at age 30
# eduyear: age at education, mostly 30, but other age if missing
# income: income percentil at age 30 (or closest valid age) among same sex, age, and birth year
# incomeyear: age from which income at age 30 chosen
# height
# bmi
# eapgs
# eapgsresid
# depressionpgs
# depressionpgsresid
# bmipgs
# bmipgsresid
# heightpgs
# heightpgsresid
# smokingpgs
# smkoingpgsresid

# note: iq and overall self-rated health not available
