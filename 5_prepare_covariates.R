library(tidyverse)
library(zoo)
library(purrr)
library(lavaan)


source('helper_functions.R')

# main dataset
data_all <- readRDS('vars_pseudoscales.Rds')


# define time period used
# c(2015) for main analysis or c(2004, 2005, 2006) for sensitivity analysis
relevant_years <- c(2015) 
suff <- '2015'












    #################################################
  ##                                                 ##
 ###                                                 ###
####                 helper variables                ####
 ###                                                 ###
  ##                                                 ##
    #################################################

meds <- read.csv('output_files/meds_de-branded.csv', sep = '|')
meds$date <- as.Date(meds$date, format = '%d/%m/%Y')
meds$year <- as.numeric(format(as.Date(meds$date), "%Y"))

# remove invalid application route
meds <- filter(meds, application == 0 & year %in% relevant_years)

scales <- read.csv('aas_combined.csv')
# only retain those that appear in the sample
scales <- filter(scales, drug %in% unique(meds$name_generic))
# add number of scales that include a drug and export
scales <- scales %>%
  mutate(n = rowSums(select(., -drug) > 0, na.rm = TRUE))
write.csv(scales, 'output_files/aas_in_sample.csv', row.names = FALSE)

scales$n <- NULL

# to numeric
scales[, 2:ncol(scales)] <- sapply(scales[, 2:ncol(scales)], as.character)
# transpose so that scales are rows
scales <- data.frame(t(scales))
# remove row with drug names
scales <- scales[2:nrow(scales), ]
# create data frame with counts of different scores for each scale
scales$scale <- rownames(scales)
rownames(scales) <- seq(1,23)

scales <- scales %>%
  rowwise() %>%
  mutate(
    count_4 = sum(c_across(X1:X242) == 4),
    count_3 = sum(c_across(X1:X242) == 3),
    count_2 = sum(c_across(X1:X242) == 2),
    count_1 = sum(c_across(X1:X242) == 1),
    count_05 = sum(c_across(X1:X242) == 0.5),
    count_any = sum(c_across(X1:X242) > 0)
  ) %>%
  ungroup()

scales <- as.data.frame(subset(scales, select = c(scale, count_4, count_3, 
                                                  count_2, count_1, count_05, 
                                                  count_any)))
write.csv(scales, 'output_files/scales_summary.csv', row.names = FALSE)




invalid_dates <- as.Date(c('01/01/1900', '01/01/1901', '02/02/1902', 
                           '03/03/1903', '07/07/2037'), format = '%d/%m/%Y')

## This code combines data fields 41270, 41271, 41280, and 41281, 
# which refer to inpatient diagnoses and dates of admission.
diagnoses <- data_all %>%
  select(eid, starts_with(c('X41270.', 'X41280.', 'X41271.', 'X41281.'))) %>%
  rename(id = eid)

diagnoses[diagnoses=='']  <- NA 
diagnoses <- as.data.frame(sapply(diagnoses, as.character))

# separate sources of diagnoses (ICD9 vs 10) and dates vs. diagnosis codes
icd9 <- diagnoses %>% select(c('id', starts_with('X41271')))
icd9_date <- diagnoses %>% select(c('id', starts_with('X41281')))
icd10 <- diagnoses %>% select(c('id', starts_with('X41270')))
icd10_date <- diagnoses %>% select(c('id', starts_with('X41280')))

# keep only rows without NAs
icd9 <- icd9[rowSums(is.na(icd9))!=ncol(icd9)-1,]
icd9_date <- icd9_date[rowSums(is.na(icd9_date))!=ncol(icd9_date)-1,]
icd10 <- icd10[rowSums(is.na(icd10))!=ncol(icd10)-1,]
icd10_date <- icd10_date[rowSums(is.na(icd10_date))!=ncol(icd10_date)-1,]

# transform to long-type format
icd9_long <- icd9 %>%  pivot_longer(-id, names_to = 'diagnosis', 
                                    values_drop_na=TRUE)
colnames(icd9_long) <- c('id', 'column', 'diagnosis')
icd9_long$column <- sub('X41271.', '', icd9_long$column)

icd9_date_long <- icd9_date %>%  pivot_longer(-id, names_to = 'diagnosis', 
                                              values_drop_na=TRUE)
colnames(icd9_date_long) <- c('id', 'column', 'date')
icd9_date_long$column <- sub('X41281.', '', icd9_date_long$column)

icd10_long <- icd10 %>%  pivot_longer(-id, names_to = 'diagnosis', 
                                      values_drop_na=TRUE)
colnames(icd10_long) <- c('id', 'column', 'diagnosis')
icd10_long$column <- sub('X41270.', '', icd10_long$column)

icd10_date_long <- icd10_date %>%  pivot_longer(-id, names_to = 'diagnosis', 
                                                values_drop_na=TRUE)
colnames(icd10_date_long) <- c('id', 'column', 'date')
icd10_date_long$column <- sub('X41280.', '', icd10_date_long$column)

# combine all diagnoses
icd9 <- merge(icd9_long, icd9_date_long, by = c('id', 'column'))
icd9$column <- NULL; icd9$version <- 'icd9'
icd10 <- merge(icd10_long, icd10_date_long, by = c('id', 'column'))
icd10$column <- NULL; icd10$version <- 'icd10'

inpatient <- rbind(icd9, icd10)
saveRDS(inpatient, 'output_files/inpatient.Rds')
rm(icd10, icd10_date, icd10_date_long, icd10_long, 
   icd9, icd9_date, icd9_date_long, icd9_long, diagnoses)
gc()












    #################################################
  ##                                                 ##
 ###                                                 ###
####        demographic and lifestyle factors        ####
 ###                                                 ###
  ##                                                 ##
    #################################################



## assessment dates, age at assessments and sex
dems <- data_all %>% 
  select(c(eid, starts_with(c('X31.', 'X34.', 'X52.', 'X48.')))) %>%
  rename(sex = X31.0.0, birth_year = X34.0.0, birth_month = X52.0.0,
         waist_0 = 'X48.0.0') %>%
  select(-starts_with('X'))
dems$birth_year <- as.character(dems$birth_year); dems$birth_month <- 
  as.character(dems$birth_month)
dems$birth_date <- as.Date(paste0('01/', dems$birth_month, '/' ,dems$birth_year), 
                           format = '%d/%m/%Y')
dems <- rename(dems, id = eid)
write.csv(select(dems, -waist_0), 'output_files/age_sex_formatted.csv')




## smoking
smoking <- data_all %>%
  select(eid, X20116.0.0)
colnames(smoking) <- c('id', 'smoking_0')
smoking[smoking == -3] <- NA




## alcohol
alcohol <- data_all %>%
  select(eid, X1558.0.0) %>%
  rename(id = eid, alc_freq_0 = X1558.0.0)
alcohol[alcohol == -3] <- NA





## level of physical activity in the past 4 weeks:
phys_act <- data_all %>%
  select(eid, X6164.0.0) %>%
  mutate(phys_act_0 = apply(select(., X6164.0.0), 1, 
                            function(x) phys_act_classify(x))) %>%
  rename(id = eid) %>%
  select(id, phys_act_0)









    #################################################
  ##                                                 ##
 ###                                                 ###
####             social and environmental            ####
 ###                                                 ###
  ##                                                 ##
    #################################################

## education
# create subsets of the education frame for each assessment; change -7 to NA, 
# and set to 1 if graduate degree present
education <- data_all %>% 
  select(eid, X6138.0.0) %>% 
  rename(id = eid)
education[education == -3] <- NA
education <- education %>%
  filter(rowSums(is.na(select(., -id))) != ncol(.) - 1) %>%
  mutate(education_0 = as.integer(rowSums(select(., X6138.0.0) == 1, na.rm = TRUE) > 0)) %>%
  select(id, education_0)




## deprivation
deprivation <- data_all %>% 
  select(c(eid, starts_with(c('X189')))) %>%
  rename(id = eid, deprivation = X189.0.0)




## air pollution
# create means across years and rename columns
pollution <- data_all %>%
  select(eid, starts_with(c('X24003.', 'X24004.', 'X24006.', 'X24016.', 
                            'X24017.', 'X24018.'))) %>%
  rename(id = eid, 
         pollution_nox = X24004.0.0, 
         pollution_25 = X24006.0.0) %>%
  mutate(pollution_no2 = rowMeans(across(c(
    X24016.0.0, X24017.0.0, X24018.0.0, X24003.0.0)), na.rm = TRUE)) %>%
  select(id, pollution_nox, pollution_no2, pollution_25)
# scale the variables
pollution[, c('pollution_nox', 'pollution_no2', 'pollution_25')] <- 
  lapply(pollution[, c('pollution_nox', 'pollution_no2', 'pollution_25')], 
         function(x) as.vector(scale(x)))
# principal component of nox, no2, and PM2.5, since they were all associated with 
# dementia in a recent systematic review (Peters et al., 2019)
pollution$pollution_pc <- psych::principal(
  r = select(pollution, pollution_no2, pollution_nox, pollution_25), 
  nfactors=1, 
  rotate='none', 
  scores=T, 
  covar=FALSE, 
  missing=F)$scores[, 1]
pollution <- pollution %>%
  select(id, pollution_pc)




## social isolation, loneliness, depressed mood
social <- data_all %>% 
  select(eid, X709.0.0, X1031.0.0, X6160.0.0, X2020.0.0, X2050.0.0) %>%
  rename(id = eid, lonely_0 = X2020.0.0, depressed_0 = X2050.0.0)
social$X709.0.0[social$X709.0.0 == -3 | social$X709.0.0 == -1] <- NA
social$X709.0.0[social$X709.0.0 > 1] <- 0
social$X1031.0.0[social$X1031.0.0 == -3 | social$X1031.0.0 == -1] <- NA
social$X1031.0.0[social$X1031.0.0 <= 4] <- 0
social$X1031.0.0[social$X1031.0.0 > 4] <- 1
social$X6160.0.0[social$X6160.0.0 == -3] <- NA
social$X6160.0.0[social$X6160.0.0 > 0] <- 0
social$X6160.0.0[social$X6160.0.0 == -7] <- 1
social <- social %>% 
  mutate(soc_isol_0 = rowSums(across(c(X709.0.0 , X1031.0.0, X6160.0.0)), na.rm = TRUE))
social$soc_isol_0[social$soc_isol_0 >= 2] <- 1
# among those that have no social isolation, check for NAs in any of the three questions
social <- social %>% 
  mutate(na_count_0 = rowSums(across(c(X709.0.0, X1031.0.0, X6160.0.0), is.na)))
# those with a score of 0 or 1 could have such a low score because of NAs, so we have to make sure
# if only one is NA and the other two are 0, keep categorised as not isolated
# if two or more are missing, categorise as NA
social$soc_isol_0[social$na_count_0 >= 2] <- NA
# loneliness
social$lonely_0[social$lonely_0 == -3 | social$lonely_0 == -1] <- NA
# depressed mood
social$depressed_0[social$depressed_0 == -3 | social$depressed_0 == -1] <- NA
social <- social %>% 
  select(id, lonely_0, depressed_0, soc_isol_0)




## cognition
cognition <- data_all %>% 
  select(eid, X20016.0.0, X20018.0.0, X20023.0.0, X399.0.2, X4282.0.0)
colnames(cognition) <- c('id', 'VNR_0', 'ProsMem_0', 'RT_0','VisMem_0', 'NM_0')
# remove outliers 
cognition[, c('VNR_0', 'ProsMem_0', 'RT_0','VisMem_0', 'NM_0')] <- 
  lapply(cognition[, c('VNR_0', 'ProsMem_0', 'RT_0','VisMem_0', 'NM_0')], 
         outliers, var_metric=4, method='SD')
# calculate the latent G
cognition$RT_0 = log(cognition$RT_0)
cognition$VisMem_0 = log(cognition$VisMem_0+1)
cognition$ProsMem_0[which(cognition$ProsMem_0 != 1)] = 0
cognition$ProsMem_0[which(cognition$ProsMem_0 == 1)] = 1
model_0 <- '
            # Structural relation
            g =~ VNR_0 + RT_0 + VisMem_0 + ProsMem_0 + NM_0
            
'
fit_0 <- sem(model_0, data=cognition, missing='fiml.x')
# extract the g-values
cognition$g_0 <- as.vector(predict(fit_0, cognition)) 
cognition <- cognition %>%
  select(id, g_0)
rm(fit_0)








    #################################################
  ##                                                 ##
 ###                                                 ###
####                relevant disorders               ####
 ###                                                 ###
  ##                                                 ##
    #################################################


## dementia
dementia <- data_all %>% select(c(eid, starts_with(c('X42018.'))))
colnames(dementia) <- c('id', 'dementia_date')
dementia$dementia_date <- as.Date(dementia$dementia_date, format <- '%Y-%m-%d')
dementia$dementia <- 0; dementia$dementia[!is.na(dementia$dementia_date)] <- 1
dementia$dementia[dementia$dementia_date == as.Date('1900-01-01', format = '%Y-%m-%d')] <- NA
dementia$dementia_date[dementia$dementia_date == as.Date('1900-01-01', format = '%Y-%m-%d')] <- NA




## death
death <- data_all %>% 
  select(c(eid, X40000.0.0))
colnames(death) <- c('id', 'death_date')
death$death_date <- as.Date(death$death_date, format = '%Y-%m-%d')
death$death <- 0; death$death[!is.na(death$death_date)] <- 1




## mood disorders
mood_ado <- data_all %>% 
  select(c(eid, starts_with(c('X130890.', 'X130892.', 'X130894.', 'X130896.', 
                              'X130898.', 'X130900.', 'X130902.')))) %>%
  mutate(across(starts_with('X'), ~ as.Date(., format = '%Y-%m-%d'))) %>%
  mutate(mood_dis_date = reduce(across(starts_with('X')), pmin, na.rm = TRUE)) %>%
  rename(id = eid) %>%
  select(id, mood_dis_date)
mood_ado$mood_dis <- 0; mood_ado$mood_dis[!is.na(mood_ado$mood_dis_date)] <- 1
mood_ado$mood_dis[mood_ado$mood_dis_date == as.Date('1902-02-02', format = '%Y-%m-%d')] <- NA
mood_ado$mood_dis_date[mood_ado$mood_dis_date == as.Date('1902-02-02', format = '%Y-%m-%d')] <- NA




## psychotic disorders
psych_ado <- data_all %>%
  select(c(eid, starts_with(c('X130874.', 'X130876.', 'X130878.', 'X130880.', 
                              'X130882.', 'X130884.', 'X130886.', 'X130888.')))) %>%
  mutate(across(starts_with('X'), ~ as.Date(., format = '%Y-%m-%d'))) %>%
  mutate(psych_dis_date = reduce(across(starts_with('X')), pmin, na.rm = TRUE)) %>%
  rename(id = eid) %>%
  select(id, psych_dis_date)
psych_ado$psych_dis <- 0; psych_ado$psych_dis[!is.na(psych_ado$psych_dis_date)] <- 1




## diabetes
diabetes <- data_all %>%
  select(eid, starts_with(c('X130706.', 'X130708.', 'X130710.', 'X130712.', 'X130714.'))) %>%
  mutate(across(starts_with('X'), ~as.Date(., format = '%Y-%m-%d'))) %>%
  mutate(diabetes_date = reduce(across(starts_with('X')), pmin, na.rm = TRUE)) %>%
  rename(id = eid) %>%
  select(id, diabetes_date)
diabetes$diabetes <- 0; diabetes$diabetes[!is.na(diabetes$diabetes_date)] <- 1
diabetes$diabetes[diabetes$diabetes_date == as.Date('1902-02-02', format = '%Y-%m-%d')] <- NA
diabetes$diabetes_date[diabetes$diabetes_date == as.Date('1902-02-02', format = '%Y-%m-%d')] <- NA




## hypertension
hypertension <- data_all %>%
  select(eid, X131286.0.0) %>%
  rename(id = eid, hypertension_date = X131286.0.0)
hypertension$hypertension_date <- as.Date(hypertension$hypertension_date, format = '%Y-%m-%d')
hypertension$hypertension <- 0
hypertension$hypertension[!is.na(hypertension$hypertension_date)] <- 1
hypertension$hypertension[hypertension$hypertension_date == 
                            as.Date('1902-02-02', format = '%Y-%m-%d')] <- NA
hypertension$hypertension_date[hypertension$hypertension_date == 
                                 as.Date('1902-02-02', format = '%Y-%m-%d')] <- NA




## hyperlipeademia
hyperchol <- data_all %>%
  select(eid, X130814.0.0) %>% 
  rename(id = eid, hyperlip_date = X130814.0.0)
hyperchol$hyperlip_date <- as.Date(hyperchol$hyperlip_date, format = '%Y-%m-%d')
hyperchol$hyperlip <- 0; hyperchol$hyperlip[!is.na(hyperchol$hyperlip_date)] <- 1
hyperchol$hyperlip[hyperchol$hyperlip_date == 
                     as.Date('1902-02-02', format = '%Y-%m-%d')] <- NA
hyperchol$hyperlip_date[hyperchol$hyperlip_date == 
                          as.Date('1902-02-02', format = '%Y-%m-%d')] <- NA




## visual impairment, sleep disorders
sleep_vision <- data_all %>%
  select(eid, starts_with(c('X131212.', 'X131060.', 'X130920.'))) %>%
  filter(rowSums(across(starts_with('X')) == '') < 3) %>%
  mutate(across(starts_with('X'), ~as.Date(., format = '%Y-%m-%d'))) %>%
  mutate(sleep_dis_any_date = reduce(across(c('X131060.0.0', 'X130920.0.0')), pmin, na.rm = TRUE)) %>%
  rename(id = eid, vision_problem_date = X131212.0.0) %>%
  select(-starts_with('X'))
sleep_vision$vision_problem <- 0
sleep_vision$vision_problem[!is.na(sleep_vision$vision_problem_date)] <- 1
sleep_vision$vision_problem[sleep_vision$vision_problem_date == 
                              as.Date('1902-02-02', format = '%Y-%m-%d')] <- 99
sleep_vision$vision_problem_date[sleep_vision$vision_problem_date == 
                                   as.Date('1902-02-02', format = '%Y-%m-%d')] <- NA
sleep_vision$sleep_dis_any <- 0
sleep_vision$sleep_dis_any[!is.na(sleep_vision$sleep_dis_any_date)] <- 1




## CNS injury/disorders except dementia
cns_dis <- data_all %>%
  select(eid, starts_with(c('X42028.', 'X42006.',  'X131000.', 'X131114.', 
                            'X130992.', 'X130994.', 'X130996.', 'X130998.', 
                            'X131002.', 'X131004.', 'X131006.', 'X131008.', 
                            'X131010.', 'X131012.', 'X131014.', 'X131016.', 
                            'X131018.', 'X131020.', 'X131022.', 'X131024.', 
                            'X131026.', 'X131028.', 'X131030.', 'X131032', 
                            'X131038.', 'X131040.', 'X131042.', 'X131044.', 
                            'X131046.', 'X131048.', 'X131050.', 'X131056.', 
                            'X131058.', 'X131100.', 'X131110.', 'X131112.', 
                            'X131116.', 'X131120.', 'X131370.', 'X131372.', 
                            'X131374.', 'X131376.', 'X131378.'))) %>%
  mutate(across(starts_with('X'), ~as.Date(., format = '%Y-%m-%d'))) %>%
  mutate(cns_vasc_date = reduce(across(c(X42006.0.0, X131370.0.0, X131372.0.0, 
                                         X131374.0.0, X131376.0.0, X131378.0.0, )), 
                                pmin, na.rm = TRUE)) %>%
  mutate(cns_infl_date = reduce(across(c(X131000.0.0, X130992.0.0, X130994.0.0, 
                                         X130996.0.0, X130998.0.0, X131002.0.0,
                                         X131004.0.0, X131006.0.0, X131008.0.0, 
                                         X131010.0.0)), pmin, na.rm = TRUE)) %>%
  mutate(cns_atroph_date = reduce(across(c(X131012.0.0, X131014.0.0, X131016.0.0, 
                                           X131018.0.0, X131020.0.0,
                                           X42028.0.0)), 
                                  pmin, na.rm = TRUE)) %>%
  mutate(cns_mov_date = reduce(across(c(X131022.0.0, X131024.0.0, X131026.0.0, 
                                        X131028.0.0, X131030.0.0,
                                        X131032.0.0)), pmin, na.rm = TRUE)) %>%
  mutate(cns_demyel_date = reduce(across(c(X131042.0.0, X131044.0.0, X131046.0.0)), 
                                  pmin, na.rm = TRUE)) %>%
  mutate(cns_parox_date = reduce(across(c(X131048.0.0, X131050.0.0, X131056.0.0, 
                                          X131058.0.0)), pmin, na.rm = TRUE)) %>%
  mutate(cns_other_date = reduce(across(c(X131038.0.0, X131040.0.0, X131100.0.0, 
                                          X131110.0.0, X131112.0.0, X131114.0.0, 
                                          X131116.0.0, X131120.0.0)), 
                                 pmin, na.rm = TRUE)) %>%
  select(eid, cns_vasc_date, cns_infl_date, cns_atroph_date, cns_mov_date, 
         cns_demyel_date, cns_parox_date, cns_other_date) %>%
  rename(id = eid)
cns_dis$cns_vasc <- 0; cns_dis$cns_vasc[!is.na(cns_dis$cns_vasc_date)] <- 1
cns_dis$cns_vasc[cns_dis$cns_vasc_date %in% invalid_dates] <- NA
cns_dis$cns_vasc_date[cns_dis$cns_vasc_date %in% invalid_dates] <- NA
cns_dis$cns_infl <- 0; cns_dis$cns_infl[!is.na(cns_dis$cns_infl_date)] <- 1
cns_dis$cns_infl[cns_dis$cns_infl_date %in% invalid_dates] <- NA
cns_dis$cns_infl_date[cns_dis$cns_infl_date %in% invalid_dates] <- NA
cns_dis$cns_atroph <- 0; cns_dis$cns_atroph[!is.na(cns_dis$cns_atroph_date)] <- 1
cns_dis$cns_atroph[cns_dis$cns_atroph_date %in% invalid_dates] <- NA
cns_dis$cns_atroph_date[cns_dis$cns_atroph_date %in% invalid_dates] <- NA
cns_dis$cns_mov <- 0; cns_dis$cns_mov[!is.na(cns_dis$cns_mov_date)] <- 1
cns_dis$cns_demyel <- 0; cns_dis$cns_demyel[!is.na(cns_dis$cns_demyel_date)] <- 1
cns_dis$cns_demyel[cns_dis$cns_demyel_date %in% invalid_dates] <- NA
cns_dis$cns_demyel_date[cns_dis$cns_demyel_date %in% invalid_dates] <- NA
cns_dis$cns_parox <- 0; cns_dis$cns_parox[!is.na(cns_dis$cns_parox_date)] <- 1
cns_dis$cns_parox[cns_dis$cns_parox_date %in% invalid_dates] <- NA
cns_dis$cns_parox_date[cns_dis$cns_parox_date %in% invalid_dates] <- NA
cns_dis$cns_other <- 0; cns_dis$cns_other[!is.na(cns_dis$cns_other_date)] <- 1
cns_dis$cns_other[cns_dis$cns_other_date %in% invalid_dates] <- NA
cns_dis$cns_other_date[cns_dis$cns_other_date %in% invalid_dates] <- NA




## cancer: import and clean; will be used below
cancer <- data_all %>%
  select(c(eid, starts_with(c('X40013.', 'X40006.', 'X40005.'))))
cancer <- data.frame(sapply(cancer, as.character))
cancer[cancer == ''] <- NA
cancer <- cancer %>%
  mutate(across(starts_with('X40005'), ~as.Date(., format = '%Y-%m-%d')))

## CNS tumours: benign meningioma, meningeal cancer, brain cancer,
icd9_brain <- c('225', '2258', '2259', as.character(seq(2250, 2254)), 
                '1921', '191', '1910', '1911', '1912', as.character(seq(1916, 1919)))
icd10_brain <- c('D32', 'D320', 'D321', 'D329', 'D33', 'D331', 'D332', 'D333', 
                 'D334', 'D337', 'C70', 'C700', 'C709', 'C71', 'C710', 'C711',
                 'C712', 'C713', 'C714', 'C715', 'C716', 'C717', 'C718', 'C719', 'C793')
cancer_brain <- cancer_subtype(cancer, icd9_brain, icd10_brain, 'brain') %>%
  rename(cns_cancer_date = cancer_brain_date)
cancer_brain$cns_cancer <- 0
cancer_brain$cns_cancer[!is.na(cancer_brain$cns_cancer_date)] <- 1




## TBI
# there is no UKB variable, so we have to manually search for the codes
codes_tbi <- read.csv('tbi_codes.csv', header=TRUE, sep = ',')
# remove potential white space and convert to lower case
codes_tbi <- data.frame(sapply(codes_tbi, trimws))
codes_tbi$source <- tolower(codes_tbi$source)
codes_tbi$code <- as.character(codes_tbi$code)
codes_tbi$n <- NA

# get inpatient diagnoses
colnames(inpatient)[colnames(inpatient) == 'diagnosis'] <- 'code'
inpatient <- filter(inpatient, 
                    code %in% codes_tbi$code[codes_tbi$source == 'icd9'] | 
                      code %in% codes_tbi$code[codes_tbi$source == 'icd10'])  
inpatient$date <- as.Date(inpatient$date, format = '%Y-%m-%d')

# GP diagnoses
meds_diagnoses <- data.table::fread('gp_clinical.txt', sep='\t', header=TRUE, quote='')
meds_diagnoses <- as.data.frame(meds_diagnoses)
meds_diagnoses <- subset(meds_diagnoses, select = c(eid, data_provider, event_dt, read_2, read_3))
colnames(meds_diagnoses) <- c('id', 'data_provider', 'date_primary', 'read2', 'read3')
meds_diagnoses <- filter(meds_diagnoses, 
                         (read2 %in% codes_tbi$code[codes_tbi$source == 'read2']) | 
                           (read3 %in% codes_tbi$code[codes_tbi$source == 'read3']))
meds_diagnoses$date_primary <- as.Date(meds_diagnoses$date_primary, format = '%d/%m/%Y')
meds_diagnoses[meds_diagnoses == ''] <- NA
diagnoses_dates_inv <- filter(meds_diagnoses, date_primary %in% invalid_dates) %>%
  select(id, date_primary) %>% rename(date_inv = date_primary)
meds_diagnoses <- filter(meds_diagnoses, !date_primary %in% invalid_dates)

# remove duplicate codes and match codes with descriptions
inpatient <- inpatient %>% arrange(date)
inpatient <- distinct(inpatient, id, code, version, .keep_all = TRUE)
for (d in c('icd9', 'icd10')){
  for (diagnosis in codes_tbi$code[codes_tbi$source == d]){
    inpatient$description[inpatient$version == d & inpatient$code == diagnosis] <- 
      codes_tbi$description[codes_tbi$source == d & codes_tbi$code == diagnosis]
    inpatient$diagnosis[inpatient$version == d & inpatient$code == diagnosis] <- 
      codes_tbi$simple[codes_tbi$source == d & codes_tbi$code == diagnosis]
    codes_tbi$n[codes_tbi$source == d & codes_tbi$code == diagnosis] <- 
      length(inpatient$diagnosis[inpatient$version == d & inpatient$code == diagnosis])
  }
}

meds_diagnoses <- meds_diagnoses %>% arrange(date_primary)
meds_diagnoses <- distinct(meds_diagnoses, id, read2, .keep_all = TRUE)
meds_diagnoses <- distinct(meds_diagnoses, id, read3, .keep_all = TRUE)

for (d in c('read2', 'read3')){
  for (diagnosis in codes_tbi$code[codes_tbi$source == d]){
    meds_diagnoses$description[!is.na(meds_diagnoses[[d]]) & meds_diagnoses[[d]] == diagnosis] <- 
      codes_tbi$description[codes_tbi$source == d & codes_tbi$code == diagnosis]
    meds_diagnoses$diagnosis[!is.na(meds_diagnoses[[d]]) & meds_diagnoses[[d]] == diagnosis] <- 
      codes_tbi$simple[codes_tbi$source == d & codes_tbi$code == diagnosis]
    codes_tbi$n[codes_tbi$source == d & codes_tbi$code == diagnosis] <- 
      length(meds_diagnoses$diagnosis[!is.na(meds_diagnoses[[d]]) & meds_diagnoses[[d]] == diagnosis])
  }
}

# combine inpatient and gp diagnoses
read2 <- filter(meds_diagnoses, !is.na(read2))
read3 <- filter(meds_diagnoses, !is.na(read3))
read2 <- subset(read2, select = -c(read3))
read3 <- subset(read3, select = -c(read2))
read2$diag_source <- 'read2'; read3$diag_source <- 'read3'
colnames(read2) <- c('id', 'diag_data_provider', 'cns_tbi_date', 'diag_code', 
                     'diag_desc', 'diag_source')
colnames(read3) <- c('id', 'diag_data_provider', 'cns_tbi_date', 'diag_code', 
                     'diag_desc', 'diag_source')
meds_diagnoses <- rbind(subset(read2, select = c(id, cns_tbi_date, diag_code, 
                                                 diag_desc, diag_source, diag_data_provider)),
                        subset(read3, select = c(id, cns_tbi_date, diag_code, 
                                                 diag_desc, diag_source, diag_data_provider)))
inpatient$data_provider <- NA
colnames(inpatient) <- c('id', 'diag_code', 'cns_tbi_date', 'diag_source', 
                         'diag_desc', 'diag_data_provider')
tbi <- rbind(meds_diagnoses,
             subset(inpatient, select = c(id, cns_tbi_date, diag_code, diag_desc, 
                                          diag_source, diag_data_provider)))
tbi$diag <- 1
tbi <- tbi %>% 
  arrange(cns_tbi_date) %>%
  distinct(id, diag, .keep_all = TRUE) %>%
  rename(cns_tbi = diag) %>%
  select(id, cns_tbi, cns_tbi_date)
tbi <- merge(tbi, diagnoses_dates_inv, by = 'id', all = TRUE)
tbi$cns_tbi[is.na(tbi$cns_tbi_date) & !is.na(tbi$date_inv)] <- 99
tbi$date_inv <- NULL
rm(codes_tbi)

# create a single variable for all CNS disorders and determine the earliest date
cns_dis <- merge(cns_dis, cancer_brain, by = 'id', all = TRUE)
cns_dis <- merge(cns_dis, tbi, by = 'id', all = TRUE)
cns_dis$cns_tbi[is.na(cns_dis$cns_tbi)] <- 0
cns_dis$cns_tbi[cns_dis$cns_tbi == 99] <- NA
cns_dis <- cns_dis %>% 
  mutate(cns_any_date = reduce(across(c(cns_vasc_date, cns_infl_date, 
                                        cns_atroph_date, cns_mov_date,
                                        cns_demyel_date, cns_parox_date, 
                                        cns_other_date, cns_cancer_date, cns_tbi_date)), 
                                                    pmin, na.rm = TRUE)) %>%
  mutate(cns_any = ifelse((cns_vasc == 1 & !is.na(cns_vasc)) |
                            (cns_infl == 1 & !is.na(cns_infl)) |
                            (cns_atroph == 1 & !is.na(cns_atroph)) |
                            (cns_mov == 1 & !is.na(cns_mov)) |
                            (cns_demyel == 1 & !is.na(cns_demyel)) |
                            (cns_parox == 1 & !is.na(cns_parox)) |
                            (cns_other == 1 & !is.na(cns_other)) |
                            (cns_cancer == 1 & !is.na(cns_cancer)) |
                            (cns_tbi == 1 & !is.na(cns_tbi)), 1, 0))
# if all disorders are 0 and at least one is NA, set to NA
cols_to_check <- c('cns_vasc', 'cns_infl', 'cns_atroph', 'cns_mov', 'cns_demyel', 
                   'cns_parox', 'cns_other', 'cns_cancer', 'cns_tbi')
# create a new column that indicates whether a row contains a 1 among the selected columns
cns_dis$contains_1 <- apply(cns_dis[cols_to_check], 1, function(row) any(row == 1, na.rm = TRUE))
# create a new column that indicates whether a row contains an NA among the selected columns
cns_dis$contains_NA <- apply(cns_dis[cols_to_check], 1, function(row) any(is.na(row)))
# if a row contains a 1 among the selected columns, 'col_new' will keep its value;
# if a row contains an NA and no 1s among the selected columns, 'col_new' will be changed to NA
cns_dis$cns_any[cns_dis$contains_NA == TRUE & cns_dis$contains_1 == FALSE] <- NA
cns_dis$contains_1 <- NULL
cns_dis$contains_NA <- NULL
# merge with TBI
rm(meds_diagnoses, inpatient, read2, read3, diagnoses_dates_inv, cancer_brain,
   tbi)
gc()




## endocrine disorders
endocrine_dis <- data_all %>%
  select(c(eid, starts_with(c('X130690.', 'X130692.', 'X130694.', 'X130696.', 
                              'X130698.', 'X130700.', 'X130702.', 'X130704.', 
                              'X130718.', 'X130720.', 'X130722.', 'X130724.', 
                              'X130726.', 'X130728.', 'X130730.', 'X130732.', 
                              'X130734.', 'X130736.', 'X130738',  'X130742.', 
                              'X130744.', 'X130746.', 'X130748')))) %>%
  mutate(across(starts_with('X'), ~as.Date(., format = '%Y-%m-%d'))) %>%
  mutate(endocrine_dis_date = reduce(across(starts_with('X')), pmin, na.rm = TRUE)) %>%
  select(eid, endocrine_dis_date) %>%
  rename(id = eid)
endocrine_dis$endocrine_dis <- 0
endocrine_dis$endocrine_dis[!is.na(endocrine_dis$endocrine_dis_date)] <- 1
endocrine_dates_inv <- filter(endocrine_dis, endocrine_dis_date %in% invalid_dates) %>%
  select(id, endocrine_dis_date) %>% rename(date_inv = endocrine_dis_date)
endocrine_dis <- filter(endocrine_dis, !endocrine_dis_date %in% invalid_dates)
# carcinoid syndrome - cancers of the endocrine system:
icd10_endocrine <- c('C25', 'C250', 'C251', 'C252', 'C253', 'C254', 'C257', 
                     'C258', 'C259', 'C73', 'C74', 'C740', 'C741', 'C749', 
                     'C75', 'C750', 'C751', 'C753', 'C754', 'C755', 'D093', 
                     'D34', 'D35', 'D351', 'D352', 'D44', 'D440', 'D441', 'D442', 
                     'D443', 'D444', 'D445', 'D446', 'D447', 'D448', 'D449')
icd9_endocrine <- c('157', '1570', '1572', '1574', '1579', '193', '194', '1940', 
                    '1941', '1943', '1944', '1949', '226', '227', '2273', '2279', 
                    '2370')
cancer_endocrine <- cancer_subtype(cancer, icd9_endocrine, icd10_endocrine, 'endocrine')
cancer_endocrine$cancer_endocrine <- 0
cancer_endocrine$cancer_endocrine[!is.na(cancer_endocrine$cancer_endocrine_date)] <- 1
# combine endocrine disorders with endocrine cancers
endocrine_dis <- merge(endocrine_dis, cancer_endocrine, by = 'id', all = TRUE)
endocrine_dis$cancer_endocrine[is.na(endocrine_dis$cancer_endocrine)] <- 0
endocrine_dis$endocrine_dis[endocrine_dis$cancer_endocrine == 1] <- 1
endocrine_dis <- endocrine_dis %>%
  mutate(endocrine_dis_date = reduce(across(c(endocrine_dis_date, cancer_endocrine_date)), 
                                     pmin, na.rm = TRUE))
endocrine_dis <- subset(endocrine_dis, select = -c(cancer_endocrine_date, cancer_endocrine))
endocrine_dis <- merge(endocrine_dis, endocrine_dates_inv, by = 'id', all = TRUE)
endocrine_dis$endocrine_dis[is.na(endocrine_dis$endocrine_dis_date) & !is.na(endocrine_dis$date_inv)] <- NA
endocrine_dis$date_inv <- NULL
rm(cancer_endocrine, endocrine_dates_inv)




## metabolic disorders
metabolic_dis <- data_all %>%
  select(c(eid, starts_with(c('X130798.', 'X130800.', 'X130802.', 'X130806.', 
                              'X130808.', 'X130810.', 'X130812.', 'X130816.', 
                              'X130818.', 'X130820.', 'X130822.', 'X130824.', 
                              'X130826.', 'X130828.', 'X130830.', 'X130832.')))) %>%
  mutate(across(starts_with('X'), ~as.Date(., format = '%Y-%m-%d'))) %>%
  mutate(metabolic_dis_date = reduce(across(starts_with('X')), pmin, na.rm = TRUE)) %>%
  select(eid, metabolic_dis_date) %>%
  rename(id = eid)
metabolic_dis$metabolic_dis <- 0
metabolic_dis$metabolic_dis[!is.na(metabolic_dis$metabolic_dis_date)] <- 1
metabolic_dates_inv <- filter(metabolic_dis, metabolic_dis_date %in% invalid_dates) %>%
  select(id, metabolic_dis_date) %>% rename(date_inv = metabolic_dis_date)
metabolic_dis <- filter(metabolic_dis, !metabolic_dis_date %in% invalid_dates)
# include hyperlipidaema
metabolic_dis <- merge(metabolic_dis, hyperchol, by = 'id', all.x = TRUE) %>%
  mutate(metabolic_dis_date = reduce(across(c(metabolic_dis_date, hyperlip_date)), 
                                     pmin, na.rm = TRUE))
metabolic_dis$metabolic_dis[metabolic_dis$hyperlip == 1] <- 1
metabolic_dis <- subset(metabolic_dis, select = -c(hyperlip_date, hyperlip))
metabolic_dis <- merge(metabolic_dis, metabolic_dates_inv, by = 'id', all = TRUE)
metabolic_dis$metabolic_dis[is.na(metabolic_dis$metabolic_dis_date) & !is.na(metabolic_dis$date_inv)] <- NA
metabolic_dis$date_inv <- NULL
rm(metabolic_dates_inv)



## various other disorders
outcomes <- data_all %>% 
  # respiratory disease
  select(c(eid, starts_with(c('X131484.', 'X131486.', 'X131488.', 'X131490.', 
                              'X131492.', 'X131494.', 'X131496.',
                              # cerebrovascular disease
                              'X131360.', 'X131362.', 'X131364.', 'X131366.', 
                              'X131368.', 'X131370.', 'X131372.', 'X131374.', 
                              'X131376.', 'X131378.',
                              # liver disease
                              'X131498.', 'X131658.', 'X131660.', 'X131662.', 
                              'X131664.', 'X131666.', 'X131668.', 'X131670.',
                              # influenza
                              'X131438.', 'X131440.', 'X131442.', 'X131444.', 
                              'X131446.', 'X131448.', 'X131450.', 'X131452.', 
                              'X131454.', 'X131456.', 'X131296.', 'X131298.', 
                              'X131300.', 'X131302.', 'X131304.', 'X131306.')))) %>%
  mutate(across(starts_with('X'), ~ as.Date(., format = '%Y-%m-%d')))  %>%
  mutate(respiratory_date = 
           reduce(across(c('X131484.0.0', 'X131486.0.0', 'X131488.0.0', 
                           'X131490.0.0', 'X131492.0.0', 'X131494.0.0', 
                           'X131496.0.0', 'X131498.0.0')), 
                  pmin, na.rm = TRUE)) %>%
  mutate(cerebrovascular_date = reduce(across(c('X131360.0.0', 'X131362.0.0', 
                                                'X131364.0.0', 'X131366.0.0', 
                                                'X131368.0.0', 'X131370.0.0', 
                                                'X131372.0.0', 'X131374.0.0', 
                                                'X131376.0.0', 'X131378.0.0')),
                                       pmin, na.rm = TRUE)) %>%
  mutate(hepatic_date = reduce(across(c('X131658.0.0', 'X131660.0.0', 'X131662.0.0', 
                                        'X131664.0.0', 'X131666.0.0', 'X131668.0.0', 
                                        'X131670.0.0')), 
                               pmin, na.rm = TRUE)) %>%
  mutate(flu_date = reduce(across(c('X131438.0.0', 'X131440.0.0', 'X131442.0.0', 
                                    'X131444.0.0', 'X131446.0.0', 'X131448.0.0', 
                                    'X131450.0.0', 'X131452.0.0', 'X131454.0.0', 
                                    'X131456.0.0')), pmin, na.rm = TRUE)) %>%
  mutate(heart_date = reduce(across(c('X131296.0.0', 'X131298.0.0', 'X131300.0.0', 
                                      'X131302.0.0', 'X131304.0.0', 'X131306.0.0')), 
                             pmin, na.rm = TRUE)) %>%
  rename(id = eid) %>%
  select(id, respiratory_date, cerebrovascular_date, hepatic_date, flu_date, heart_date) %>%  
  filter(rowSums(is.na(select(., -id))) != ncol(.) - 1)
outcomes$respiratory <- 0; outcomes$respiratory[!is.na(outcomes$respiratory_date)] <- 1
outcomes$respiratory[outcomes$respiratory_date %in% invalid_dates] <- 99
outcomes$respiratory_date[outcomes$respiratory_date %in% invalid_dates] <- NA
outcomes$cerebrovascular <- 0; outcomes$cerebrovascular[!is.na(outcomes$cerebrovascular_date)] <- 1
outcomes$cerebrovascular[outcomes$cerebrovascular_date %in% invalid_dates] <- 99
outcomes$cerebrovascular_date[outcomes$cerebrovascular_date %in% invalid_dates] <- NA
outcomes$hepatic <- 0; outcomes$hepatic[!is.na(outcomes$hepatic_date)] <- 1
outcomes$hepatic[outcomes$hepatic_date %in% invalid_dates] <- 99
outcomes$hepatic_date[outcomes$hepatic_date %in% invalid_dates] <- NA
outcomes$flu <- 0; outcomes$flu[!is.na(outcomes$flu_date)] <- 1
outcomes$flu[outcomes$flu_date %in% invalid_dates] <- 99
outcomes$flu_date[outcomes$flu_date %in% invalid_dates] <- NA
outcomes$heart <- 0; outcomes$heart[!is.na(outcomes$heart_date)] <- 1
outcomes$heart[outcomes$heart_date %in% invalid_dates] <- 99
outcomes$heart_date[outcomes$heart_date %in% invalid_dates] <- NA




## nutritional deficiencies
nutr_dis <- data_all %>%
  select(c(eid, starts_with(c('X130750.', 'X130752.', 'X130756.', 'X130758.', 
                              'X130760.', 'X130762.', 'X130764.', 'X130766.', 
                              'X130768.', 'X130770.', 'X130772.', 'X130774.', 
                              'X130776.', 'X130778.', 'X130780.', 'X130782.', 
                              'X130784.', 'X130786.', 'X130788')))) %>%
  mutate(across(starts_with('X'), ~as.Date(., format = '%Y-%m-%d'))) %>%
  mutate(nutr_dis_date = reduce(across(starts_with('X')), pmin, na.rm = TRUE)) %>%
  select(eid, nutr_dis_date) %>%
  rename(id = eid)
nutr_dis$nutr_dis <- 0; nutr_dis$nutr_dis[!is.na(nutr_dis$nutr_dis_date)] <- 1




## colon cancer, prostate-, lung-, breast-, and ovarian cancer
icd9_colon <- c('153', as.character(seq(1530, 1539)), '154', 
                as.character(seq(1540, 1543)), '1548', '2303', '2304')
icd10_colon <- c('C18', 'C180', 'C181', 'C182', 'C183', 'C184', 'C185', 'C186', 
                 'C187', 'C188', 'C189', 'C19', 'C20', 'C21', 'C210', 'C211', 
                 'C212', 'C218')
cancer_colon <- cancer_subtype(cancer, icd9_colon, icd10_colon, 'colon')
icd9_prostate <- c('185', '2365')
icd10_prostate <- c('C61')
cancer_prostate <- cancer_subtype(cancer, icd9_prostate, icd10_prostate, 'prostate')
icd9_lung <- c('162', '1620', '1622', '1623', '1624', '1625', '1628', '1629', '2357')
icd10_lung <- c('C33', 'C34', 'C340', 'C341', 'C342', 'C343', 'C348', 'C349')
cancer_lung <- cancer_subtype(cancer, icd9_lung, icd10_lung, 'lung')
icd9_breast <- c('174', as.character(seq(1740, 1746)), '1748', '1749', '2330')
icd10_breast <- c('C50', 'C500', 'C501', 'C502', 'C503', 'C504', 'C505', 'C506', 
                  'C508', 'C509')
cancer_breast <- cancer_subtype(cancer, icd9_breast, icd10_breast, 'breast')
icd9_ovary <- c('183', '1830', '1832', '184', '1840', '1841', '1844', '1848', 
                '1849', '1986', '2362')
icd10_ovary <- c('C56')
cancer_ovary <- cancer_subtype(cancer, icd9_ovary, icd10_ovary, 'ovary')
cancer_d <- merge(cancer_colon, cancer_prostate, by = 'id', all = TRUE)
cancer_d <- merge(cancer_d, cancer_lung, by = 'id', all = TRUE)
cancer_d <- merge(cancer_d, cancer_breast, by = 'id', all = TRUE)
cancer_d <- merge(cancer_d, cancer_ovary, by = 'id', all = TRUE)
cancer_d$cancer_colon <- 0
cancer_d$cancer_colon[!is.na(cancer_d$cancer_colon_date)] <- 1
cancer_d$cancer_prostate <- 0
cancer_d$cancer_prostate[!is.na(cancer_d$cancer_prostate_date)] <- 1
cancer_d$cancer_lung <- 0
cancer_d$cancer_lung[!is.na(cancer_d$cancer_lung_date)] <- 1
cancer_d$cancer_breast <- 0
cancer_d$cancer_breast[!is.na(cancer_d$cancer_breast_date)] <- 1
cancer_d$cancer_ovary <- 0
cancer_d$cancer_ovary[!is.na(cancer_d$cancer_ovary_date)] <- 1
cancer_d$cancer_prostate_ovary <- 0
cancer_d$cancer_prostate_ovary[cancer_d$cancer_prostate == 1 |
                                 cancer_d$cancer_ovary == 1] <- 1
rm(cancer_colon, icd9_colon, icd10_colon, cancer_prostate, icd9_prostate, 
   icd10_prostate, cancer_lung, icd9_lung, icd10_lung, cancer_breast, icd9_breast, 
   icd10_breast, cancer_ovary, icd9_ovary, icd10_ovary, cancer)
gc()



## hearing loss (adapted from Mur et al., 2024;
# GitHub: https://github.com/JuM24/HA-and-dementia-in-UKBB)
hear <- data_all %>%
  select(eid, starts_with(c('X2247.', 'X2257.', 'X20019.', 'X20021.', 
                            'X131258.', 'X131259.', 'X131260.', 'X131261.',
                            'X53.')))
hear$eid <- as.character(hear$eid)

# change colnames and code emtpy strings as NAs
colnames(hear) <- c('id', 'hear_dif_0', 'hear_dif_1', 'hear_dif_2', 'hear_dif_3', 
                    'hear_difn_0', 'hear_difn_1', 'hear_difn_2', 'hear_difn_3',
                    'srt_r_0', 'srt_r_1', 'srt_r_2', 'srt_r_3', 'srt_l_0', 
                    'srt_l_1', 'srt_l_2', 'srt_l_3', 'date_hear_loss_a',
                    'source_hear_loss_a', 'date_hear_loss_b', 'source_hear_loss_b',
                    'date_0', 'date_1', 'date_2', 'date_3')
hear[hear == ''] <- NA


# people who didn't know or didn't want to answer are coded as NA for all hearing-related self-report data
hear$hear_dif_0[hear$hear_dif_0 == -1 | hear$hear_dif_0 == -3] <- NA
hear$hear_dif_1[hear$hear_dif_1 == -1 | hear$hear_dif_1 == -3] <- NA
hear$hear_dif_2[hear$hear_dif_2 == -1 | hear$hear_dif_2 == -3] <- NA
hear$hear_dif_3[hear$hear_dif_3 == -1 | hear$hear_dif_3 == -3] <- NA

hear$hear_difn_0[hear$hear_difn_0 == -1 | hear$hear_difn_0 == -3] <- NA
hear$hear_difn_1[hear$hear_difn_1 == -1 | hear$hear_difn_1 == -3] <- NA
hear$hear_difn_2[hear$hear_difn_2 == -1 | hear$hear_difn_2 == -3] <- NA
hear$hear_difn_3[hear$hear_difn_3 == -1 | hear$hear_difn_3 == -3] <- NA

# This imports the medical hospital and primary-care data and searches 
# for the relevant codes in each. It then combines both sources to create 
# a single source of diagnoses, keeping just the earliest date of diagnosis.

# diagnosis codes
diagnosis_codes <- read.csv('hearing_codes.csv') # EHR codes for HL and HA
# remove duplicate codes within one source (so within ICD10, or ICD9, etc.)
diagnosis_codes <- distinct(diagnosis_codes, code, source, .keep_all = TRUE)

# hearing loss and hearing aid ascertainment
# inpatient diagnoses
inpatient <- readRDS('output_files/inpatient.Rds') # field IDs 41270, 41271, 41280, and 41281
colnames(inpatient)[colnames(inpatient) == 'diagnosis'] <- 'code'
# keep only relevant diagnoses
inpatient <- filter(inpatient, code %in% diagnosis_codes$code[diagnosis_codes$source == 'icd9'] | 
                      code %in% diagnosis_codes$code[diagnosis_codes$source == 'icd10']) %>%
  # replace invalid dates with NAs
  mutate(across(date, ~replace(., . %in% c('1900-01-01', '1901-01-01', 
                                           '1902-02-02', '1903-03-03', '2037-07-07'), NA))) 
inpatient$date <- as.Date(inpatient$date, format = '%Y-%m-%d')

# GP diagnoses
gp_diagnoses <- data.table::fread('gp_clinical.txt', sep='\t', header=TRUE, quote='') # field ID 42040
gp_diagnoses <- as.data.frame(gp_diagnoses)
gp_diagnoses <- subset(gp_diagnoses, select = c(eid, data_provider, event_dt, 
                                                read_2, read_3))
colnames(gp_diagnoses) <- c('id', 'data_provider', 'date_primary', 'read2', 'read3')
gp_diagnoses <- filter(gp_diagnoses, 
                       (read2 %in% diagnosis_codes$code[diagnosis_codes$source == 'read2']) | 
                         (read3 %in% diagnosis_codes$code[diagnosis_codes$source == 'read3'])) %>%
  # replace invalid dates with NAs
  mutate(across(date_primary, ~replace(., . %in% c('01/01/1900', '01/01/1901', 
                                                   '02/02/1902', '03/03/1903', 
                                                   '07/07/2037'), NA))) 
gp_diagnoses$date_primary <- as.Date(gp_diagnoses$date_primary, format = '%d/%m/%Y')
gp_diagnoses$year <- as.numeric(format(as.Date(gp_diagnoses$date_primary), '%Y'))
gp_diagnoses[gp_diagnoses == ''] <- NA

# remove duplicate codes and match codes with descriptions
inpatient <- inpatient %>% arrange(date)
inpatient <- distinct(inpatient, id, code, version, .keep_all = TRUE)
for (d in c('icd9', 'icd10')){
  for (diagnosis in diagnosis_codes$code[diagnosis_codes$source == d]){
    inpatient$description[inpatient$version == d & inpatient$code == diagnosis] <-
      diagnosis_codes$description[diagnosis_codes$source == d & diagnosis_codes$code == diagnosis]
    inpatient$diagnosis[inpatient$version == d & inpatient$code == diagnosis] <- 
      diagnosis_codes$simple[diagnosis_codes$source == d & diagnosis_codes$code == diagnosis]
    diagnosis_codes$n[diagnosis_codes$source == d & diagnosis_codes$code == diagnosis] <- 
      length(inpatient$diagnosis[inpatient$version == d & inpatient$code == diagnosis])
  }
}

# repeat for primary care
gp_diagnoses <- gp_diagnoses %>% arrange(date_primary)
gp_diagnoses <- distinct(gp_diagnoses, id, read2, .keep_all = TRUE)
gp_diagnoses <- distinct(gp_diagnoses, id, read3, .keep_all = TRUE)
for (d in c('read2', 'read3')){
  for (diagnosis in diagnosis_codes$code[diagnosis_codes$source == d]){
    gp_diagnoses$description[!is.na(gp_diagnoses[[d]]) & gp_diagnoses[[d]] == diagnosis] <- 
      diagnosis_codes$description[diagnosis_codes$source == d & diagnosis_codes$code == diagnosis]
    gp_diagnoses$diagnosis[!is.na(gp_diagnoses[[d]]) & gp_diagnoses[[d]] == diagnosis] <- 
      diagnosis_codes$simple[diagnosis_codes$source == d & diagnosis_codes$code == diagnosis]
    diagnosis_codes$n[diagnosis_codes$source == d & diagnosis_codes$code == diagnosis] <- 
      length(gp_diagnoses$diagnosis[!is.na(gp_diagnoses[[d]]) & gp_diagnoses[[d]] == diagnosis])
  }
}

# combine inpatient and gp diagnoses and export
# create two separate files for read2 and read3 code from the GP record
read2 <- filter(gp_diagnoses, !is.na(read2)) 
read3 <- filter(gp_diagnoses, !is.na(read3))
read2 <- subset(read2, select = -c(read3))
read3 <- subset(read3, select = -c(read2))
read2$diag_source <- 'read2'; read3$diag_source <- 'read3'
colnames(read2) <- c('id', 'diag_data_provider', 'diag_date', 'diag_code', 
                     'diag_year', 'diag_desc', 'diag', 'diag_source')
colnames(read3) <- c('id', 'diag_data_provider', 'diag_date', 'diag_code', 
                     'diag_year', 'diag_desc', 'diag', 'diag_source')
# re-combine read2 and read3
gp_diagnoses <- rbind(subset(read2, select = c(id, diag, diag_date, diag_year, diag_code, 
                                               diag_desc, diag_source, diag_data_provider)), 
                      subset(read3, select = c(id, diag, diag_date, diag_year, diag_code, 
                                               diag_desc, diag_source, diag_data_provider)))
inpatient$year <- as.numeric(format(as.Date(inpatient$date), '%Y'))
inpatient$data_provider <- NA
colnames(inpatient) <- c('id', 'diag_code', 'diag_date', 'diag_source', 'diag_desc', 
                         'diag', 'diag_year', 'diag_data_provider')
diagnoses <- rbind(gp_diagnoses,
                   subset(inpatient, select = c(id, diag, diag_date, diag_year, diag_code, 
                                                diag_desc, diag_source, diag_data_provider)))

# remove duplicate diagnoses of the same type and check for multiple types for same participant
diagnoses <- arrange(diagnoses, diag_date)
diagnoses <- distinct(diagnoses, id, diag, .keep_all = TRUE)

# separate into data frames of distinct diagnoses, remove duplicates, and bind again; 
# then merge with main data frame; tag the ones without dates (they are going to be removed later)
# do this for hearing loss, hearing aid use, hearing aid use cessation, and cochlear implants
diagnoses <- distinct(diagnoses, id, .keep_all = TRUE) %>%
  rename(hear_loss_code = diag_code, hear_loss_diag = diag, hear_loss_date = diag_date, 
         hear_loss_year = diag_year, hear_loss_desc = diag_desc,
         hear_loss_source = diag_source, hear_loss_data_provider = diag_data_provider)
diagnoses$hear_loss <- 1
diagnoses$hear_loss_nodate <- 0
diagnoses$hear_loss_nodate[is.na(diagnoses$hear_loss_date)] <- 1

# merge with the hearing masterfile
hear <- merge(hear, diagnoses, by = 'id', all.x = TRUE)

# for now, all participants without an explicit diagnosis are considered undiagnosed; 
# diagnoses without dates will be removed later using the 'nodate' variables created above
hear$hear_loss[is.na(hear$hear_loss)] <- 0 

hear$hear_dif_both_0 <- rowSums(hear[ , c('hear_dif_0', 'hear_difn_0')], na.rm = TRUE)
# NA if both of them are NA
hear$hear_dif_both_0[is.na(hear$hear_dif_0) & is.na(hear$hear_difn_0)] <- NA
# NA if one of them is 1 and the other one is NA 
# (because they would be classified as having hearing loss if the other one was 1)
hear$hear_dif_both_0[is.na(hear$hear_dif_0) & hear$hear_difn_0 == 1] <- NA 
hear$hear_dif_both_0[hear$hear_dif_0 == 1 & is.na(hear$hear_difn_0)] <- NA

hear$hear_dif_both_1 <- rowSums(hear[ , c('hear_dif_1', 'hear_difn_1')], na.rm = TRUE)
hear$hear_dif_both_1[is.na(hear$hear_dif_1) & is.na(hear$hear_difn_1)] <- NA
hear$hear_dif_both_1[is.na(hear$hear_dif_1) & hear$hear_difn_1 == 1] <- NA
hear$hear_dif_both_1[hear$hear_dif_1 == 1 & is.na(hear$hear_difn_1)] <- NA

hear$hear_dif_both_2 <- rowSums(hear[ , c('hear_dif_2', 'hear_difn_2')], na.rm = TRUE)
hear$hear_dif_both_2[is.na(hear$hear_dif_2) & is.na(hear$hear_difn_2)] <- NA
hear$hear_dif_both_2[is.na(hear$hear_dif_2) & hear$hear_difn_2 == 1] <- NA
hear$hear_dif_both_2[hear$hear_dif_2 == 1 & is.na(hear$hear_difn_2)] <- NA

hear$hear_dif_both_3 <- rowSums(hear[ , c('hear_dif_3', 'hear_difn_3')], na.rm = TRUE)
hear$hear_dif_both_3[is.na(hear$hear_dif_3) & is.na(hear$hear_difn_3)] <- NA
hear$hear_dif_both_3[is.na(hear$hear_dif_3) & hear$hear_difn_3 == 1] <- NA
hear$hear_dif_both_3[hear$hear_dif_3 == 1 & is.na(hear$hear_difn_3)] <- NA

# add the diagnoses identified by UKB 'first occurrences ID fields' 
# but unidentified by my search of the EHR; get earliest of both UKB dates
hear$hear_loss_a <- 0; hear$hear_loss_a[!is.na(hear$date_hear_loss_a)] <- 1
hear$hear_loss_b <- 0; hear$hear_loss_b[!is.na(hear$date_hear_loss_b)] <- 1
hear$date_hear_loss_a[as.character(hear$date_hear_loss_a) %in% 
                        c('1900-01-01', '1901-01-01', '2037-07-07', 
                          '1902-02-02', '1903-03-03')] <- NA # invalid dates to NA
hear$date_hear_loss_b[as.character(hear$date_hear_loss_b) %in% 
                        c('1900-01-01', '1901-01-01', '2037-07-07', 
                          '1902-02-02', '1903-03-03')] <- NA
hear$date_hear_loss_a <- as.Date(hear$date_hear_loss_a, format = '%Y-%m-%d')
hear$date_hear_loss_b <- as.Date(hear$date_hear_loss_b, format = '%Y-%m-%d')
# create 'nodate' variable that indicates lack of date of diagnosis
hear$hear_loss_a_nodate[hear$hear_loss_a == 1] <- 0 
hear$hear_loss_b_nodate[hear$hear_loss_b == 1] <- 0
hear$hear_loss_a_nodate[hear$hear_loss_a == 1 & is.na(hear$date_hear_loss_a)] <- 1
hear$hear_loss_b_nodate[hear$hear_loss_b == 1 & is.na(hear$date_hear_loss_b)] <- 1

# combine both types of hearing loss from first occurrences; 
# determine if in those with loss at least one date is not missing
# if at least one of them has a date, that it's fine
hear$hear_loss_ab_nodate[hear$hear_loss_a == 1 | hear$hear_loss_b == 1] <- 0 
# if the one we have a diagnosis for doesn't have a date, that's not fine
hear$hear_loss_ab_nodate[(hear$hear_loss_a_nodate == 1 & is.na(hear$hear_loss_b_nodate)) | 
                           (hear$hear_loss_b_nodate == 1 & is.na(hear$hear_loss_a_nodate)) |
                           (hear$hear_loss_a_nodate == 1 & hear$hear_loss_b_nodate == 1)] <- 1

# get earliest of both UKB dates and set that as the source of the diagnosis
hear <- hear %>%
  mutate(date_hear_loss_ab = pmin(date_hear_loss_a, date_hear_loss_b, na.rm = TRUE),
         source_hear_loss_ab = ifelse((date_hear_loss_a < date_hear_loss_b & 
                                         !is.na(date_hear_loss_a)) | 
                                        is.na(date_hear_loss_b), 
                                      source_hear_loss_a, source_hear_loss_b))

# now let's include the objective hearing assessment (SRT)
# create a new variable indicating SRT for 'better' ear for each visit (lower score means better hearing)
hear <- transform(hear, srt_min_0 = pmin(srt_r_0, srt_l_0, na.rm = TRUE))
hear <- transform(hear, srt_min_1 = pmin(srt_r_1, srt_l_1, na.rm = TRUE))
hear <- transform(hear, srt_min_2 = pmin(srt_r_2, srt_l_2, na.rm = TRUE))
hear <- transform(hear, srt_min_3 = pmin(srt_r_3, srt_l_3, na.rm = TRUE))

# new variable indicating hearing problems according to any of our criteria: 
# (1) self-report (2 indicates affirmative answers to both questions, 
# 99 indicates deafness),(2) hearing loss acc. to our search of the EHR, 
# (3) SRT, and (4) first-occurrences variables in UKB
hear <- hear %>%
  mutate(hear_loss_any = ifelse((!is.na(hear_dif_both_0) & hear_dif_both_0 %in% c(2, 99)) | 
                                  (!is.na(hear_dif_both_1) & hear_dif_both_1 %in% c(2, 99)) | 
                                  (!is.na(hear_dif_both_2) & hear_dif_both_2 %in% c(2, 99)) | 
                                  (!is.na(hear_dif_both_3) & hear_dif_both_3 %in% c(2, 99)) | 
                                  (!is.na(hear_loss) & hear_loss == 1) | 
                                  (!is.na(srt_min_0) & srt_min_0 > -5.5) | 
                                  (!is.na(srt_min_1) & srt_min_1 > -5.5) | 
                                  (!is.na(srt_min_2) & srt_min_2 > -5.5) | 
                                  (!is.na(srt_min_3) & srt_min_3 > -5.5) | 
                                  (!is.na(hear_loss_a) & hear_loss_a == 1) |
                                  (!is.na(hear_loss_b) & hear_loss_b == 1), 1, 0))

# assign date to hearing loss corresponding to date of assessment if hearing loss 
# was established at that assessment (i.e., through self-report or SRT)
# assessment 0
hear$date_hear_loss_0[((hear$hear_dif_both_0 == 2 | hear$hear_dif_both_0 == 99) & 
                         !is.na(hear$hear_dif_both_0)) | (hear$srt_min_0 > -5.5  
                                                          & !is.na(hear$srt_min_0))] <- 
  hear$date_0[((hear$hear_dif_both_0 == 2 | hear$hear_dif_both_0 == 99) & 
                 !is.na(hear$hear_dif_both_0)) | (hear$srt_min_0 > -5.5 
                                                  & !is.na(hear$srt_min_0))]
hear$date_hear_loss_0 <- as.Date(hear$date_hear_loss_0)
# assessment 1
hear$date_hear_loss_1[((hear$hear_dif_both_1 == 2 | hear$hear_dif_both_1 == 99) & 
                         !is.na(hear$hear_dif_both_1)) | (hear$srt_min_1 > -5.5  
                                                          & !is.na(hear$srt_min_1))] <- 
  hear$date_1[((hear$hear_dif_both_1 == 2 | hear$hear_dif_both_1 == 99) & 
                 !is.na(hear$hear_dif_both_1)) | (hear$srt_min_1 > -5.5 & !is.na(hear$srt_min_1))]
hear$date_hear_loss_1 <- as.Date(hear$date_hear_loss_1)
# assessment 2
hear$date_hear_loss_2[((hear$hear_dif_both_2 == 2 | hear$hear_dif_both_2 == 99) & 
                         !is.na(hear$hear_dif_both_2)) | (hear$srt_min_2 > -5.5  
                                                          & !is.na(hear$srt_min_2))] <- 
  hear$date_2[((hear$hear_dif_both_2 == 2 | hear$hear_dif_both_2 == 99) & 
                 !is.na(hear$hear_dif_both_2)) | (hear$srt_min_2 > -5.5 & !is.na(hear$srt_min_2))]
hear$date_hear_loss_2 <- as.Date(hear$date_hear_loss_2)
# assessment 3 
hear$date_hear_loss_3[((hear$hear_dif_both_3 == 2 | hear$hear_dif_both_3 == 99) & 
                         !is.na(hear$hear_dif_both_3)) | (hear$srt_min_3 > -5.5  
                                                          & !is.na(hear$srt_min_3))] <- 
  hear$date_3[((hear$hear_dif_both_3 == 2 | hear$hear_dif_both_3 == 99) & 
                 !is.na(hear$hear_dif_both_3)) | (hear$srt_min_3 > -5.5 & !is.na(hear$srt_min_3))]
hear$date_hear_loss_3 <- as.Date(hear$date_hear_loss_3)

# find earliest hearing loss date among all sources
hear <- transform(hear, date_hear_loss_any = 
                    pmin(date_hear_loss_0, date_hear_loss_1, date_hear_loss_2, 
                         date_hear_loss_3, hear_loss_date, date_hear_loss_ab, na.rm = TRUE)) %>%
  select(id, hear_loss_any, date_hear_loss_any)
rm(diagnoses, diagnosis_codes, gp_diagnoses, read2, read3, inpatient)
gc()




## merge all variables
covs <- Reduce(function(x, y) merge(x, y, by = 'id', all = TRUE), 
               list(alcohol, cognition, dems, deprivation,
                    education, phys_act, pollution, smoking,
                    social, cancer_d, cns_dis, death, dementia, diabetes, hear, 
                    hyperchol, hypertension,  metabolic_dis, mood_ado, 
                    nutr_dis, outcomes, psych_ado, sleep_vision,
                    endocrine_dis)) %>%
  # NAs to 0s
  mutate(across(c(vision_problem, sleep_dis_any, cerebrovascular,
                  respiratory, hepatic, flu, heart, cancer_colon, 
                  cancer_prostate, cancer_lung, cancer_breast, cancer_ovary,
                  cns_cancer, cns_any), 
                ~ replace_na(., 0)))
# 99 dummy value to NA
covs[covs == 99] <- NA

covs$id <- as.character(covs$id)

# export and clear environment
saveRDS(covs, 'output_files/covariates.Rds')
rm(list = ls())
gc()