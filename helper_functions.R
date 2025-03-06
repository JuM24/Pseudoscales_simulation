# outlier-removal function
outliers <- function(x, var_metric, method) {
  if (method == 'SD'){
    maximum <- (mean(x, na.rm=T)) + (var_metric * sd(x, na.rm=T))
    minimum <- (mean(x, na.rm=T)) - (var_metric * sd(x, na.rm=T))
  }
  else if (method == 'IQR'){
    maximum <- (quantile(x, 0.75, na.rm=T)) + (var_metric * IQR(x, na.rm=T))
    minimum <- (quantile(x, 0.25, na.rm=T)) - (var_metric * IQR(x, na.rm=T))
  }
  x[(x > maximum) | (x < minimum)] <- NA
  return(x)
}



# for each row, moves all NAs to the right of the data frame
move_nas_to_right <- function(x) {
  non_nas <- x[!is.na(x)]
  n_nas <- sum(is.na(x))
  return(c(non_nas, rep(NA, n_nas)))
}



# replace values not in the icd codes list with 'DUMMY'
replace_not_in_icd9 <- function(col, icd9) {
  col[!col %in% icd9 & !is.na(col)] <- 'DUMMY'
  return(col)
}
replace_not_in_icd10 <- function(col, icd10) {
  col[!col %in% icd10 & !is.na(col)] <- 'DUMMY'
  return(col)
}



## classification of physical activity
# 0 - (none, low [light household tasks only], 
# 1 - medium [heavy household tasks and/or walking for pleasure and/or other exercise], 
# 2 - high [strenuous sports]).
phys_act_classify <- function(x){
  if (any(x == 3 & !is.na(x), na.rm = TRUE)){
    return(3)
  } else if (any((x == 1 | x == 2 | x == 5) & !is.na(x), na.rm = TRUE)){
    return(2)
  } else if(any(x == 4 & !is.na(x), na.rm = TRUE)){
    return(1)
  } else if (any (x == -7  & !is.na(x), na.rm = TRUE)){
    return(0)
  } else {
    return(NA)
  }
}




# takes in main disorder data frame (produced below), a list of icd9 and icd10 
# codes for the relevant cancer subtype, and the name of the subtype (e.g., 'breast')
# it outputs the data frame with the earliest diagnoses of that subtype
cancer_subtype <- function(df, icd9, icd10, subtype){
  # separate X40013 and X40006 from X40005
  subtype_dates <- df %>% 
    select(eid, starts_with('X40005')) 
  subtype_codes <- df %>% 
    select(eid, starts_with(c('X40006', 'X40013')))
  # remove rows with all NAs from each of the two data frames
  subtype_codes <- subtype_codes %>%  
    filter(rowSums(is.na(select(., -eid))) != ncol(.) - 1) %>% arrange(eid)
  subtype_dates <- subtype_dates %>%
    filter(rowSums(is.na(select(., -eid))) != ncol(.) - 1) %>%
    # remove the participants with dates but without codes (n=4)
    filter(eid %in% subtype_codes$eid) %>% arrange(eid)
  # replace irrelevant codes with the placeholder 'DUMMY'
  subtype_codes <- subtype_codes %>%
    mutate(across(starts_with('X40013'), ~ replace_not_in_icd9(.x, icd9) )) %>%
    mutate(across(starts_with('X40006'), ~ replace_not_in_icd10(.x, icd10) ))
  # move NAs to the right
  subtype_codes <- as.data.frame(t(apply(subtype_codes, 1, move_nas_to_right)), 
                                 stringsAsFactors = FALSE)
  # remove empty columns (all NAs)
  subtype_codes <- Filter(function(x)!all(is.na(x)), subtype_codes)
  # retain only dates that refer to the relevant diagnosis
  subtype_dates[is.na(subtype_codes) | subtype_codes == 'DUMMY'] <- NA
  # remove rows with all NAs
  subtype_dates <- subtype_dates %>%
    filter(rowSums(is.na(select(., -eid))) != ncol(.) - 1) %>%
    # calculate the earliest date and select relevant columns
    mutate(subtype_date = reduce(across(starts_with('X40005')), pmin, na.rm = TRUE)) %>%
    select(eid, subtype_date)
  colnames(subtype_dates) <- c('id', paste0('cancer_', subtype, '_date'))
  return(subtype_dates)
}


# function to calculate statistical mode
Mode <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



# import and bind the pseudoscales
# This code bind all pseudo-scales into a single data frame and exports as .RDS. 
# It also exports a list with the n's for the individual scales
# The argument 'version' determines the type of simulation that was performed.
combine_scales <- function(version, j_max = 20,  years){
  library(tidyverse)
  scale_names <- c()
  if (version == 'within_achb' | version == 'within_all'){
    scales_achb <- read.csv('output_files/scales_summary.csv')
    for (i in seq(1: nrow(scales_achb))){
      scale_names[length(scale_names) + 1] <- 
        paste0('score_', scales_achb$scale[i], '.csv')
    }
  } else if (version == 'across_achb' | version == 'across_all'){
    for (i in seq(0, j_max-1)){
      scale_names[length(scale_names) + 1] <- 
        paste0('pseudo_scale_', as.character(i), '.csv')
    }
  }
  # add the registration periods (merge by id and year) and remove years for which the year_present == 0
  scales <- read.csv(paste0('output_files/', version, '/', scale_names[1], sep = '|'))
  data_period <- read.csv('output_files/data_period_long.csv')
  colnames(data_period)[colnames(data_period) == 'eid'] <- 'id'
  data_period <- filter(data_period, year_present != 0 & year >= years[1] & 
                          year <= years[length(years)]) # remove irrelevant years
  data_period$year_present <- NULL
  # by merging on data_period, we will include some previously excluded id's; 
  # this is desired, as those id's were previously removed when pseudo_scales() was run
  # (the id's in question were id's that were prescribed only drugs with invalid 
  # administration routes on given years)
  # some id's are also going to get removed from scales; this are the id's that 
  # do not have more than one year of continuous EHR; since we're calculating
  # annual drug burden, we consider these invalid
  scales <- merge(scales, data_period, by = c('id', 'year'), all.y = TRUE)
  rm(data_period)
  # merge all pseudo-scales
  for (scale_name in scale_names[2:length(scale_names)]){
    print(scale_name)
    scale_new <- read.csv(paste0(paste0('output_files/', version, '/', scale_name)), sep = '|')
    scales <- merge(scales, scale_new, by = c('id', 'year'), all.x = TRUE)
    rm(scale_new)
  }
  # replace NAs with 0s ("data_period" introduced NAs because it contains 
  # participants-years that had been removed from the precursors to "scales" 
  # when invalid administration routes were removed)
  scales[is.na(scales)] <- 0
  # also tidy up the scale sizes
  if (version == 'across_all' | version == 'across_achb'){
    version_abr <- 'across'
  } else if (version == 'within_all' | version == 'within_achb'){
    version_abr <- 'within'
  }
  scale_siz <- read.csv(paste0('output_files/pseudo_scale_size_', version_abr, '.csv'), 
                        header = TRUE)
  colnames(scale_siz) <- c('n')
  scale_siz$scale_name <- colnames(select(scales, -ends_with('_alt') & starts_with('score')))
  
  # export
  write.csv(scale_siz, paste0('output_files/pseudo_scale_size_', version_abr, '.csv'), 
            row.names = FALSE)
  saveRDS(scales, file = paste0('output_files/pseudo_scales_', version, '.Rds'))
  return(scales)
}





### This function takes the participant-ID data frame and transforms it into a 
### participant data frame; the medication burden is averaged.

## The function below takes as inputs: 
#   - `how`: the type of simulation/sampling ('across_all', 'across_achb', 
#     'within_all', 'within_achb')
#   - `file_name`: name of the file with the outcomes to be read; 
#     it should have two columns, where the first one indicates the presence (1) 
#     or absence (0) of the outcome, and the second column indicates the earliest 
#     date of diagnosis
#   - 'outcome_name': 'death', 'dementia', or 'delirium'
#   - `year_range`: range (inclusive) for which the mean of the burden score will 
#     be calculated

## It outputs:
#   - a data frame where each ID is an individual observation and the mean 
#     prescribing burden for that individual is given;
#     it also includes various demographic and other variables that will be used 
#     as covariates
#   - it saves the output to disk in the current directory


prepare_scales <- function(version, outcome_name, year_range){
  library(tidyverse)
  library(lubridate)
  
  start_year <- year_range[1]
  stop_year <- year_range[length(year_range)]
  
  date_outcome <- paste0(outcome_name, '_date')
  year_outcome <- paste0(outcome_name, '_year')
  
  # get file with outcome
  covs <- readRDS('output_files/covariates.Rds')
  outcome <- covs[, c('id', date_outcome, outcome_name)]
  outcome[[year_outcome]] <- as.numeric(format(as.Date(outcome[[date_outcome]]), '%Y'))
  
  # sex, age
  sex_age <- covs %>% 
    select(id, sex, birth_year, birth_date)
  # year first present in prescription sample
  data_period <- read.csv('output_files/data_period_long.csv')
  colnames(data_period)[colnames(data_period) == 'eid'] <- 'id'
  year_first <- filter(data_period, year_present !=0) %>% 
    group_by(id) %>% 
    summarise(year_first = min(year))
  # merge demographic variables and outcomes
  df <- merge(year_first, sex_age, by = 'id')
  df <- merge(df, outcome, by = 'id', all.x = TRUE)
  
  ## choose first year to include for AChB averaging
  # for those for whom sampling started in 2015 or later, choose that year 
  # as first year of achb sampling
  df$year_achb_first[df$year_first >= start_year] <- 
    df$year_first[df$year_first >= start_year]
  # for those for whom sampling started before 2015, choose 2015
  df$year_achb_first[df$year_first < start_year] <- start_year
  
  # subsample with only participants for whom there is data before sampling period 
  # (i.e., remove participants that were first recorded after 2015)
  df <- filter(df, year_first <= stop_year, )
  # calculate the number of sampling years, so we know how much to average over; 
  # also calculate the last year of AChB sampling
  sampling_time <- 
    filter(data_period, year_present !=0 & year >= start_year & year <= stop_year) %>% 
    group_by(id) %>% summarise(sampling_time = n())
  year_achb_last <- 
    filter(data_period, year_present !=0 & year >= start_year & year <= stop_year) %>% 
    group_by(id) %>% summarise(year_achb_last = max(year))
  df <- merge(df, sampling_time, by = 'id')
  df <- merge(df, year_achb_last, by = 'id')
  
  # add inpatient data provider to calculate censoring dates
  data_provider_inpatient <- covs %>%
    select(id, data_provider_inpatient_last)
  df <- merge(df, data_provider_inpatient, by = 'id', all.x = TRUE)
  
  # calculate censoring year
  loss_to_followup <- covs %>% 
    select(c(id, loss_to_followup, death_date))
  if (outcome_name == 'death'){
    df <- merge(df, # if outcome is death, we already have the date
                subset(loss_to_followup, select = -death_date), 
                by = 'id', all.x = TRUE)
    df$censor_date <- pmin(df[[date_outcome]], df$loss_to_followup, na.rm = TRUE)
  } else{
    df <- merge(df, loss_to_followup, by = 'id', all.x = TRUE)
    df$censor_date <- pmin(df[[date_outcome]], df$loss_to_followup, df$death_date, na.rm = TRUE)
  }
  # set the correct censoring date for non-events based on data provider
  df$censor_date[df$data_provider_inpatient_last == 'HES' & is.na(df$censor_date)] <- 
    as.Date('31.10.2022', format = '%d.%m.%Y')
  df$censor_date[df$data_provider_inpatient_last == 'SMR' & is.na(df$censor_date)] <- 
    as.Date('31.8.2022', format = '%d.%m.%Y')
  df$censor_date[df$data_provider_inpatient_last == 'PEDW' & is.na(df$censor_date)] <- 
    as.Date('31.5.2022', format = '%d.%m.%Y')
  
  # remove cases lost to follow-up before or during the sampling period)
  df$censor_year <- as.numeric(format(as.Date(df$censor_date), '%Y'))
  df <- df[(df$censor_year > stop_year) | (is.na(df$censor_year)) ,]
  # remove NA cases
  df <- df[ !is.na(df[[outcome_name]]), ]
  
  # calculate age at the last year of sampling
  df$date_achb_last <- ymd(paste0(df$year_achb_last, '-12-31'))
  df$age <- as.numeric(difftime(df$date_achb_last, 
                                as.Date(df$birth_date, format = '%Y-%m-%d'), 
                                units = 'days'))/365.25
  
  # calculate follow-up
  df[[paste0('follow_up_', outcome_name)]] <- 
    as.numeric(difftime(df$censor_date, df$date_achb_last, units = 'days')/365.25)
  
  # add additional column for competing events
  if (outcome_name != 'death'){
    df$death <- ifelse(!is.na(df$death_date), 1, 0)
    df$status <- 0
    df$status[df$death == 1] <- 2
    df$status[df[[outcome_name]] == 1] <- 1
  } else{
    df$status <- df$death
  }
  
  # add data provider
  data_provider <- read.csv('output_files/meds_de-branded.csv', sep = '|') %>%
    select(id, data_provider, date)
  data_provider$date <- as.Date(data_provider$date, format = '%d/%m/%Y')
  data_provider$year <- as.numeric(format(as.Date(data_provider$date), '%Y'))
  data_provider <- data_provider %>% 
    group_by(id, year) %>% 
    summarise(data_provider = Mode(data_provider))
  
  
  # import scales, select only relevant prescribing period, and merge AChB scales with pseudoscales
  meds <- readRDS(paste0('output_files/pseudo_scales_', version, '.Rds'))
  data.table::setDT(meds) # more efficient as data table
  meds <- meds[year >= start_year & year <= stop_year]
  meds <- as.data.frame(meds) # back to data frame
  achb <- read.csv('output_files/achb_scales.csv', sep = '|')
  achb_alt <- read.csv('output_files/achb_scales_poly.csv', sep = '|')
  colnames(achb)[3:ncol(achb)] <- 
    paste('score_', colnames(achb)[3:ncol(achb)], sep='')
  colnames(achb_alt)[3:ncol(achb_alt)] <- 
    paste('score_', colnames(achb_alt)[3:ncol(achb_alt)], sep='')
  achb <- filter(achb, year >= start_year & year <= stop_year)
  achb_alt <- filter(achb_alt, year >= start_year & year <= stop_year)
  achb <- cbind(achb, subset(achb_alt, select = -c(id, year))); rm(achb_alt)
  meds <- merge(achb, meds, by = c('id', 'year'), all.y = TRUE); rm(achb)
  meds[is.na(meds)] <- 0
  meds <- merge(data_provider, meds, by = c('id', 'year'), all.y = TRUE)
  
  # data providers for each participant; for participants for which there is no 
  # data provider info for the given years, impute data provider from other years
  data_provider_sum <- meds %>% 
    select(id, data_provider) %>% 
    group_by(id) %>% 
    summarise(data_provider = Mode(data_provider, na.rm = TRUE))
  dp_missing <- filter(data_provider_sum, is.na(data_provider)) %>% 
    select(id)
  dp_missing <- merge(dp_missing, data_provider, by = 'id', all.x = TRUE)
  # calculate absolute distance from start and finish of relevant prescribing period
  dp_missing$start_distance <- abs(dp_missing$year - start_year) 
  dp_missing$stop_distance <- abs(dp_missing$year - stop_year)
  dp_missing <- dp_missing %>% 
    transform(min_distance = pmin(start_distance, stop_distance, na.rm = TRUE)) %>% 
    arrange(min_distance)
  dp_missing <- distinct(dp_missing, id, .keep_all = TRUE) # keep the row with the lowest distance
  colnames(dp_missing)[colnames(dp_missing) == 'data_provider'] <- 'data_provider_imputed'
  data_provider_sum <- merge(data_provider_sum, 
                             subset(dp_missing, select = c(id, data_provider_imputed)), 
                             by = 'id', 
                             all = TRUE)
  data_provider_sum$data_provider_imputed[is.na(data_provider_sum$data_provider_imputed)] <- 
    data_provider_sum$data_provider[is.na(data_provider_sum$data_provider_imputed)]
  data_provider_sum$data_provider_imputed[is.na(data_provider_sum$data_provider_imputed)] <- 0
  df <- merge(df, data_provider_sum, by = 'id')
  rm(data_period, data_provider, data_provider_sum, dp_missing, outcome, sex_age, 
     year_achb_last, year_first)
  # this has imputed data providers for 30,638 participants; 1,941 participants 
  # remain without data providers and were assigned dp 0
  
  # calculate the average
  meds$year <- NULL
  meds$data_provider <- NULL
  meds$birth_date <- NULL
  
  data.table::setDT(meds)
  # in cases of multiple years: sum across the years
  meds <- meds[, lapply(.SD, sum, na.rm = TRUE), by = id] 
  meds <- data.frame(meds)
  meds <- merge(df[, c('id', 'sex', 'age', 'data_provider', 
                       'data_provider_imputed', 'sampling_time', outcome_name,
                       'status', 'censor_date',
                       paste0('follow_up_', outcome_name))], 
                meds, 
                by = 'id')
  # in cases of multiple years: average across the years
  meds <- meds %>% 
    mutate(across(starts_with('score'), ~ . / sampling_time))
  
  # add the covariates
  covs <- covs[, c('id', setdiff(colnames(covs), colnames(meds)))]
  meds <- merge(meds, covs, by = 'id', all.x = TRUE)
  rm(covs)
  
  # identify columns to mutate
  columns_to_mutate <- setdiff(grep('_date$', names(meds), value = TRUE), 
                               c('birth_date', 'death_date', date_outcome,
                                 'censor_date'))
  # change first dates of diagnosis to years
  meds[columns_to_mutate] <- lapply(meds[columns_to_mutate], function(col) {
    year(as.Date(col, format = '%d/%m/%Y'))
  })
  # if occurrence after end of sampling, change date to NA
  meds[columns_to_mutate] <- lapply(meds[columns_to_mutate], function(col) {
    if_else(col > stop_year, NA_real_, col)
  })
  # if the date of occurrence is NA, change disease coding to 0, 
  # effectively removing cases after end of sampling
  disease_cols <- str_remove_all(names(meds)[str_detect(names(meds), '_date$')], '_date')
  disease_cols <- disease_cols[!disease_cols %in% c('birth', 'censor', 'death', outcome_name)]
  meds <- meds %>% 
    mutate(across(all_of(disease_cols), 
                  ~if_else(is.na(meds[[paste0(cur_column(), '_date')]]), 0, .)))
  
  # change some variables to factors
  meds <- meds %>% 
    mutate(across(c(id, sex, mood_dis, diabetes, hypertension, hyperlip, psych_dis, 
                    hear_loss_any, cns_vasc, cns_infl, cns_atroph, cns_mov, 
                    cns_demyel, cns_parox, cns_other, cns_tbi, cns_cancer,
                    vision_problem, sleep_dis_any, endocrine_dis, nutr_dis, 
                    metabolic_dis, cerebrovascular, respiratory, hepatic, flu, 
                    heart, dementia, delirium, cancer_colon, cancer_prostate_ovary, 
                    cancer_lung, cancer_breast, cancer_ovary, status,
                    starts_with(c('data_prov', 'education', 'alc_freq', 'smoking', 
                                  'phys_act', 'depressed', 'lonely', 'soc_isol'))), 
                  as.factor))
  
  # export
  saveRDS(meds, file = paste0('output_files/pseudo_scales_summarised_', version, '_', outcome_name, '.Rds'))
  return(meds)
}







###  Calculate the RR, the OR, and the CIs for logistic models predicting binary outcomes with medication burden and the inclusion of covariates

## Takes as inputs:
#   - `version`: the type of simulation/sampling ('across_all', 'across_achb', 
#     'within_all', 'within_achb')
#   - `outcome_name`: 'death', 'dementia', or 'delirium'
#   - `control`: 'basic' (i.e., none) vs. 'full' adjustment
#   - `smote`: whether SMOTE should be run to adjust group imbalance;
#     WARNING: implemented only for logistic regression
#   - `model_type`: 'logistic' vs. 'survival'
#   - `competing_death`: TRUE or FALSE; whether subdistribution hazard for the 
#     competing event of death should also be calculated 
#   - `other_predictors`: list of column names; whether the effects of other 
#     variables should be included in the output;
#     WARNING: implemented only for numerical or binary variables in logistic regression
#   - `file_path`: the path to the folder where the 'pseudo_scales_summarised'
#     file with the medication burden according to each scale is saved
#   - `output_file_name`: complete path and file name of the output which will be saved to disk
#   - `core_number`: the number of cores to dedicate to the computation

## Returns as output:
#   - a data frame with the effect sizes for each medication burden scale, the sample sizes, and the number of drugs included in the scales
#   - also saves the output to the current directory

outcome_effect_parallel <- function(version, 
                                    outcome_name, 
                                    control, 
                                    smote = FALSE,
                                    model_type,
                                    competing_death = FALSE,
                                    other_predictors = c(''),
                                    file_path = getwd(), 
                                    output_file_name,
                                    core_number = 1){
  
  
  # error messages for combinations of arguments that have not been implemented
  if (smote == TRUE & model_type == 'survival') {
    stop('Function not implemented to perform SMOTE for time-to-event modelling.
           Please run logistic model or disable SMOTE.')
  }
  
  if (other_predictors != '' & model_type == 'survival') {
    stop('Function not implemented to estimate time-to-event for multiple variables.
           Please run logistic model or leave `other_predictors` at default value.')
  }
  
  
  library(tidyverse)
  library(caret)
  library(bigmemory)
  
  scales <- readRDS(file = paste0(file_path, '/output_files/pseudo_scales_summarised_', version, '_', outcome_name, '.Rds'))
  
  
  # clean the data for analysis
  output <- prep_for_analysis(df = scales, 
                              model_type = model_type, 
                              outcome_name = outcome_name)
  scales <- output[[1]]
  
  # set aside big matrix with burden scores
  burden_scales <- scales %>%
    select(starts_with('score_'))
  burden_scales_big_mat <- as.big.matrix(burden_scales)
  
  # vector with burden scale names
  burden_scale_colnames <- colnames(burden_scales)
  
  # covariates used for adjustment and SMOTE
  if (outcome_name == 'death'){
    covs_compl <- c('sex', 'age', 'data_provider_imputed', 'education_0', 'deprivation', 
                    'alc_freq_0', 'waist_0', 'smoking_0', 'phys_act', 'diabetes', 
                    'cerebrovascular', 'respiratory', 'hepatic', 'flu', 'heart', 'cancer_colon', 
                    'cancer_prostate_ovary', 'cancer_lung', 'cancer_breast')
  } else if (outcome_name == 'dementia'){
    covs_compl <- c('sex', 'age', 'data_provider_imputed', 'education_0', 'deprivation', 
                    'g_0', 'pollution_pc', 'alc_freq_0', 'waist_0', 'smoking_0', 
                    'phys_act', 'mood_dis', 'diabetes', 'hyperlip', 'hear_loss_any', 
                    'cns_infl', 'cns_atroph', 'cns_mov', 'cns_demyel', 'cns_parox', 
                    'cns_other', 'cns_cancer', 'cns_tbi', 'hypertension', 'heart', 
                    'soc_isol_0', 'cerebrovascular', 'lonely_0', 'depressed_0')
  } else if (outcome_name == 'delirium'){
    covs_compl <- c('sex', 'age', 'data_provider_imputed', 'education_0', 'deprivation', 
                    'g_0', 'alc_freq_0', 'waist_0', 'smoking_0', 'phys_act', 
                    'mood_dis', 'psych_dis',  'sleep_dis_any', 'vision_problem', 
                    'hear_loss_any', 'cns_any', 'endocrine_dis', 'nutr_dis', 
                    'metabolic_dis', 'cerebrovascular', 'soc_isol_0', 'lonely_0')
  }
  
  # save other, non-scale variables separately
  analytic_df <- scales[, c(outcome_name, covs_compl, 'follow_up', 'status')]
  
  # save indices of the burden scale columns within the big matrix
  exposure_colnames <- burden_scale_colnames[!grepl('alt$', burden_scale_colnames)]
  cols_relevant <- which(burden_scale_colnames %in% exposure_colnames)
  
  # create and export the descriptor for the big matrix
  big_desc <- describe(burden_scales_big_mat)
  saveRDS(big_desc, file = 'temp/big_matrix_desc.rds')
  
  # clear unnecessary objects so they are not copied to the parallel workers
  rm(output, scales, burden_scales); gc()    
  
  library(foreach)
  library(doFuture)
  
  set.seed(6)
  
  plan(multisession, workers = core_number)
  
  # access the indices saved above
  outcome <- foreach(i = cols_relevant,
                     .combine = 'rbind', 
                     .options.future = list(packages = c('tidyverse', 
                                                         'caret', 
                                                         'marginaleffects', 
                                                         'smotefamily', 
                                                         'bigmemory'), 
                                            globals = c('analytic_df', 
                                                        'burden_scale_colnames', 
                                                        'cols_relevant',
                                                        'version', 
                                                        'outcome_name', 
                                                        'control', 
                                                        'smote', 
                                                        'model_type', 
                                                        'competing_death',
                                                        'other_predictors', 
                                                        'file_path', 
                                                        'output_file_name'), 
                                            seed = TRUE)) %dofuture% {
                                              
                                              source('helper_functions.R')
                                              
                                              # load the descriptor file and attach the big matrix
                                              big_desc <- readRDS('temp/big_matrix_desc.rds')
                                              burden_scales_big_mat <- attach.big.matrix(big_desc)
                                              
                                              # the name of the burden scale (exposure) in this iteration
                                              scale_name <- burden_scale_colnames[i]
                                              
                                              # the name and index of the scale-specific non-scale 
                                              # polypharmacy variable for this iteration
                                              poly_0 <- paste0(scale_name, '_alt')
                                              i_alt <- which(burden_scale_colnames == poly_0)
                                              
                                              
                                              # create data frame for this iteration by combining
                                              # the outcome of interest, covariates, follow-up info,
                                              # and scale-specific values from the big matrix
                                              analytic_df[[scale_name]] <- burden_scales_big_mat[, i]
                                              analytic_df[[poly_0]] <- burden_scales_big_mat[, i_alt]
                                              
                                              # just keep non-missing rows
                                              analytic_df <- analytic_df[complete.cases(analytic_df), ]
                                              analytic_df$id <- NULL
                                              
                                              # define the covariates depending on the adjustment argument
                                              if (control == 'basic'){
                                                predictor_names <- ''
                                              } else if (control == 'full'){
                                                predictor_names <- setdiff(names(analytic_df), 
                                                                           c(outcome_name, scale_name,
                                                                             'follow_up', 'status'))
                                              }
                                              # collapse into a single string with variables separated by "+"
                                              predictor_str <- paste(predictor_names, collapse = ' + ')
                                              
                                              
                                              # run models
                                              model <- run_models(df = analytic_df, 
                                                                  model_type = model_type, 
                                                                  exposure = scale_name, 
                                                                  covariates = predictor_str,
                                                                  outcome_name = outcome_name,
                                                                  control = control,
                                                                  competing_death = competing_death)
                                              model_output <- model[[1]]
                                              estimates <- model[[2]]
                                              
                                              
                                              
                                              # potentially apply SMOTE
                                              if (smote == TRUE){
                                                analytic_df_smote <- smote_df(df = analytic_df, 
                                                                              outcome_name = outcome_name)
                                                
                                                # run models with SMOTE data
                                                model_smote <- run_models(df = analytic_df_smote,
                                                                          model_type = model_type, 
                                                                          exposure = scale_name, 
                                                                          covariates = predictor_str,
                                                                          outcome_name = outcome_name,
                                                                          control = control,
                                                                          competing_death = competing_death)
                                                model_output_smote <- model_smote[[1]]
                                                estimates_smote <- model_smote[[2]]
                                                
                                                
                                                
                                              } else if (smote == FALSE) {
                                                # if SMOTE is FALSE, set the results to the values
                                                # of the non-SMOTE analysis
                                                analytic_df_smote <- analytic_df[]
                                                model_output_smote <- model_output
                                                estimates_smote <- estimates
                                              }
                                              
                                              # extract effect sizes and create temporary data frame
                                              all_outcomes <- create_temp_df(model_type = model_type,
                                                                             other_predictors = other_predictors,
                                                                             model_output = model_output,
                                                                             model_output_smote = model_output_smote,
                                                                             estimates = estimates,
                                                                             estimates_smote = estimates_smote,
                                                                             predictor_var = scale_name,
                                                                             scale_name = scale_name,
                                                                             scale_df = analytic_df,
                                                                             scale_df_smote = analytic_df_smote)
                                              
                                              # change estimates to NA to reduce size on disk
                                              if (model_type == 'logistic'){
                                                estimates <- c(NA)
                                              }
                                              
                                              return(list(all_outcomes, estimates))
                                            }
  
  # shut down all workers
  plan(sequential); gc()
  
  # extract the two objects
  all_outcomes <- outcome[, 1] %>%
    bind_rows(.)
  
  if (model_type == 'logistic'){
    rm(outcome)
  } else if (model_type == 'survival'){
    survival_data <- outcome[, 2]
    rm(outcome)
  }
  
  # remove nonsensical columns for SMOTE if SMOTE not run
  if (smote == FALSE){
    all_outcomes <- all_outcomes %>% 
      select(-ends_with('_smote'))
  }
  
  # same for other predictors
  if (other_predictors == '' & model_type == 'logistic'){
    all_outcomes <- all_outcomes %>% 
      select(-X_OR)
  }
  
  
  # same for competing risk estimates
  if (model_type == 'survival' & competing_death == FALSE){
    all_outcomes <- all_outcomes %>% 
      select(-starts_with('HR_comp'))
  }
  
  # add the numbers of drugs included in the scales
  if (version == 'across_all' | version == 'across_achb'){
    version_abr <- 'across'
  } else if (version == 'within_all' | version == 'within_achb'){
    version_abr <- 'within'
  }
  
  scale_length <- read.csv(paste0('output_files/pseudo_scale_size_', version_abr, '.csv'))
  scale_length$type <- 'pseudo'
  scale_length_achb <- read.csv('output_files/scale_size.csv')
  scale_length_achb$type <- 'achb'
  scale_length_achb$scale_name <- paste0('score_', scale_length_achb$scale_name)
  scale_length <- rbind(scale_length_achb, scale_length)
  all_outcomes <- merge(scale_length, all_outcomes, by = 'scale_name')
  # if survival modelling, also output the data needed for curve construction
  if (model_type == 'survival'){
    saveRDS(all_outcomes, file = paste0('output_files/', 
                                        output_file_name))
    saveRDS(survival_data, file = paste0('output_files/', 
                                         'survival_data_', 
                                         output_file_name))
    return(list(all_outcomes, survival_data))
  } else{
    saveRDS(all_outcomes, file = paste0('output_files/',
                                        output_file_name))
    return(all_outcomes)
  }
}






# preparation for analysis: removal of outliers, scaling, etc.
prep_for_analysis <- function(df, 
                              model_type, 
                              outcome_name){
  
  # change names and positions of some variables to that code below works well
  df <- df %>% 
    rename(drug_number_unique = score_drug_number_unique, drug_number = score_drug_number) %>% 
    relocate(drug_number_unique, .before = sampling_time)
  
  # remove Scottish records for dementia (due to lack of pollution data)
  if (outcome_name == 'dementia' | length(unique(df$sampling_time)) == 3){
    df$data_provider_imputed <- as.character(df$data_provider_imputed)
    df <- filter(df, data_provider_imputed != '2')
    df$data_provider_imputed <- as.factor(df$data_provider_imputed)
  }
  
  # these columns are the predictors
  cols_relevant_0 <- grep('score_', colnames(df))
  cols_relevant_1 <- grep('_alt', colnames(df), invert = TRUE)
  cols_relevant <- intersect(cols_relevant_0, cols_relevant_1)
  
  # potentially outcome to factor
  if (model_type == 'logistic'){
    df[[outcome_name]] <- as.factor(df[[outcome_name]])
  } else if(model_type == 'survival'){
    df[[outcome_name]] <- as.numeric(df[[outcome_name]])
  }
  # outliers to NAs
  df <- df %>% 
    mutate(across(starts_with('score_'), ~ outliers(.x, var_metric = 4, method = 'SD')))
  
  # normalise all numerical variables
  numeric_columns <- sapply(df, is.numeric)
  numeric_columns[which(colnames(df) %in% c(outcome_name, 'status',
                                            paste0('follow_up_', outcome_name)))] <- FALSE
  df[numeric_columns] <- lapply(df[numeric_columns], function(x) as.vector(scale(x)))
  
  # set follow-up variable to follow-up value of appropriate outcome
  df$follow_up <- df[[paste0('follow_up_', outcome_name)]]
  
  return(list(df, cols_relevant))
}




# create temporary data frame for each parallel worker
create_temp_df <- function(model_type, 
                           other_predictors,
                           model_output, 
                           model_output_smote,
                           estimates, 
                           estimates_smote, 
                           predictor_var,
                           scale_name, 
                           scale_df, 
                           scale_df_smote){
  
  if (model_type == 'logistic'){
    
    coefs_model <- model_output$coefficients
    RR <- estimates[['estimate']]
    OR <- exp(coefs_model[predictor_var, 'Estimate'])
    OR_SE <- coefs_model[predictor_var, 'Std. Error']
    
    coefs_model_smote <- model_output_smote$coefficients
    RR_smote <- estimates_smote[['estimate']]
    OR_smote <- exp(coefs_model_smote[predictor_var, 'Estimate'])
    OR_SE_smote <- coefs_model_smote[predictor_var, 'Std. Error']
    
    other_predictor_results <- setNames(vector('list', 
                                               length(other_predictors)), 
                                        paste0(other_predictors, '_OR'))
    # compute values for other_predictors
    # WARNING: works only for binary or continuous predictors
    if (other_predictors != ''){
      for (predictor in other_predictors){
        predictor_names <- names(coefs_model[, 1])
        name_index <- which(grepl(paste0('^', predictor), predictor_names))
        other_predictor_results[[paste0(predictor, '_OR')]] <- 
          exp(coefs_model[name_index, 1])
        
        if (smote == TRUE){
          other_predictor_results[[paste0(predictor, '_OR_smote')]] <- 
            exp(coefs_model_smote[name_index, 1])
        }
      } 
    } else{
      other_predictor_results$`_OR` <- NA
    }          
  } else if (model_type == 'survival'){
    HR <- model_output[predictor_var, 'exp(coef)']
    HR_SE <- model_output[predictor_var, 'se(coef)']
    HR_comp <- model_output[predictor_var, 'HR_comp']
    HR_comp_low <- model_output[predictor_var, 'HR_comp_low']
    HR_comp_high <- model_output[predictor_var, 'HR_comp_high']
    
    HR_smote <- model_output_smote[predictor_var, 'exp(coef)']
    HR_SE_smote <- model_output_smote[predictor_var, 'se(coef)']
  }
  
  
  if (model_type == 'logistic'){
    # data frame of estimates of other predictors
    other_predictor_df <- as.data.frame(other_predictor_results)
    
    # output data frame of each worker
    temp_outcome <- data.frame(
      scale_name = scale_name,
      RR = RR,
      OR = OR,
      OR_SE = OR_SE,
      n_id = nrow(scale_df),
      RR_smote = RR_smote,
      OR_smote = OR_smote,
      OR_SE_smote = OR_SE_smote,
      n_id_smote = nrow(scale_df_smote)
    )
    
    # combined data frame across workers
    combined_outcome <- cbind(temp_outcome, other_predictor_df)
    
  } else if (model_type == 'survival'){
    combined_outcome <- data.frame(
      scale_name = scale_name,
      HR = HR,
      HR_SE = HR_SE,
      HR_comp = HR_comp,
      HR_comp_low = HR_comp_low,
      HR_comp_high = HR_comp_high,
      n_id = nrow(scale_df),
      HR_smote = HR_smote,
      HR_SE_smote = HR_SE_smote,
      n_id_smote = nrow(scale_df_smote)
    )
  }
  return(combined_outcome)
}






# run SMOTE
smote_df <- function(df, 
                     outcome_name){
  
  # remove variables used for survival modelling
  df_new <- df %>% 
    select(-c(follow_up, status))
  # calculate the necessary factor by which to multiply cases 
  # to reach the same number as non-cases
  dup_size <- floor((sum(df_new[[outcome_name]] == '0') / 
                       sum(df_new[[outcome_name]] == '1'))/3)
  # determine ordinal variables (as opposed to nominal) - required for KNN
  ordinals <- c('alc_freq_0', 'phys_act', 'depressed_0', 'lonely_0')
  ordinals <- ordinals[ordinals %in% names(df_new)] 
  # ordinal variables to numeric (required for KNN)
  df_new[ordinals] <- lapply(df_new[ordinals], as.numeric)
  # one-hot encoding for nominal variables
  dmy <- dummyVars(' ~ .', data = df_new[, -which(names(df_new) == outcome_name)])
  X <- data.frame(predict(dmy, newdata = df_new[, -which(names(df_new) == outcome_name)]))
  target <- df_new[[outcome_name]]
  # apply SMOTE
  data_smote <- smotefamily::SMOTE(X = X, target = target, K = 5, dup_size = dup_size)
  # create new data frame
  df_smote <- data_smote$data
  # remove potential dots inserted into the column names
  colnames(df_smote) <- gsub('\\.', '', colnames(df_smote))
  # first, select the new one-hot encoded variables and round them (because some values are now between 0 and 1)
  categ_vars <- colnames(df_smote %>% 
                           select_if((is.numeric)))
  categ_vars <- categ_vars[!categ_vars %in% colnames(df_new)]
  df_smote <- df_smote %>% 
    mutate(across(starts_with(c(categ_vars, ordinals)), round)) %>%
    mutate(across(starts_with(c(categ_vars, ordinals)), as.factor))
  # second, reverse the hot-one encoding
  # data provider is present everywhere, so we can hard-code it
  df_smote$data_provider_imputed <- '0'
  df_smote$data_provider_imputed[df_smote$data_provider_imputed1 == 1] <- '1'
  df_smote$data_provider_imputed[df_smote$data_provider_imputed2 == 1] <- '2'
  df_smote$data_provider_imputed[df_smote$data_provider_imputed3 == 1] <- '3'
  df_smote$data_provider_imputed[df_smote$data_provider_imputed4 == 1] <- '4'
  categ_vars <- categ_vars[categ_vars != 'data_provider_imputed0' & 
                             categ_vars != 'data_provider_imputed1' & 
                             categ_vars != 'data_provider_imputed2' &
                             categ_vars != 'data_provider_imputed3' & 
                             categ_vars != 'data_provider_imputed4']
  
  # for the binary variables, we need to undo one-hot encoding
  for (old_col in categ_vars){
    new_col <- stringr::str_sub(old_col, end = -2) # remove last character
    df_smote[[new_col]] <- '0'
    df_smote[[new_col]][df_smote[[old_col]] == 1] <- '1'
    df_smote[[old_col]] <- NULL
  }
  df_smote <- subset(df_smote, 
                     select = -c(data_provider_imputed0, 
                                 data_provider_imputed1, 
                                 data_provider_imputed3, 
                                 data_provider_imputed4))
  # rename back outcome variable
  colnames(df_smote)[colnames(df_smote) == 'class'] <- outcome_name
  df_smote[[outcome_name]] <- as.factor(df_smote[[outcome_name]])
  df_smote <- df_smote %>% 
    mutate_if(is.character, as.factor)
  
  return(df_smote)
}







## function to run logistic/survival models and save the estimates
run_models <- function(df, 
                       model_type, 
                       exposure, 
                       covariates, 
                       outcome_name, 
                       control,
                       competing_death){
  
  if (control == 'full'){
    covariates_plus <- paste0(' + ', covariates)
  } else if (control == 'basic'){
    covariates_plus <- covariates
  }
  
  if (model_type == 'logistic'){
    # formula for the models
    formula <- as.formula(paste0(outcome_name, ' ~ ', exposure, covariates_plus))
    
    # apply model
    model <- glm(formula, data = df, family = 'binomial')
    estimates <- marginaleffects::avg_comparisons(model,
                                                  variables = c(exposure),
                                                  vcov = FALSE,
                                                  comparison = 'lnratioavg',
                                                  transform = 'exp')
    model_summary <- summary(model)
    # for Cox regression
  } else if (model_type == 'survival'){
    # dichotomise for the creation of survival curve
    low_cutoff <- quantile(df[[exposure]], 0.25)
    high_cutoff <- quantile(df[[exposure]], 0.75)
    lowest_25 <- df[df[[exposure]] <= low_cutoff, ]
    lowest_25$exposure_binary <- 'low'
    highest_25 <- df[df[[exposure]] >= high_cutoff, ]
    highest_25$exposure_binary <- 'high'
    df_dich <- rbind(lowest_25, highest_25)
    df_dich$exposure_binary <- as.factor(df_dich$exposure_binary)
    
    # basic formula to create survival curve
    formula_simple <- as.formula(paste0('survival::Surv(follow_up,',
                                        outcome_name, ')',
                                        ' ~ exposure_binary'))
    if (control == 'full'){
      # if we are adjusting for covariates, we need weights to later 
      # create adjusted survival curves
      formula_matching <-  as.formula(paste0('exposure_binary ~ ', 
                                             covariates))
      df_matched <- MatchIt::matchit(formula_matching, 
                                     data = df_dich,
                                     method = 'quick',
                                     distance = 'glm')
      # formula to create the survival object
      formula = as.formula(paste0('survival::Surv(follow_up,',
                                  outcome_name, ')',
                                  ' ~ ', exposure, 
                                  covariates_plus))
      estimates <- broom::tidy(survival::survfit(formula_simple, 
                                                 data = df_dich,
                                                 weights = df_matched$weights))
    } else if (control == 'basic'){
      formula <- as.formula(paste0('survival::Surv(follow_up,',
                                   outcome_name, ')',
                                   ' ~ ',
                                   exposure))
      estimates <- broom::tidy(survival::survfit(formula_simple, 
                                                 data = df_dich))
    }
    # run Cox model
    model <- survival::coxph(formula, data = df)
    
    # complementary log-log transformation to later pool survival curves
    estimates <- estimates %>%
      mutate(cloglog = log(-log(1-estimate)))
    estimates$scale_name <- exposure
    model_summary <- summary(model)
    # calculate subdistribution hazard, accounting for competing risk of death
    if (outcome_name != 'death' & competing_death == TRUE){
      crr_formula <- as.formula(paste0('survival::Surv(follow_up, status) ~ ',
                                       exposure, covariates_plus))
      comp_model <- tidycmprsk::crr(crr_formula, data = df)
      HR_comp <- exp(comp_model$tidy$estimate)
      HR_comp_low <- exp(comp_model$tidy$estimate - comp_model$tidy$std.error*1.97)
      HR_comp_high <- exp(comp_model$tidy$estimate + comp_model$tidy$std.error*1.97)
      # create output exception for death
    } else {
      HR_comp <- NA
      HR_comp_low <- NA
      HR_comp_high <- NA
    }
    model_summary <- cbind(model_summary$coefficients, 
                           HR_comp, HR_comp_low, HR_comp_high)
  }
  return(list(model_summary, estimates))
}





## returns the inferential statistics for across-sampling
# - outcome_name: name of the outcome (used as column name in the data frame
#   and as part of the file name)
# - model_type: 'logistic' vs. 'survival'
# - adjustment: 'unadjusted' vs. 'adjusted'
# - smote: can be TRUE only for model_type = 'logistic'
# - competing_death: can be TRUE only for model_type = 'survival'
# - period: which prescription period is used (string)

output_results <- function(outcome_name, model_type, adjustment, smote, 
                           competing_death, period){
  
  # create data frame for outputing results
  output_df <- data.frame(matrix(ncol = 9, nrow = 12))
  colnames(output_df) <- c('outcome', 'effect_type', 'adjusted', 'smote', 'comp_risk',
                           'metric', 'type', 'period', 'estimate')
  
  output_df$metric <- c(rep('range', 3), rep('SI', 3), rep('median', 3),
                        rep('n_sig', 3))
  output_df$type <- rep(c('general', 'anticholinergic', 'ABS'), 4)
  output_df$period <- period
  
  
  # set name of outcome column
  output_df$outcome <- outcome_name
  
  # set adjustment column
  if (adjustment == 'unadjusted'){
    output_df$adjusted <- 0
  } else if (adjustment == 'adjusted'){
    output_df$adjusted <- 1
  }
  
  
  if (model_type == 'logistic'){
    
    # set effect type column
    output_df$effect_type <- 'OR'
    # in logistic regression there is no competing risk adjustment
    output_df$comp_risk <- 0
    
    # read in logistic model output
    outcome_all <- readRDS(paste0('output_files/across_all_', 
                                  outcome_name, 
                                  '_', 
                                  adjustment, 
                                  '.Rds'))
    outcome_achb <- readRDS(paste0('output_files/across_achb_', 
                                   outcome_name, 
                                   '_', 
                                   adjustment, 
                                   '.Rds'))
    
    # determine which effect estimate to use
    if (smote == FALSE){
      effect_type <- 'OR'
      se_type <- 'OR_SE'
      # set SMOTE column
      output_df$smote <- 0
    } else{
      effect_type <- 'OR_smote'
      se_type <- 'OR_SE_smote'
      output_df$smote <- 1
    }
    outcome_all$effect <- outcome_all[[effect_type]]
    outcome_achb$effect <- outcome_achb[[effect_type]]
    outcome_all$se <- outcome_all[[se_type]]
    outcome_achb$se <- outcome_achb[[se_type]]
    
    # calculate CIs
    outcome_all$CI_low <- outcome_all$effect - 1.96*outcome_all$se
    outcome_all$CI_high <- outcome_all$effect + 1.96*outcome_all$se
    outcome_achb$CI_low <- outcome_achb$effect - 1.96*outcome_achb$se
    outcome_achb$CI_high <- outcome_achb$effect + 1.96*outcome_achb$se
    
    # read in survival models
  } else if (model_type == 'survival'){
    output_df$effect_type <- 'HR'
    output_df$smote <- 0
    
    outcome_all <- readRDS(paste0('output_files/across_all_', 
                                  outcome_name, 
                                  '_', 
                                  adjustment, 
                                  '_Cox.Rds'))
    outcome_achb <- readRDS(paste0('output_files/across_achb_', 
                                   outcome_name, 
                                   '_', 
                                   adjustment, 
                                   '_Cox.Rds'))
    
    surv_all <- readRDS(paste0('output_files/survival_data_across_all_', 
                               outcome_name, 
                               '_', 
                               adjustment, 
                               '_Cox.Rds'))
    surv_achb <- readRDS(paste0('output_files/survival_data_across_achb_', 
                                outcome_name, 
                                '_', 
                                adjustment, 
                                '_Cox.Rds'))
    
    # if competing risks set to TRUE, we use cause-specific hazard
    if (competing_death == FALSE | outcome_name == 'death'){
      output_df$comp_risk <- 0
      
      outcome_all$effect <- outcome_all$HR
      outcome_achb$effect <- outcome_achb$HR
      outcome_all$se <- outcome_all$HR_SE
      outcome_achb$se <- outcome_achb$HR_SE
      
      outcome_all$CI_low <- outcome_all$effect - 1.96*outcome_all$se
      outcome_all$CI_high <- outcome_all$effect + 1.96*outcome_all$se
      outcome_achb$CI_low <- outcome_achb$effect - 1.96*outcome_achb$se
      outcome_achb$CI_high <- outcome_achb$effect + 1.96*outcome_achb$se
      
      # if competing risks set to TRUE, we use subdistribution hazard
    } else if (competing_death == TRUE){
      output_df$comp_risk <- 1
      
      outcome_all$effect <- outcome_all$HR_comp
      outcome_achb$effect <- outcome_achb$HR_comp
      
      outcome_all$CI_low <- outcome_all$HR_comp_low
      outcome_all$CI_high <- outcome_all$HR_comp_high
      
      outcome_achb$CI_low <- outcome_achb$HR_comp_low
      outcome_achb$CI_high <- outcome_achb$HR_comp_high
    }
  }
  
  
  # subset
  pseudo_all <- outcome_all %>%
    filter(type == 'pseudo')
  pseudo_achb <- outcome_achb %>%
    filter(type == 'pseudo')
  
  scales_all <- outcome_all %>%
    filter(type == 'achb' & scale_name != 'score_achb_poly')
  scales_achb <- outcome_achb %>%
    filter(type == 'achb' & scale_name != 'score_achb_poly')
  
  # minima and maxima for effects
  output_df$estimate[output_df$metric == 'range' &
                       output_df$type == 'general'] <- 
    paste0(as.character(round(range(pseudo_all$effect)[1], 3)),
           ' - ',
           as.character(round(range(pseudo_all$effect)[2], 3)))
  
  output_df$estimate[output_df$metric == 'range' &
                       output_df$type == 'anticholinergic'] <- 
    paste0(as.character(round(range(pseudo_achb$effect)[1], 3)),
           ' - ',
           as.character(round(range(pseudo_achb$effect)[2], 3)))
  
  output_df$estimate[output_df$metric == 'range' &
                       output_df$type == 'ABS'] <- 
    paste0(as.character(round(range(scales_all$effect)[1], 3)),
           ' - ',
           as.character(round(range(scales_all$effect)[2], 3)))
  
  
  # median of effects
  output_df$estimate[output_df$metric == 'median' &
                       output_df$type == 'general'] <- median(pseudo_all$effect)
  
  output_df$estimate[output_df$metric == 'median' &
                       output_df$type == 'anticholinergic'] <- median(pseudo_achb$effect)
  
  output_df$estimate[output_df$metric == 'median' &
                       output_df$type == 'ABS'] <- median(scales_all$effect)
  
  
  
  # calculate the simulation interval
  effects_all <- sort(pseudo_all$effect)
  effects_achb <- sort(pseudo_achb$effect)
  # calculate the 2.5th and 97.5th percentiles
  low_bound_all <- quantile(effects_all, 0.025)
  up_bound_all <- quantile(effects_all, 0.975)
  low_bound_achb <- quantile(effects_achb, 0.025)
  up_bound_achb <- quantile(effects_achb, 0.975)
  
  output_df$estimate[output_df$metric == 'SI' &
                       output_df$type == 'general'] <- 
    paste0(as.character(round(low_bound_all[[1]], 3)), 
           ' - ',
           as.character(round(up_bound_all[[1]], 3)))
  
  output_df$estimate[output_df$metric == 'SI' &
                       output_df$type == 'anticholinergic'] <- 
    paste0(as.character(round(low_bound_achb[[1]], 3)), 
           ' - ',
           as.character(round(up_bound_achb[[1]], 3)))
  
  output_df$estimate[output_df$metric == 'SI' &
                       output_df$type == 'ABS'] <- NA
  
  
  # calculate 95% CI to see which proportion would be significant
  output_df$estimate[output_df$metric == 'n_sig' &
                       output_df$type == 'general'] <- 
    nrow(filter(pseudo_all, CI_low > 1))
  
  output_df$estimate[output_df$metric == 'n_sig' &
                       output_df$type == 'anticholinergic'] <- 
    nrow(filter(pseudo_achb, CI_low > 1))
  
  output_df$estimate[output_df$metric == 'n_sig' &
                       output_df$type == 'ABS'] <- 
    nrow(filter(scales_all, CI_low > 1))
  
  if (model_type == 'logistic'){
    return(list(output_df, pseudo_all, pseudo_achb, scales_all))
  } else if (model_type == 'survival'){
    return(list(output_df, pseudo_all, pseudo_achb, scales_all, surv_all, surv_achb))
  }
}