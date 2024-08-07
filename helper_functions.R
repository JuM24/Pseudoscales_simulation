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
    mutate(across(starts_with("X40013"), ~ replace_not_in_icd9(.x, icd9) )) %>%
    mutate(across(starts_with("X40006"), ~ replace_not_in_icd10(.x, icd10) ))
  # move NAs to the right
  subtype_codes <- as.data.frame(t(apply(subtype_codes, 1, move_nas_to_right)), 
                                 stringsAsFactors = FALSE)
  # remove empty columns (all NAs)
  subtype_codes <- Filter(function(x)!all(is.na(x)), subtype_codes)
  # retain only dates that refer to the relevant diagnosis
  subtype_dates[is.na(subtype_codes) | subtype_codes == "DUMMY"] <- NA
  # remove rows with all NAs
  subtype_dates <- subtype_dates %>%
    filter(rowSums(is.na(select(., -eid))) != ncol(.) - 1) %>%
    # calculate the earliest date and select relevant columns
    mutate(subtype_date = reduce(across(starts_with("X40005")), pmin, na.rm = TRUE)) %>%
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
combine_scales <- function(version = '', j_max = 20,  years = c(2015, 2015)){
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


prepare_scales <- function(version, file_name, outcome_name, year_range = c(2000, 2005)){
  library(tidyverse)
  library(lubridate)
  
  start_year <- year_range[1]
  stop_year <- year_range[length(year_range)]
  
  date_outcome <- paste('date_', outcome_name, sep = '')
  year_outcome <- paste('year_', outcome_name, sep = '')
  subsample_outcome <- paste('subsample_', outcome_name, sep = '')
  
  outcome <- read.csv(file_name)[, 1:2]
  colnames(outcome) <- c('id', date_outcome)
  outcome[outcome == ''] <- NA
  
  # add indicator variable and year variable for outcome
  outcome[[outcome_name]] <- 0
  outcome[[outcome_name]][!is.na(outcome[[date_outcome]])] <- 1
  outcome[[outcome_name]][outcome[[date_outcome]] %in% 
                            c("1900-01-01", "1901-01-01", 
                              "1902-02-02", "1903-03-03", "2037-07-07")] <- NA
  outcome[[date_outcome]][outcome[[date_outcome]] %in% 
                            c("1900-01-01", "1901-01-01", "1902-02-02", 
                              "1903-03-03", "2037-07-07")] <- NA
  
  if (outcome_name %in% c('dementia', 'delirium')){
    outcome[[date_outcome]] <- as.Date(outcome[[date_outcome]], format = '%Y-%m-%d')
  } else if (outcome_name %in% c('death')){
    outcome[[date_outcome]] <- as.Date(outcome[[date_outcome]], format = '%d/%m/%Y')
  }
  
  outcome[[year_outcome]] <- as.numeric(format(as.Date(outcome[[date_outcome]]), "%Y"))
  
  # sex, age
  sex_age <- read.csv('age_sex_formatted.csv', sep = '|') %>% 
    select(id, sex, birth_year, birth_date)
  # year first present in prescription sample
  data_period <- read.csv('data_period_long.csv')
  colnames(data_period)[colnames(data_period) == 'eid'] <- 'id'
  year_first <- filter(data_period, year_present !=0) %>% 
    group_by(id) %>% 
    summarise(year_first = min(year))
  # merge demographic variables and outcomes
  df <- merge(year_first, sex_age, by = 'id')
  df <- merge(df, outcome, by = 'id', all.x = TRUE)
  
  # remove those that don't want to participate in the study anymore
  opt_outs <- read.csv('participant_opt_out.csv')
  df <- filter(df, !id %in% opt_outs$X1005679)
  
  ## choose first year to include for AChB averaging
  # for those for whom sampling started in 2015 or later, choose that year 
  # as first year of achb sampling
  df$year_achb_first[df$year_first >= start_year] <- 
    df$year_first[df$year_first >= start_year]
  # for those for whom sampling started before 2015, choose 2015
  df$year_achb_first[df$year_first < start_year] <- start_year
  
  ## since we'll ideally be averaging over 6 years, we want :
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
  # sampling time 2000 - 2005 (remove cases before or during this period)
  df <- df[(df[[year_outcome]] > stop_year) | (is.na(df[[year_outcome]])) ,]
  # remove NA cases
  df <- df[ !is.na(df[[outcome_name]]), ]
  
  # calculate age at the last year of sampling
  df$date_achb_last <- ymd(paste0(df$year_achb_last, "-12-31"))
  df$age <- as.numeric(difftime(df$date_achb_last, 
                                as.Date(df$birth_date, format = '%Y-%m-%d'), 
                                units = 'days'))/365.25
  
  
  # add data provider
  data_provider <- read.csv('meds_de-branded.csv', sep = '|') %>%
    select(id, data_provider, date)
  data_provider$date <- as.Date(data_provider$date, format = '%d/%m/%Y')
  data_provider$year <- as.numeric(format(as.Date(data_provider$date), "%Y"))
  data_provider <- data_provider %>% 
    group_by(id, year) %>% 
    summarise(data_provider = Mode(data_provider))
  
  
  # import scales, select only relevant prescribing period, and merge AChB scales with pseudoscales
  meds <- readRDS(paste0('pseudo_scales_', version, '.Rds'))
  data.table::setDT(meds) # more efficient as data table
  meds <- meds[year >= start_year & year <= stop_year]
  meds <- as.data.frame(meds) # back to data frame
  achb <- read.csv('achb_scales.csv', sep = '|')
  achb_alt <- read.csv('achb_scales_poly.csv', sep = '|')
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
                       'data_provider_imputed', 'sampling_time', outcome_name)], 
                meds, 
                by = 'id')
  # in cases of multiple years: average across the years
  meds <- meds %>% 
    mutate(across(starts_with('score'), ~ . / sampling_time))
  
  # add the covariates
  cov_dem <- readRDS('covs_demographic.Rds')
  cov_lif <- readRDS('covs_lifestyle.Rds')
  cov_soc <- readRDS('covs_social.Rds')
  cov_bio <- readRDS('covs_biology.Rds')
  cov_dis <- readRDS('covs_disorders.Rds')
  
  meds <- merge(meds, cov_dem %>% select(-sex), by = 'id', all.x = TRUE)
  meds <- merge(meds, cov_lif, by = 'id', all.x = TRUE)
  meds <- merge(meds, cov_soc, by = 'id', all.x = TRUE)
  meds <- merge(meds, cov_bio, by = 'id', all.x = TRUE)
  # remove potentially duplicate columns before merging with disorders
  if (outcome_name %in% c('dementia')) {
    meds <- merge(meds, cov_dis %>% 
                    select(-dementia, -dementia_date, -death, -death_date), 
                  by = 'id', all.x = TRUE)
  } else if (outcome_name %in% c('death')) {
    meds <- merge(meds, cov_dis %>% 
                    select(-death, -death_date), by = 'id', all.x = TRUE)
  } else{
    meds <- merge(meds, cov_dis, by = 'id', all.x = TRUE)
  }
  rm(cov_dem, cov_lif, cov_soc, cov_bio, cov_dis)
  
  # change first dates of diagnosis to years
  meds$birth_date <- NULL
  meds <- meds %>%
    mutate(across(ends_with('_date'), ~year(as.Date(., format = "%d/%m/%Y")))) %>%
    # if occurrence after end of sampling, change date to NA
    mutate(across(ends_with("_date"), ~if_else(. > stop_year, NA_real_, .)))
  # if the date of occurrence is NA, change disease coding to 0, 
  # effectively removing cases after end of sampling
  disease_cols <- str_remove_all(names(meds)[str_detect(names(meds), "_date$")], "_date")
  disease_cols <- disease_cols[disease_cols != "birth"]
  meds <- meds %>% 
    mutate(across(all_of(disease_cols), 
                  ~if_else(is.na(meds[[paste0(cur_column(), "_date")]]), 0, .)))
  
  # change some variables to factors
  meds <- meds %>% 
    mutate(across(c(id, sex, mood_dis, diabetes, hypertension, hyperlip, psych_dis, 
                    hear_loss_any, cns_vasc, cns_infl, cns_atroph, cns_mov, 
                    cns_demyel, cns_parox, cns_other, cns_tbi, cns_cancer,
                    vision_problem, sleep_dis_any, endocrine_dis, nutr_dis, 
                    metabolic_dis, cerebrovascular, respiratory, hepatic, flu, 
                    heart, dementia, cancer_colon, cancer_prostate_ovary, cancer_lung, 
                    cancer_breast, cancer_ovary, 
                    starts_with(c('data_prov', 'education', 'alc_freq', 'smoking', 
                                  'phys_act', 'depressed', 'lonely', 'soc_isol'))), 
                  as.factor))
  
  # export
  saveRDS(meds, file = paste0('pseudo_scales_summarised_', version, '_', outcome_name, '.Rds'))
  return(meds)
}







###  Calculate the RR, the OR, and the CIs for logistic models predicting binary outcomes with medication burden and the inclusion of covariates

## Takes as inputs:
#   - `type`: the type of simulation/sampling ('across_all', 'across_achb', 'within_all', 'within_achb')
#   - `outcome_name`: 'death', 'dementia', or 'delirium'
#   - `control`: 'basic' vs. 'full' adjustment (i.e., only sex and age or all covariates)
#   - `smote`: whether SMOTE should be run to adjust group imbalance
#   - `other_predictors`: whether the effects of other predictors should be included in the output (WORKS ONLY for continuous or binary categorical predictors)
#   - `output_file_name`: complete path and file name of the output which will be saved to disk
#   - `core_number`: the number of cores to dedicate to the loop; all available cores are included by default

## Returns as output:
#   - a data frame with the effect sizes for each medication burden scale, the sample sizes, and the number of drugs included in the scales
#   - also saves the output to the current directory




outcome_effect_parallel <- 
  function(version, outcome_name, control = 'full', smote = TRUE, 
           other_predictors = c(), file_path = file_path, output_file_name, 
           core_number = parallel::detectCores() - 1){
  
  library(tidyverse)
  library(caret)
  
  scales <- readRDS(file = file_path)
  
  # change names and positions of some variables to that code below works well
  scales <- scales %>% 
    rename(drug_number_unique = score_drug_number_unique, drug_number = score_drug_number) %>% 
    relocate(drug_number_unique, .before = sampling_time)
  
  # remove Scottish records for dementia
  #  if (outcome_name == 'dementia'){
  scales$data_provider_imputed <- as.character(scales$data_provider_imputed)
  scales <- filter(scales, data_provider_imputed != '2')
  scales$data_provider_imputed <- as.factor(scales$data_provider_imputed)
  #  }
  
  # these columns are the predictors
  cols_relevant_0 <- grep("score_", colnames(scales))
  cols_relevant_1 <- grep("_alt", colnames(scales), invert = TRUE)
  cols_relevant <- intersect(cols_relevant_0, cols_relevant_1)
  
  # outcome to factor
  scales[[outcome_name]] <- as.factor(scales[[outcome_name]])
  
  # outliers to NAs
  scales <- scales %>% 
    mutate(across(starts_with('score_'), ~outliers(.x, var_metric = 4, method = 'SD')))
  
  
  # normalize all numerical variables except outcome
  numeric_columns <- sapply(scales, is.numeric)
  numeric_columns[which(colnames(scales) == outcome_name)] <- FALSE
  scales[numeric_columns] <- lapply(scales[numeric_columns], function(x) as.vector(scale(x)))
  
  # initialise vector with scale names
  scale_names <- colnames(scales[cols_relevant])
  
  
  
  library(foreach)
  library(doParallel)
  
  set.seed(6)
  registerDoParallel(cores = core_number) # initialise parallel workers
  
  outcome <- foreach(i = seq_along(cols_relevant), 
                     .combine = 'rbind', 
                     .packages = c('tidyverse', 'caret', 'marginaleffects', 'smotefamily')) %dopar% {
    
    scale_name <- scale_names[i]
    
    poly_0 <- paste0(scale_names[i], '_alt') # non-scale polypharmacy
    covs_basic <- c('sex', 'age', 'data_provider_imputed')
    if (outcome_name == 'death'){
      covs_compl <- c('sex', 'age', 'data_provider_imputed', 'education_0', 'deprivation', 
                      'alc_freq_0', 'waist_0', 'smoking_0', 'phys_act_0', 'diabetes', 
                      'cerebrovascular', 'respiratory', 'hepatic', 'flu', 'heart', 'cancer_colon', 
                      'cancer_prostate_ovary', 'cancer_lung', 'cancer_breast', poly_0)
    } else if (outcome_name == 'dementia'){
      covs_compl <- c('sex', 'age', 'data_provider_imputed', 'education_0', 'deprivation', 
                      'g_0', 'pollution_pc', 'alc_freq_0', 'waist_0', 'smoking_0', 
                      'phys_act_0', 'mood_dis', 'diabetes', 'hyperlip', 'hear_loss_any', 
                      'cns_infl', 'cns_atroph', 'cns_mov', 'cns_demyel', 'cns_parox', 
                      'cns_other', 'cns_cancer', 'cns_tbi', 'hypertension', 'heart', 
                      'soc_isol_0', 'cerebrovascular', 'lonely_0', 'depressed_0', poly_0)
    } else if (outcome_name == 'delirium'){
      covs_compl <- c('sex', 'age', 'data_provider_imputed', 'education_0', 'deprivation', 
                      'g_0', 'alc_freq_0', 'waist_0', 'smoking_0', 'phys_act_0', 
                      'mood_dis', 'psych_dis',  'sleep_dis_any', 'vision_problem', 
                      'hear_loss_any', 'cns_any', 'endocrine_dis', 'nutr_dis', 
                      'metabolic_dis', 'cerebrovascular', 'soc_isol_0', 'lonely_0', poly_0)
    }
    
    
    predictor_var <- colnames(scales[cols_relevant[i]])

    # temporary file used to remove observations with missing data
    temp <- scales[, c('id', outcome_name, predictor_var, covs_compl)]
    
    # file with the desired covariates and the scale that is to be assessed
    if (control == 'basic'){
      scale_df <- scales[, c('id', outcome_name, predictor_var, covs_basic)] 
    } else if (control == 'full'){
      scale_df <- scales[, c('id', outcome_name, predictor_var, covs_compl)]
    }
    
    # just keep non-missing rows
    temp <- temp[complete.cases(temp), ]
    scale_df <- filter(scale_df, id %in% temp$id)
    temp$id <- NULL; scale_df$id <- NULL
    
    # get all variable names in the dataframe excluding the outcome variable
    predictor_names <- setdiff(names(scale_df), c(outcome_name, predictor_var))
    # collapse into a single string with variables separated by "+"
    predictor_str <- paste(predictor_names, collapse = " + ")
    # formula for the models
    formula <- as.formula(paste0(outcome_name, " ~ ", predictor_var, '+', predictor_str))
    
    # apply model
    model <- glm(formula, data = scale_df, family = 'binomial')
    estimates_RR <- marginaleffects::avg_comparisons(model,
                                                     variables = c(predictor_var),
                                                     vcov = FALSE,
                                                     comparison = "lnratioavg",
                                                     transform = "exp")
    model <- summary(model)
    
    # potentially apply SMOTE
    if (smote == TRUE){
      df_new <- temp[]
      # calculate the necessary factor by which to multiply cases 
      # to reach the same number as non-cases
      dup_size <- floor((sum(df_new[[outcome_name]] == '0') / 
                           sum(df_new[[outcome_name]] == '1'))/3)
      # determine ordinal variables (as opposed to nominal) - required for KNN
      ordinals <- c('alc_freq_0', 'phys_act_0', 'depressed_0', 'lonely_0', 'sol_isol_0')
      ordinals <- ordinals[ordinals %in% names(df_new)] 
      # ordinal variables to numeric (required for KNN)
      df_new[ordinals] <- lapply(df_new[ordinals], as.numeric)
      # one-hot encoding for nominal variables
      dmy <- dummyVars(" ~ .", data = df_new[, -which(names(df_new) == outcome_name)])
      X <- data.frame(predict(dmy, newdata = df_new[, -which(names(df_new) == outcome_name)]))
      target <- df_new[[outcome_name]]
      # apply SMOTE
      data_smote <- smotefamily::SMOTE(X = X, target = target, K = 5, dup_size = dup_size)
      # create new data frame
      scale_df_smote <- data_smote$data
      # remove potential dots inserted into the column names
      colnames(scale_df_smote) <- gsub('\\.', '', colnames(scale_df_smote))
      # first, select the new one-hot encoded variables and round them (because some values are now between 0 and 1)
      categ_vars <- colnames(scale_df_smote %>% 
                               select_if((is.numeric)))
      categ_vars <- categ_vars[!categ_vars %in% colnames(df_new)]
      scale_df_smote <- scale_df_smote %>% 
        mutate(across(starts_with(c(categ_vars, ordinals)), round)) %>%
        mutate(across(starts_with(c(categ_vars, ordinals)), as.factor))
      # second, reverse the hot-one encoding
      # data provider is present everywhere, so we can hard-code it
      scale_df_smote$data_provider_imputed <- '0'
      scale_df_smote$data_provider_imputed[scale_df_smote$data_provider_imputed1 == 1] <- '1'
      scale_df_smote$data_provider_imputed[scale_df_smote$data_provider_imputed2 == 1] <- '2'
      scale_df_smote$data_provider_imputed[scale_df_smote$data_provider_imputed3 == 1] <- '3'
      scale_df_smote$data_provider_imputed[scale_df_smote$data_provider_imputed4 == 1] <- '4'
      categ_vars <- categ_vars[categ_vars != 'data_provider_imputed0' & 
                                 categ_vars != 'data_provider_imputed1' & 
                                 categ_vars != 'data_provider_imputed2' &
                                 categ_vars != 'data_provider_imputed3' & 
                                 categ_vars != 'data_provider_imputed4']
      # for the binary variables, we need to undo one-hot encoding
      for (old_col in categ_vars){
        new_col <- stringr::str_sub(old_col, end = -2) # remove last character
        scale_df_smote[[new_col]] <- '0'
        scale_df_smote[[new_col]][scale_df_smote[[old_col]] == 1] <- '1'
        scale_df_smote[[old_col]] <- NULL
      }
      scale_df_smote <- subset(scale_df_smote, 
                               select = -c(data_provider_imputed0, 
                                           data_provider_imputed1, 
                                           data_provider_imputed3, 
                                           data_provider_imputed4))
      # rename back outcome variable
      colnames(scale_df_smote)[colnames(scale_df_smote) == 'class'] <- outcome_name
      scale_df_smote[[outcome_name]] <- as.factor(scale_df_smote[[outcome_name]])
      scale_df_smote <- scale_df_smote %>% 
        mutate_if(is.character, as.factor)
      
      # remove variables used for SMOTE if the control is basic
      if (control == 'basic'){
        scale_df_smote <- scale_df_smote[, c(outcome_name, predictor_var, covs_basic)] 
      } else if (control == 'full'){
        scale_df_smote <- scale_df_smote[, c(outcome_name, predictor_var, covs_compl)]
      }
      
      # run the model with SMOTE      
      model_smote <- glm(formula, data = scale_df_smote, family = 'binomial')
      estimates_RR_smote <- marginaleffects::avg_comparisons(model_smote,
                                                             variables = c(predictor_var),
                                                             vcov = FALSE,
                                                             comparison = "lnratioavg",
                                                             transform = "exp")
      model_smote <- summary(model_smote)
      
    } else{
      scale_df_smote <- scale_df[]
      model_smote <- model
      estimates_RR_smote <- estimates_RR
    }
    
    
    tryCatch(
      {
        coefs_model <- model$coef
        RR <- estimates_RR[[3]]
        OR <- exp(coefs_model[2,1])
        OR_SE <- coefs_model[2,2]
        
        coefs_model_smote <- model_smote$coef
        RR_smote <- estimates_RR_smote[[3]]
        OR_smote <- exp(coefs_model_smote[2,1])
        OR_SE_smote <- coefs_model_smote[2,2]
        
        # works only for binary or continuous predictors
        other_predictor_results <- setNames(vector('list', 
                                                   length(other_predictors)), 
                                            paste0(other_predictors, '_OR'))
        # compute values for other_predictors
        if (length(other_predictors) != 0){
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
        }
        
        
      },
      error=function(e) {
        message('An error occurred.')
        print(e)
      }
    )
    
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
    
    return(combined_outcome)
  }
  
  # shut down all workers
  stopImplicitCluster()
  
  # remove nonsensical columns for SMOTE if SMOTE not run
  if (smote == FALSE){
    outcome <- outcome %>% 
      select(-ends_with('_smote'))
  }
  
  
  
  # add the numbers of drugs included in the scales
  if (version == 'across_all' | version == 'across_achb'){
    version_abr <- 'across'
  } else if (version == 'within_all' | version == 'within_achb'){
    version_abr <- 'within'
  }
  
  scale_length <- read.csv(paste0('pseudo_scale_size_', version_abr, '.csv'))
  scale_length$type <- 'pseudo'
  scale_length_achb <- read.csv('scale_size.csv')
  scale_length_achb$type <- 'achb'
  scale_length_achb$scale_name <- paste0('score_', scale_length_achb$scale_name)
  scale_length <- rbind(scale_length_achb, scale_length)
  outcome <- merge(scale_length, outcome, by = 'scale_name')
  saveRDS(outcome, file = output_file_name)
  return(outcome)
}