# select the period for which the cumulative drug burden is to be calculated

years_used <- c(2015, 2015)
#years_used <- c(2004, 2005, 2006)

# vector with names of all four sampling appraoches
versions_all <- c('across_all', 'across_achb', 'within_all', 'within_achb')

### combine_scales ###
### if the data with the scales was exported in several chunks, combine them into a single file
for (current_version in versions_all){
  combine_scales(version = current_version, j_max = 20,  years = years_used)
}



### prepare_scales ###
### prepare the data frame for analysis (remove cases before time 0, add covariates, etc.)
for (current_version in versions_all){
  death <- prepare_scales(version = current_version, 'death_date.csv', 'death', years_used)
  rm(death); gc()
  dementia <- prepare_scales(version = current_version, 'ADOs_dementia.csv', 'dementia', years_used)
  rm(dementia); gc()
  delirium <- prepare_scales(version = current_version, 'FOs_delirium.csv', 'delirium', years_used)
  rm(delirium); gc()
}



### outcome_effect_parallel ###
### Running logistic regression models for each scale to predict the outcome using the burden score

### WARNING: this is very memory-hungry when run on several cores; 
### the `core_number` used below was the maximum number that could be stably run
### on a machine with 96 Gb RAM
death_all <- outcome_effect_parallel(version = 'across_all', outcome_name = 'death',
                                 file_path = getwd(), 
                                 output_file_name = 'across_all_death.Rds', 
                                 core_number = 10)
dementia_all <- outcome_effect_parallel(version = 'across_all', outcome_name = 'dementia',
                                 file_path = getwd(), 
                                 output_file_name = 'across_all_dementia.Rds', 
                                 core_number = 10)
delirium_all <- outcome_effect_parallel(version = 'across_all', outcome_name = 'delirium',
                                 file_path = getwd(), 
                                 output_file_name = 'across_all_delirium.Rds', 
                                 core_number = 10)
rm(list = ls())





## modelling the effects of scales sampled "across achb"; i.e., by sampling from anticholinergic drugs prescribed in the period of interest

death_achb <- outcome_effect_parallel(version = 'across_achb', outcome_name = 'death',
                                 file_path = getwd(),
                                 output_file_name = 'across_achb_death.Rds', 
                                 core_number = 10)
dementia_achb <- outcome_effect_parallel(version = 'across_achb', outcome_name = 'dementia',
                                    file_path = getwd(),
                                    output_file_name = 'across_achb_dementia.Rds', 
                                    core_number = 10)
delirium_achb <- outcome_effect_parallel(version = 'across_achb', outcome_name = 'delirium',
                                    file_path = getwd(),
                                    output_file_name = 'across_achb_delirium.Rds', 
                                    core_number = 10)
rm(list = ls())




## modelling the effects of scales sampled "within all"; i.e., by sampling from all drugs prescribed in the period of interest separately for each scale
death_all <- outcome_effect_parallel(version = 'within_all', outcome_name = 'death',
                                     file_path = getwd(),
                                     output_file_name = 'within_all_death.Rds', 
                                     core_number = 5)
dementia_all <- outcome_effect_parallel(version = 'within_all', outcome_name = 'dementia',
                                        file_path = getwd(),
                                        output_file_name = 'within_all_dementia.Rds', 
                                        core_number = 5)
delirium_all <- outcome_effect_parallel(version = 'within_all', outcome_name = 'delirium',
                                        file_path = getwd(),
                                        output_file_name = 'within_all_delirium.Rds', 
                                        core_number = 4)
rm(list = ls())




## modelling the effects of scales sampled "within achb"; i.e., by sampling from anticholinergic drugs prescribed in the period of interest separately for each scale
death_achb <- outcome_effect_parallel(version = 'within_achb', outcome_name = 'death',
                                     file_path = getwd(),
                                     output_file_name = 'within_achb_death.Rds', 
                                     core_number = 5)
dementia_achb <- outcome_effect_parallel(version = 'within_achb', outcome_name = 'dementia',
                                        file_path = getwd(),
                                        output_file_name = 'within_achb_dementia.Rds', 
                                        core_number = 5)
delirium_achb <- outcome_effect_parallel(version = 'within_achb', outcome_name = 'delirium',
                                        file_path = getwd(),
                                        output_file_name = 'within_achb_delirium.Rds', 
                                        core_number = 4)
rm(list = ls())