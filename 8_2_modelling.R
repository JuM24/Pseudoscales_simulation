years_used <- c(2015, 2015)

  
#years_used <- c(2004, 2005, 2006)


### combine_scales ###

### if the data with the scales was exported in several chunks, combine them into a single file
combine_scales(version = 'across_all', j_max = 20,  years = years_used)
combine_scales(version = 'across_achb', j_max = 20,  years = years_used)
combine_scales(version = 'within_all', j_max = 20,  years = years_used)
combine_scales(version = 'within_achb', j_max = 20,  years = years_used)



### prepare_scales ###

### prepare the data frame for analysis (remove cases before time 0, add covariates, etc.)

## across_all
death <- prepare_scales(version = 'across_all', 'death_date.csv', 'death', years_used)
rm(death); gc()
dementia <- prepare_scales(version = 'across_all', 'ADOs_dementia.csv', 'dementia', years_used)
rm(dementia); gc()
delirium <- prepare_scales(version = 'across_all', 'FOs_delirium.csv', 'delirium', years_used)
rm(delirium); gc()

# across_achb
death <- prepare_scales(version = 'across_achb', 'death_date.csv', 'death', years_used)
rm(death); gc()
dementia <- prepare_scales(version = 'across_achb', 'ADOs_dementia.csv', 'dementia', years_used)
rm(dementia); gc()
delirium <- prepare_scales(version = 'across_achb', 'FOs_delirium.csv', 'delirium', years_used)
rm(delirium); gc()

## within_all
death <- prepare_scales(version = 'within_all', 'death_date.csv', 'death', years_used)
rm(death); gc()
dementia <- prepare_scales(version = 'within_all', 'ADOs_dementia.csv', 'dementia', years_used)
rm(dementia); gc()
delirium <- prepare_scales(version = 'within_all', 'FOs_delirium.csv', 'delirium', years_used)
rm(delirium); gc()

# within_achb
death <- prepare_scales(version = 'within_achb', 'death_date.csv', 'death', years_used)
rm(death); gc()
dementia <- prepare_scales(version = 'within_achb', 'ADOs_dementia.csv', 'dementia', years_used)
rm(dementia); gc()
delirium <- prepare_scales(version = 'within_achb', 'FOs_delirium.csv', 'delirium', years_used)
rm(delirium); gc()






### outcome_effect_parallel ###
### Running logistic regression models for each scale to predict the outcome using the burden score

## modelling the effects of scales sampled "across all"; i.e., by sampling from all drugs prescribed in the period of interest

death_all <- outcome_effect_parallel(version = 'across_all', outcome_name = 'death', control = 'full', smote = TRUE, 
                                 other_predictors = c('age', 'sex'), 
                                 file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_across_all_death.Rds', 
                                 output_file_name = 'across_all_death_full.Rds', 
                                 core_number = 10)
dementia_all <- outcome_effect_parallel(version = 'across_all', outcome_name = 'dementia', control = 'full', smote = TRUE, 
                                 other_predictors = c('age', 'sex'), 
                                 file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_across_all_dementia.Rds', 
                                 output_file_name = 'across_all_dementia_full.Rds', 
                                 core_number = 10)
delirium_all <- outcome_effect_parallel(version = 'across_all', outcome_name = 'delirium', control = 'full', smote = TRUE, 
                                 other_predictors = c('age', 'sex'), 
                                 file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_across_all_delirium.Rds', 
                                 output_file_name = 'across_all_delirium_full.Rds', 
                                 core_number = 10)
rm(list = ls())





## modelling the effects of scales sampled "across achb"; i.e., by sampling from anticholinergic drugs prescribed in the period of interest

death_achb <- outcome_effect_parallel(version = 'across_achb', outcome_name = 'death', control = 'full', smote = TRUE, 
                                 other_predictors = c('age', 'sex'), 
                                 file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_across_achb_death.Rds', 
                                 output_file_name = 'across_achb_death_full.Rds', 
                                 core_number = 10)
dementia_achb <- outcome_effect_parallel(version = 'across_achb', outcome_name = 'dementia', control = 'full', smote = TRUE, 
                                    other_predictors = c('age', 'sex'), 
                                    file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_across_achb_dementia.Rds', 
                                    output_file_name = 'across_achb_dementia_full.Rds', 
                                    core_number = 10)
delirium_achb <- outcome_effect_parallel(version = 'across_achb', outcome_name = 'delirium', control = 'full', smote = TRUE, 
                                    other_predictors = c('age', 'sex'), 
                                    file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_across_achb_delirium.Rds', 
                                    output_file_name = 'across_achb_delirium_full.Rds', 
                                    core_number = 10)
rm(list = ls())




## modelling the effects of scales sampled "within all"; i.e., by sampling from all drugs prescribed in the period of interest separately for each scale
death_all <- outcome_effect_parallel(version = 'within_all', outcome_name = 'death', control = 'full', smote = TRUE, 
                                     other_predictors = c('age', 'sex'), 
                                     file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_within_all_death.Rds', 
                                     output_file_name = 'within_all_death_full.Rds', 
                                     core_number = 5)
dementia_all <- outcome_effect_parallel(version = 'within_all', outcome_name = 'dementia', control = 'full', smote = TRUE, 
                                        other_predictors = c('age', 'sex'), 
                                        file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_within_all_dementia.Rds', 
                                        output_file_name = 'within_all_dementia_full.Rds', 
                                        core_number = 5)
delirium_all <- outcome_effect_parallel(version = 'within_all', outcome_name = 'delirium', control = 'full', smote = TRUE, 
                                        other_predictors = c('age', 'sex'), 
                                        file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_within_all_delirium.Rds', 
                                        output_file_name = 'within_all_delirium_full.Rds', 
                                        core_number = 4)
rm(list = ls())




## modelling the effects of scales sampled "within achb"; i.e., by sampling from anticholinergic drugs prescribed in the period of interest separately for each scale
death_achb <- outcome_effect_parallel(version = 'within_achb', outcome_name = 'death', control = 'full', smote = TRUE, 
                                     other_predictors = c('age', 'sex'), 
                                     file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_within_achb_death.Rds', 
                                     output_file_name = 'within_achb_death_full.Rds', 
                                     core_number = 5)
dementia_achb <- outcome_effect_parallel(version = 'within_achb', outcome_name = 'dementia', control = 'full', smote = TRUE, 
                                        other_predictors = c('age', 'sex'), 
                                        file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_within_achb_dementia.Rds', 
                                        output_file_name = 'within_achb_dementia_full.Rds', 
                                        core_number = 5)
delirium_achb <- outcome_effect_parallel(version = 'within_achb', outcome_name = 'delirium', control = 'full', smote = TRUE, 
                                        other_predictors = c('age', 'sex'), 
                                        file_path = 'D:/Job/Projects/Pseudo-scales/pseudoscales GitHub/test folder/pseudo_scales_summarised_within_achb_delirium.Rds', 
                                        output_file_name = 'within_achb_delirium_full.Rds', 
                                        core_number = 4)
rm(list = ls())