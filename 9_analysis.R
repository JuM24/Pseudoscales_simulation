## This script returns the inferential statistics as reported in the manuscript.
## Some of the values for plotting (e.g., ranges of the axes) may need to be manually
## adjusted as they will differ between the three outcomes.

library(tidyverse)



##### --------------------- #####
#                               #
##### -- Across-sampling -- #####
#                               #
##### --------------------- #####


## Get a table with results

# change arguments as needed
adjustment <- 'unadjusted' # 'unadjusted' or 'adjusted'
model_type <- 'logistic' # 'logistic' or 'survival'
competing_death <- FALSE 
smote <- FALSE
period <- '2015'


df_all_outcomes <- rbind(output_results(outcome_name = 'death',
                                        model_type = model_type,
                                        adjustment = adjustment,
                                        smote = smote,
                                        competing_death = competing_death,
                                        period = period)[[1]],
                         output_results(outcome_name = 'dementia',
                                        model_type = model_type,
                                        adjustment = adjustment,
                                        smote = smote,
                                        competing_death = competing_death,
                                        period = period)[[1]],
                         output_results(outcome_name = 'delirium',
                                        model_type = model_type,
                                        adjustment = adjustment,
                                        smote = smote,
                                        competing_death = competing_death,
                                        period = period)[[1]])
write.csv(df_all_outcomes, 'df_all_outcomes.csv')




## Get the figures

# Compare ORs of all drugs with those of anticholinergics.
# Change the argument values as needed

df_all_outcomes <- output_results(outcome_name = 'dementia',
                                  model_type = 'logistic',
                                  adjustment = 'unadjusted',
                                  smote = TRUE,
                                  competing_death = FALSE,
                                  period = 2015)
pseudo_all <- df_all_outcomes[[2]]
pseudo_achb <- df_all_outcomes[[3]]




# Plot overlap
tiff(paste0('output_files/overlap_', outcome_name, '.tif'), units='in', width=6, height=3, res=300)
ggplot() +
  geom_histogram(data = pseudo_all, aes(x = effect, y = after_stat(density)), 
                 fill = '#377eb8', alpha = 0.4, bins = 50) +
  geom_density(data = pseudo_all, aes(x = effect), colour = '#377eb8', 
               alpha = 0.1, adjust = 3) +
  geom_histogram(data = pseudo_achb, aes(x = effect, y = after_stat(density)), 
                 fill = '#e41a1c', alpha = 0.4, bins = 50) +
  geom_density(data = pseudo_achb, aes(x = effect), colour = '#e41a1c', 
               alpha = 0.1, adjust = 3) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.87, 1.38), 
                     breaks = seq(0.9, 1.35, by = 0.05)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 11), 
                       breaks = seq(0, 10, by = 2)) +
  labs(x=NULL, y=NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(face = 'plain', size = 12, color = 'grey15', angle=0),
        axis.title.x = element_text(face = 'bold', size = 11.5, color = 'grey15'),
        axis.text.y = element_text(face = 'plain', size = 12, color = 'grey15'),
        axis.title.y = element_text(face = 'bold', size = 12.5, color = 'grey15'), legend.position = 'none') 
dev.off()

## % overlap (using kernel density estimation)
# Calculate the KDE using the Sheather-Jones method to determine the smoothing bandwith. 
# Then calculate the Overlapping Coefficient (the integral of the smaller of 
# the two probability density functions over their entire range; 
# it ranges from 0 (no overlap) to 1 (complete overlap); interpretation: common area under the PDFs.
bayestestR::overlap(pseudo_all$effect, pseudo_achb$effect)

# central tendency, spread
median(pseudo_all$effect, na.rm = TRUE); median(pseudo_achb$effect, na.rm = TRUE)
IQR(pseudo_all$effect, na.rm = TRUE); IQR(pseudo_achb$effect, na.rm = TRUE)

# difference test between distributions
wilcox.test(pseudo_all$effect, pseudo_achb$effect, 
            alternative = 'two.sided',
            paired = TRUE)




## Relationship betweeen effect size and scale size
# set the number of bins for the histogram
bins <- 110

# calculate a common set of breaks for both datasets; 
# this ensures that both datasets are binned using the same intervals
combined_min <- min(min(pseudo_all$effect, na.rm = TRUE), min(pseudo_achb$effect, na.rm = TRUE))
combined_max <- max(max(pseudo_all$effect, na.rm = TRUE), max(pseudo_achb$effect, na.rm = TRUE))
bin_width <- (combined_max - combined_min) / bins
breaks <- seq(combined_min, combined_max, by = bin_width)

# bin the data for 'pseudo_all' using the common breaks; 
# calculate the count and average for each bin; 
# determine the left edge of each bin for plotting
binned_data_all <- pseudo_all %>%
  mutate(a_bin = cut(effect, breaks = breaks, include.lowest = TRUE)) %>%
  group_by(a_bin) %>%
  summarise(count = n(), b_avg = mean(n, na.rm = TRUE)) %>%
  mutate(a_bin_left = as.numeric(gsub('.*?([0-9.]+),.*', '\\1', as.character(a_bin))))

binned_data_achb <- pseudo_achb %>%
  mutate(a_bin = cut(effect, breaks = breaks, include.lowest = TRUE)) %>%
  group_by(a_bin) %>%
  summarise(count = n(), b_avg = mean(n, na.rm = TRUE)) %>%
  mutate(a_bin_left = as.numeric(gsub('.*?([0-9.]+),.*', '\\1', as.character(a_bin))))

# create the plot for 'binned_data_all'; 
# use geom_rect() to create the histogram; 
# add a gradient fill based on the average value; 
# add dashed lines to indicate specific effects from 'scales_all'
plot_all <- ggplot(binned_data_all) +
  geom_rect(aes(xmin = a_bin_left, xmax = a_bin_left + bin_width, ymin = 0, 
                ymax = count, fill = b_avg), color = 'grey65', alpha = 0.8) +
  scale_fill_gradient(low = '#ffffcc', high = '#800026') +
  geom_segment(data = scales_all, aes(x = effect, xend = effect, y = 0, 
                                      yend = 40, color = n), 
               linetype = 'dashed', linewidth = 0.5) +
  scale_color_gradient(low = '#ffffcc', high = '#800026', name = 'Variable B') +
  scale_x_continuous(expand = c(0, 0), limits = c(0.87, 1.38), 
                     breaks = seq(0.9, 1.35, by = 0.1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 41), 
                     breaks = seq(0, 40, by = 10)) +
  labs(x=NULL, y=NULL) +
  theme_minimal() +
  theme(legend.position='none', 
        axis.text.x = element_text(face = 'plain', size = 15, color = 'grey15', angle=0),
        axis.text.y = element_text(face = 'plain', size = 15, color = 'grey15'))

# create the plot for 'binned_data_achb'
# use geom_rect() to create the histogram; add a gradient fill based on the average value; add dashed lines to indicate specific effects from 'scales_all'
plot_achb <- ggplot(binned_data_achb) +
  geom_rect(aes(xmin = a_bin_left, xmax = a_bin_left + bin_width, ymin = 0, 
                ymax = count, fill = b_avg), color = 'grey65', alpha = 0.8) +
  scale_fill_gradient(low = '#ffffcc', high = '#800026') +
  geom_segment(data = scales_all, aes(x = effect, xend = effect, y = 0,
                                      yend = 40, color = n), linetype = 'dashed', 
               linewidth = 0.5) +
  scale_color_gradient(low = '#ffffcc', high = '#800026', name = 'Variable B') +
  scale_x_continuous(expand = c(0, 0), limits = c(0.87, 1.38), 
                     breaks = seq(0.9, 1.35, by = 0.1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 41), 
                     breaks = seq(0, 40, by = 10)) +
  labs(x=NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position='none', 
        axis.text.x = element_text(face = 'plain', size = 15, color = 'grey15', angle=0),
        axis.text.y = element_text(face = 'plain', size = 15, color = 'grey15'))
# arrange the two plots vertically
gridExtra::grid.arrange(plot_all, plot_achb, nrow=1)

tiff(paste0('output_files/heatmap_', outcome_name, '.tif'), units='in', width=12, height=2.5, res=300)
gridExtra::grid.arrange(plot_all, plot_achb, nrow=1)
dev.off()

# r between OR and n (all pseudo scales)
cor.test(pseudo_all$effect, pseudo_all$n)$estimate
cor.test(pseudo_all$effect, pseudo_all$n)$conf.int
plot(pseudo_all$effect, pseudo_all$n)
# r between OR and n (achb pseudo scales)
cor.test(pseudo_achb$effect, pseudo_achb$n)$estimate
cor.test(pseudo_achb$effect, pseudo_achb$n)$conf.int
plot(pseudo_achb$effect, pseudo_achb$n)
# r between OR and n (anticholinergic scales)
cor.test(scales_all$effect, scales_achb$n)$estimate
cor.test(scales_all$effect, scales_all$n)$conf.int
plot(scales_all$effect, scales_all$n)

# location/quantile in distribution of non-anticholinergic ORs for each anticholinergic scale
cdf_cur <- 0.5 # the CDF we're interested in
achb_scales_all <- filter(outcome_all, type == 'achb')
achb_scales_all <- achb_scales_all %>% 
  mutate(prop_less_than_effect = sapply(effect, ecdf(filter(outcome_all, type == 'pseudo')$effect)))
dim(filter(achb_scales_all, prop_less_than_effect > cdf_cur))
# median proportion of ABS with stronger effects than the pseudoscales
median(achb_scales_all$prop_less_than_effect)
# correlation between "performance" and size
cor.test(achb_scales_all$n, achb_scales_all$prop_less_than_effect) 

achb_scales_achb <- filter(outcome_achb, type == 'achb')
achb_scales_achb <- achb_scales_achb %>% 
  mutate(prop_less_than_effect = sapply(effect, ecdf(filter(outcome_achb, type == 'pseudo')$effect)))
dim(filter(achb_scales_achb, prop_less_than_effect > cdf_cur))
median(achb_scales_achb$prop_less_than_effect)
cor.test(achb_scales_achb$n, achb_scales_achb$prop_less_than_effect)







##### --------------------- #####
#                               #
##### --- Time-to-event --- #####
#                               #
##### --------------------- #####


### For survival modelling, create survival curves


surv_all <- df_all_outcomes[[5]]
surv_achb <- df_all_outcomes[[6]]

# set aside existing anticholinergic scales
scale_names <- read.csv('output_files/scale_size.csv')
scale_names <- paste0('score_', scale_names$scale_name)
scale_indices <- c()
for (i in 1:length(surv_all)){
  if(unique(surv_all[[i]]$scale_name) %in% scale_names){
    scale_indices <- c(scale_indices, i)
  }
  if (unique(surv_all[[i]]$scale_name) == 'score_achb_poly'){
    unused_index <- i
  }
}
surv_scales <- surv_all[setdiff(scale_indices, unused_index)]
surv_all <- surv_all[-scale_indices]
surv_achb <- surv_achb[-scale_indices]

# the `estimate` column is the survival probability S(t): the probability
# to remain without the event at time t

# combine "achb" and "all" curves into a single data frame with
# 4 strata so you can easily plot (only if `time` column match; if not,
# try rounding the column and then merging)
surv_all <- bind_rows(surv_all, .id = 'scale_id') %>%
  filter(estimate != 0 & estimate != 1)
surv_achb <- bind_rows(surv_achb, .id = 'scale_id') %>%
  filter(estimate != 0 & estimate != 1)


M <- 1000 # number of scales

surv_all <- surv_all %>%
  group_by(time, strata) %>%
  summarize(
    Q_bar = mean(cloglog, na.rm = TRUE),
    U_bar = mean((1 / (log(1 - estimate)) * (1 - estimate))^2 * std.error^2, na.rm = TRUE),
    B = var(cloglog, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    T = U_bar + (1 + 1 / M) * B,  # 5 = M - adjust accordingly
    CI_lower = Q_bar - 1.96 * sqrt(T),
    CI_upper = Q_bar + 1.96 * sqrt(T),
    pooled_surv = 1 - exp(-exp(Q_bar)),
    CI95_lower = 1 - exp(-exp(CI_lower)),
    CI95_upper = 1 - exp(-exp(CI_upper)),
    strata = paste0(strata, '_all')
  )

surv_achb <- surv_achb %>%
  group_by(time, strata) %>%
  summarize(
    Q_bar = mean(cloglog, na.rm = TRUE),
    U_bar = mean((1 / (log(1 - estimate)) * (1 - estimate))^2 * std.error^2, na.rm = TRUE),
    B = var(cloglog, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    T = U_bar + (1 + 1 / M) * B,  # 5 = M - adjust accordingly
    CI_lower = Q_bar - 1.96 * sqrt(T),
    CI_upper = Q_bar + 1.96 * sqrt(T),
    pooled_surv = 1 - exp(-exp(Q_bar)),
    CI95_lower = 1 - exp(-exp(CI_lower)),
    CI95_upper = 1 - exp(-exp(CI_upper)),
    strata = paste0(strata, '_achb')
  )

surv <- rbind(surv_all, surv_achb)
surv$pooled_fail <- 1 - surv$pooled_surv

# change names of strata
surv$strata[surv$strata == 'exposure_binary=high_all'] <- 'general, high'
surv$strata[surv$strata == 'exposure_binary=high_achb'] <- 'anticholinergic, high'
surv$strata[surv$strata == 'exposure_binary=low_all'] <- 'general, low'
surv$strata[surv$strata == 'exposure_binary=low_achb'] <- 'anticholinergic, low'

custom_colors <- c('anticholinergic, low' = '#FC8D62',  # Orange
                   'general, high' = '#66C2A5',  # Green
                   'general, low' = '#8DA0CB',  # Blue
                   'anticholinergic, high' = '#E78AC3') 

if (adjustment == 'unadjusted' & outcome_name == 'death'){
  y_axis <- scale_y_continuous(expand = c(0, 0), limits = c(0, 0.09), 
                               breaks = seq(0, 0.09, by = 0.01))
} else if (adjustment == 'unadjusted' & outcome_name == 'dementia'){
  y_axis <- scale_y_continuous(expand = c(0, 0), limits = c(0, 0.025), 
                               breaks = seq(0, 0.025, by = 0.005))
} else if (adjustment == 'unadjusted' & outcome_name == 'delirium'){
  y_axis <- scale_y_continuous(expand = c(0, 0), limits = c(0, 0.027), 
                               breaks = seq(0, 0.027, by = 0.005))
} else if (adjustment == 'adjusted' & outcome_name %in% c('dementia', 'delirium')){
  y_axis <- scale_y_continuous(expand = c(0, 0), limits = c(0, 0.015), 
                               breaks = seq(0, 0.015, by = 0.005))
} else if (adjustment == 'adjusted' & outcome_name == 'death'){
  y_axis <- scale_y_continuous(expand = c(0, 0), limits = c(0, 0.05), 
                               breaks = seq(0, 0.05, by = 0.01))
}


curves <- ggplot(surv, aes(x = time, y = pooled_fail, color = strata)) +
  #geom_ribbon(aes(ymin = CI95_lower, ymax = CI95_upper, fill = strata), alpha = 0.2, color = NA) + 
  geom_smooth(se = FALSE, method = 'loess', span = 0.01) +
  #geom_step() +
  scale_color_manual(values = custom_colors) + 
  labs(x = 'Years', y = paste0(stringr::str_to_title(outcome_name), ' proportion')) +
  theme_minimal() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 6.7), breaks = seq(0, 16, by = 1)) +
  y_axis +
  theme(axis.text.x = element_text(face = 'plain', size = 11, color = 'grey15', angle=0),
        axis.title.x = element_text(face = 'bold', size = 12.5, color = 'grey15'),
        axis.text.y = element_text(face = 'plain', size = 11, color = 'grey15'),
        axis.title.y = element_text(face = 'bold', size = 12.5, color = 'grey15'), 
        legend.position = 'none')

tiff(paste0('output_files/', outcome_name, '_curves_', adjustment, '.tif'), 
     units='in', width=8, height=5, res=300)
curves
dev.off()







##### --------------------- #####
#                               #
##### -- Within-sampling -- #####
#                               #
##### --------------------- #####


# effect size for anticholinergic scale within its 'polypharmacy n'
# delirium
library(stringr)
outcome_all <- readRDS(paste0('output_files/within_all_', outcome_name, '_', adjustment, '.Rds')) %>%
  filter(!scale_name %in% c('score_drug_number', 'score_drug_number_unique') & 
           type != 'poly' &
           scale_name != 'score_achb_poly') %>%
  mutate(OR_CI_low = OR - OR_SE*1.96, OR_CI_high = OR + OR_SE*1.96,
         OR_CI_low_smote = OR_smote - OR_SE_smote*1.96, 
         OR_CI_high_smote = OR_smote + OR_SE_smote*1.96)
outcome_achb <- readRDS(paste0('output_files/within_achb_', outcome_name, '_', adjustment, '.Rds')) %>%
  filter(!scale_name %in% c('score_drug_number', 'score_drug_number_unique') & 
           type != 'poly' &
           scale_name != 'score_achb_poly') %>%
  mutate(OR_CI_low = OR - OR_SE*1.96, OR_CI_high = OR + OR_SE*1.96,
         OR_CI_low_smote = OR_smote - OR_SE_smote*1.96, 
         OR_CI_high_smote = OR_smote + OR_SE_smote*1.96)

outcome_all$effect <- outcome_all$OR_smote
outcome_achb$effect <- outcome_achb$OR_smote

outcome_all$ci_low <- outcome_all$OR_CI_low_smote
outcome_achb$ci_low <- outcome_achb$OR_CI_low_smote

outcome_all$ci_high <- outcome_all$OR_CI_high_smote
outcome_achb$ci_high <- outcome_achb$OR_CI_high_smote


# change name of m-ARS so that it doesn't clash when with ARS when detecting strings below
outcome_all$scale_name <- gsub('score_rudolph_sumukadas', 'score_sumukadas', outcome_all$scale_name)
outcome_achb$scale_name <- gsub('score_rudolph_sumukadas', 'score_sumukadas', outcome_achb$scale_name)
  
# create a plot for each scale's within-sampling
scale_names <- c('score_summers', 'score_han', 'score_ancelin', 'score_carnahan', 'score_chew', 'score_cancelli', 'score_rudolph', 
                 'score_ehrt', 'score_sittironnarit', 'score_boustani', 'score_sumukadas', 'score_duran', 
                 'score_hefner', 'score_nguyen', 'score_bishara', 'score_briet', 'score_kiesel', 'score_nery',
                 'score_jun', 'score_kable', 'score_ramos', 'score_rihani', 'score_yamada')

scale_abbs <- c('DRN', 'CrAS', 'ABC', 'ADS', 'AAS', 'CABS', 'ARS', 'AAS-r', 'ALS', 'ACB', 'm-ARS', 'DS', 'DRS', 'DDS', 
                'AEC', 'AIS', 'GABS', 'BAAS', 'KABS', 'mACB', 'CALS', 'ACSBC', 'YS')

outcome <- data.frame(scale_name = scale_names, effect = NA, low_ci = NA, 
                      high_ci = NA, med_OR_all = NA, med_OR_achb= NA, CDF_all = NA, 
                      CDF_achb = NA, n = NA, overlap = NA, dif = NA)

for (i in seq(1, length(scale_names))){
  s <- scale_names[i]
  lbl <- scale_abbs[i]
  pseudo_all <- outcome_all %>%
    filter(str_detect(scale_name, s) & type == 'pseudo')  
  pseudo_achb <- outcome_achb %>%
    filter(str_detect(scale_name, s) & type == 'pseudo')

  scale_value <- filter(outcome_all, scale_name == s & type == 'achb')$effect
  ci_low <- filter(outcome_all, scale_name == s & type == 'achb')$ci_low
  ci_high <- filter(outcome_all, scale_name == s & type == 'achb')$ci_high
  
  outcome$effect[i] <- scale_value
  outcome$low_ci[i] <- ci_low
  outcome$high_ci[i] <- ci_high
  # fill in the stats for each scale in regard to its within-pseudoscales
  outcome$med_OR_all[i] <- median(pseudo_all$effect)
  outcome$med_OR_achb[i] <- median(pseudo_achb$effect)
  outcome$CDF_all[i] <- sum(pseudo_all$effect <= scale_value)/nrow(pseudo_all) # proportion pf scales with higher effect size than general pseudoscales
  outcome$CDF_achb[i] <- sum(pseudo_achb$effect <= scale_value)/nrow(pseudo_achb) # proportion pf scales with higher effect size than anticholinergic pseudoscales
  outcome$n[i] <- mean(pseudo_all$n) # scale size
  outcome$overlap[i] <- bayestestR::overlap(pseudo_all$effect, pseudo_achb$effect) # overlap between the two groups of pseudoscales
  outcome$dif[i] <- mean(pseudo_achb$effect, na.rm = TRUE) - mean(pseudo_all$effect, na.rm = TRUE) # difference between means of the two groups of pseudoscales
  
  scale_plot <-  ggplot() +
    geom_histogram(data = pseudo_all, aes(x = effect, y = after_stat(density)), fill = '#377eb8', alpha = 0.4, bins = 30) +
    geom_histogram(data = pseudo_achb, aes(x = effect, y = after_stat(density)), fill = '#e41a1c', alpha = 0.4, bins = 30) +
    geom_rect(xmin = ci_low, xmax = ci_high, aes(ymin = -Inf, ymax = Inf), 
              fill = 'grey60', alpha = 0.5) + 
    geom_vline(xintercept = scale_value, linewidth = 0.4, linetype = 'dashed') +
    scale_x_continuous(expand = c(0, 0), limits = c(0.87, 1.38), breaks = seq(0.9, 1.35, by = 0.1)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 15), breaks = seq(0, 15, by = 3)) +
    labs(x=NULL, y=NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(face = 'plain', size = 9, color = 'grey15', angle=0),
          #axis.text.y=element_blank(), 
          axis.text.y = element_text(face = 'plain', size = 9, color = 'grey15'), 
          legend.position = 'none') +
    annotate('text', x=0.92, y=13.7, label=lbl, size = 3.5)
  tiff(paste0('output_files/', s, '_', outcome_name, '.tif'), units='in', width=4, height=3, res=300)
  print(scale_plot)
  dev.off()
}
write.csv(outcome, 'output_files/outcome.csv')

# proportion of scales for different CDFs (change 'cdf' value as needed)
cdf <- 0.68
dim(filter(outcome, CDF_all >= cdf))
dim(filter(outcome, CDF_achb >= cdf))

# median proportion of pseudoscales exhibiting lower ORs than existing scales
median(outcome$CDF_all)
median(outcome$CDF_achb)
