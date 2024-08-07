## This script returns the inferential statistics as reported in the manuscript.
## Some of the values for plotting (e.g., ranges of the axes) may need to be manually
## adjusted as they will differ between the three outcomes.


setwd('D:/Job')
library(tidyverse)

outcome_name <- 'dementia'


## Predictive value of anticholinergics vs polypharmacy: 

# Compare ORs of all drugs with those of anticholinergics. 

outcome_all <- readRDS(paste0('across_all_', outcome_name, '.Rds'))
outcome_achb <- readRDS(paste0('across_achb_', outcome_name, '.Rds'))

# compare the correlation between the different effect sizes
cur_df <- outcome_all[]
cor.test(cur_df$RR, cur_df$OR); plot(cur_df$RR, cur_df$OR)
cor.test(cur_df$RR_smote, cur_df$OR_smote); plot(cur_df$RR_smote, cur_df$OR_smote)
cor.test(cur_df$RR, cur_df$RR_smote); plot(cur_df$RR, cur_df$RR_smote)
cor.test(cur_df$OR, cur_df$OR_smote); plot(cur_df$OR, cur_df$OR_smote)
cor.test(cur_df$OR, cur_df$RR_smote); plot(cur_df$OR, cur_df$RR_smote) 
cor.test(cur_df$RR, cur_df$OR_smote); plot(cur_df$RR, cur_df$OR_smote)

# choose which effect size to use
outcome_all$effect <- outcome_all$OR_smote
outcome_achb$effect <- outcome_achb$OR_smote

# subset
pseudo_all <- outcome_all %>%
  filter(type == 'pseudo')
pseudo_achb <- outcome_achb %>%
  filter(type == 'pseudo')

scales_all <- outcome_all %>%
  filter(type == 'achb')
scales_achb <- outcome_achb %>%
  filter(type == 'achb')

# minima and maxima for ORs
range(c(outcome_all$effect, outcome_achb$effect))
range(outcome_all$effect)
range(outcome_achb$effect)


# calculate 95% CI to see which proportion would be significant
pseudo_all$CI_low <- pseudo_all$OR_smote - 1.96*pseudo_all$OR_SE_smote
pseudo_all$CI_high <- pseudo_all$OR_smote + 1.96*pseudo_all$OR_SE_smote
pseudo_achb$CI_low <- pseudo_achb$OR_smote - 1.96*pseudo_achb$OR_SE_smote
pseudo_achb$CI_high <- pseudo_achb$OR_smote + 1.96*pseudo_achb$OR_SE_smote
dim(filter(pseudo_all, CI_low > 1))
dim(filter(pseudo_achb, CI_low > 1))

# calculate the simulation interval
effects_all <- sort(pseudo_all$effect)
effects_achb <- sort(pseudo_achb$effect)
# calculate the 2.5th and 97.5th percentiles
low_bound_all <- quantile(effects_all, 0.025)
up_bound_all <- quantile(effects_all, 0.975)
low_bound_achb <- quantile(effects_achb, 0.025)
up_bound_achb <- quantile(effects_achb, 0.975)
low_bound_all; up_bound_all
low_bound_achb; up_bound_achb







# Plot overlap
tiff('overlap_outcome.tif', units='in', width=6, height=3, res=300)
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

# t-test between distributions
t.test(pseudo_all$effect, pseudo_achb$effect, paired = TRUE)

# mean, variance
mean(pseudo_all$effect, na.rm = TRUE); mean(pseudo_achb$effect, na.rm = TRUE)
sd(pseudo_all$effect, na.rm = TRUE); sd(pseudo_achb$effect, na.rm = TRUE)




## Relationship betweeen effect size and scale size
# set the number of bins for the histogram
bins <- 100

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

tiff('heatmap_outcome.tif', units='in', width=12, height=2.5, res=300)
gridExtra::grid.arrange(plot_all, plot_achb, nrow=1)
dev.off()

# r between OR and n (all pseudo scales)
cor.test(pseudo_all$effect, pseudo_all$n)$estimate;
cor.test(pseudo_all$effect, pseudo_all$n)$conf.int
plot(pseudo_all$effect, pseudo_all$n)
# r between OR and n (achb pseudo scales)
cor.test(pseudo_achb$effect, pseudo_achb$n)$estimate
cor.test(pseudo_achb$effect, pseudo_achb$n)$conf.int
plot(pseudo_achb$effect, pseudo_achb$n)
# r between OR and n (anticholinergic scales)
cor.test(outcome_all$effect[outcome_all$type == 'achb'], outcome_all$n[outcome_all$type == 'achb'])$estimate
cor.test(outcome_all$effect[outcome_all$type == 'achb'], outcome_all$n[outcome_all$type == 'achb'])$conf.int
plot(outcome_all$effect[outcome_all$type == 'achb'], outcome_all$n[outcome_all$type == 'achb'])

# location/quantile in distribution of non-anticholinergic ORs for each anticholinergic scale
cdf_cur <- 0.95 # the CDF we're interested in
achb_scales_all <- filter(outcome_all, type == 'achb')
achb_scales_all <- achb_scales_all %>% 
  mutate(prop_less_than_effect = sapply(effect, ecdf(filter(outcome_all, type == 'pseudo')$effect)))
dim(filter(achb_scales_all, prop_less_than_effect > cdf_cur))
# correlation between "performance" and size
cor.test(achb_scales_all$n, achb_scales_all$prop_less_than_effect) 

achb_scales_achb <- filter(outcome_achb, type == 'achb')
achb_scales_achb <- achb_scales_achb %>% 
  mutate(prop_less_than_effect = sapply(effect, ecdf(filter(outcome_achb, type == 'pseudo')$effect)))
dim(filter(achb_scales_achb, prop_less_than_effect > cdf_cur))
cor.test(achb_scales_achb$n, achb_scales_achb$prop_less_than_effect)

## effects of covariates

# SI for sex
quantile(sort(pseudo_all$sex_OR_smote), 0.025)
quantile(sort(pseudo_all$sex_OR_smote), 0.075)
quantile(sort(pseudo_achb$sex_OR_smote), 0.025)
quantile(sort(pseudo_achb$sex_OR_smote), 0.075)

# CI for age
quantile(sort(pseudo_all$age_OR_smote), 0.025)
quantile(sort(pseudo_all$age_OR_smote), 0.075)
quantile(sort(pseudo_achb$age_OR_smote), 0.025)
quantile(sort(pseudo_achb$age_OR_smote), 0.075)








## Within-sampling
# effect size for anticholinergic scale within its 'polypharmacy n'
# delirium
library(stringr)
outcome_all <- readRDS(paste0('within_all_', outcome_name, '.Rds')) %>%
  filter(!scale_name %in% c('score_drug_number', 'score_drug_number_unique') & type != 'poly')
outcome_achb <- readRDS(paste0('within_achb_', outcome_name, '.Rds')) %>%
  filter(!scale_name %in% c('score_drug_number', 'score_drug_number_unique') & type != 'poly')

outcome_all$effect <- outcome_all$OR_smote
outcome_achb$effect <- outcome_achb$OR_smote

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

outcome <- data.frame(scale_name = scale_names, effect = NA, CDF_all = NA, CDF_achb = NA, n = NA, overlap = NA,
                      dif = NA)

for (i in seq(1, length(scale_names))){
  s <- scale_names[i]
  lbl <- scale_abbs[i]
  pseudo_all <- delirium_all %>%
    filter(str_detect(scale_name, s) & type == 'pseudo')
  pseudo_achb <- delirium_achb %>%
    filter(str_detect(scale_name, s) & type == 'pseudo')

  scale_value <- filter(delirium_all, scale_name == s & type == 'achb')$effect
  outcome$effect[i] <- scale_value
  # fill in the stats for each scale in regard to its within-pseudoscales
  outcome$CDF_all[i] <- sum(pseudo_all$effect <= scale_value)/nrow(pseudo_all) # proportion pf scales with higher effect size than general pseudoscales
  outcome$CDF_achb[i] <- sum(pseudo_achb$effect <= scale_value)/nrow(pseudo_achb) # proportion pf scales with higher effect size than anticholinergic pseudoscales
  outcome$n[i] <- mean(pseudo_all$n) # scale size
  outcome$overlap[i] <- bayestestR::overlap(pseudo_all$effect, pseudo_achb$effect) # overlap between the two groups of pseudoscales
  outcome$dif[i] <- mean(pseudo_achb$effect, na.rm = TRUE) - mean(pseudo_all$effect, na.rm = TRUE) # difference between means of the two groups of pseudoscales
  
  scale_plot <-  ggplot() +
    geom_histogram(data = pseudo_all, aes(x = effect, y = after_stat(density)), fill = '#377eb8', alpha = 0.4, bins = 30) +
    geom_histogram(data = pseudo_achb, aes(x = effect, y = after_stat(density)), fill = '#e41a1c', alpha = 0.4, bins = 30) +
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
  tiff(paste0(s, '.tif'), units='in', width=4, height=3, res=300)
  print(scale_plot)
  dev.off()
}
write.csv(outcome, 'outcome.csv')

# proportion of scales for different CDFs
cdf <- 0.95
dim(filter(outcome, CDF_all >= cdf))
dim(filter(outcome, CDF_achb >= cdf))

# median proportion of pseudoscales exhibiting lower ORs than existing scales
median(outcome$CDF_all)
median(outcome$CDF_achb)

# correlation between anticholinergic and general overlap and size of scale
cor.test(outcome$overlap, outcome$n)

tiff(paste0('within_', outcome_name, '.tif'), units='in', width=9, height=2, res=300)
gridExtra::grid.arrange(boustani, carnahan, rudolph, nrow=1)
dev.off()