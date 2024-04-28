### Extract the required variables from the main UKBB phenotype file. ###

## This code assumes that the UKBB file is a .csv file and
# that its column names have the format:
# XfieldID.assessment0... etc; so for field ID 53 (date of attending assesment centre),
# we would have X53.0.0, X53.1.0, X53.2.0, X53.3.0 for assessments
# 0, 1, 2, and 3, respectively

file_path <- '' # path to the folder where the UKKB phenotype file is located 
file_name <- '' # name of the UKBB phenotype file
output_path <- '' # path to the folder where the GitHub repo content is

setwd(fle_path)

library(tidyverse)

ukb <- read.csv(file_name)

# main raw data
main_vars <- ukb %>% 
  select(c(eid, starts_with(
    # demographic, lifestyle, and environmental
    c('X52.', 'X34.','X31.', 'X6138.', 'X189.', 'X20116', 'X24003.',  'X24004.', 
      'X24006.', 'X24016.', 'X24017.', 'X24018.', 'X6164.', 'X53.', 'X48.', 
      'X1558.', 'X709.', 'X1031.', 'X6160.', 'X2020.', 'X2050.',
      # dementia
      'X42018.',
      # delirium
      'X130846.',
      # cognitive tests
      'X20016.', 'X20018.', 'X20023.', 'X399.', 'X6351.', 'X6373.', 'X23324.', 
      'X4282.', 'X21004.',
      # date of death
      'X40000.',
      # mood disorders
      'X130890.', 'X130892.', 'X130894.', 'X130896.', 'X130898.', 'X130900.', 
      'X130902.', 
      # psychotic disorders
      'X130874.', 'X130876.', 'X130878.', 'X130880.', 'X130882.', 'X130884.', 
      'X130886.', 'X130888.',
      # respiratory disease
      'X131484.', 'X131486.', 'X131488.', 'X131490.', 'X131492.', 'X131494.', 
      'X131496.',
      # cerebrovascular disease
      'X131360.', 'X131362.', 'X131364.', 'X131366.', 'X131368.', 'X131370.', 
      'X131372.', 'X131374.', 'X131376.', 'X131378.',
      # liver disease
      'X131498.', 'X131658.', 'X131660.', 'X131662.', 'X131664.', 'X131666.', 
      'X131668.', 'X131670.',
      # influenza/pneumonia
      'X131438.', 'X131440.', 'X131442.', 'X131444.', 'X131446.', 'X131448.', 
      'X131450.', 'X131452.', 'X131454.', 'X131456.',
      # heart disease
      'X131296.', 'X131298.', 'X131300.', 'X131302.', 'X131304.', 'X131306.',
      # diabetes
      'X130706.', 'X130708.', 'X130710.', 'X130712.', 'X130714.',
      # hypertension and hyperlipeademia
      'X131286.', 'X130814.',
      # sleep or visual impairment
      'X131212', 'X131060', 'X130920',
      # inpatient data
      'X40022.' ,'X41270.', 'X41280.', 'X41271.', 'X41281.',
      # hearing-related codes
      'X2247.', 'X2257.', 'X20019.', 'X20021.', 'X131258.', 'X131260.',
      'X131259.', 'X131261.',
      # CNS disorders
      'X42028.', 'X42006.',  'X131000.', 'X131110.', 'X131114.', 'X130992.', 'X130994.', 
      'X130996.', 'X130998.', 'X131002.', 'X131004.', 'X131006.', 'X131008.', 
      'X131010.', 'X131012.', 'X131014.', 'X131016.', 'X131018.', 'X131020.', 
      'X131022.', 'X131024.', 'X131026.', 'X131028.', 'X131030.', 'X131032', 
      'X131038.', 'X131040.', 'X131042.', 'X131044.', 'X131046.', 'X131048.', 
      'X131050.', 'X131056.', 'X131058.', 'X131100.', 'X131112.', 'X131116.', 
      'X131120.',
      # cancer
      'X40005.', 'X40006.', 'X40013.',
      # endocrinopathy
      'X130690.', 'X130692.', 'X130694.', 'X130696.', 'X130698.', 'X130700.', 
      'X130702.', 'X130704.', 'X130718.', 'X130720.', 'X130722.', 'X130724.', 
      'X130726.', 'X130728.', 'X130730.', 'X130732.', 'X130734.', 'X130736.', 
      'X130738',  'X130742.', 'X130744.', 'X130746.', 'X130748.',
      # metabolic disorders
      'X130798.', 'X130800.', 'X130802.', 'X130806.', 'X130808.', 'X130810.', 
      'X130812.', 'X130816.', 'X130818.', 'X130820.', 'X130822.', 'X130824.', 
      'X130826.', 'X130828.', 'X130830.', 'X130832.',
      # nutritional deficiencies
      'X130750.', 'X130752.', 'X130756.', 'X130758.', 'X130760.', 'X130762.', 
      'X130764.', 'X130766.', 'X130768.', 'X130770.', 'X130772.', 'X130774.', 
      'X130776.', 'X130778.', 'X130780.', 'X130782.', 'X130784.', 'X130786.', 
      'X130788.',
      # air pollution
      'X24003.', 'X24004.', 'X24006.', 'X24016.', 'X24017.', 'X24018.'
      ))))
saveRDS(main_vars, file = paste0(output_path, '/main_vars.Rds'))