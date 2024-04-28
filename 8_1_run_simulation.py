# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 17:36:58 2024

@author: jurem
"""

pseudo_scales(file = 'meds_de_branded.csv', 
              version = 'across',
              complete_sample = 1,
              j_max = 20,
              i_max = 50,
              keep_years = [2015, 2015],
              achb_crit = 1,
              uniform_min = 14,
              uniform_max = 150,
              p_scores = [14/1502, 782/1502, 298/1502, 382/1502, 26/1502],
              output_folder = 'output_files/across_all/',
              random_seed = 24)

pseudo_scales(file = 'meds_de_branded.csv',
              version = 'across',
              complete_sample = 0,
              j_max = 20,
              i_max = 50,
              keep_years = [2015, 2015],
              achb_crit = 1,
              uniform_min = 14,
              uniform_max = 150,
              p_scores = [14/1502, 782/1502, 298/1502, 382/1502, 26/1502],
              output_folder = 'output_files/across_achb/',
              random_seed = 24)

pseudo_scales(file = 'meds_de_branded.csv',
              version = 'within',
              complete_sample = 1,
              j_max = 20,
              i_max = 250,
              keep_years = [2015, 2015],
              achb_crit = 1,
              output_folder = 'output_files/within_all/',
              random_seed = 24)

pseudo_scales(file = 'meds_de_branded.csv',
              version = 'within',
              complete_sample = 0,
              j_max = 20,
              i_max = 250,
              keep_years = [2015, 2015],
              achb_crit = 1, 
              output_folder = 'output_files/within_achb/',
              random_seed = 24)