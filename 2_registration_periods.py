# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 12:46:45 2023
@author: jurem

This transforms the final data frame produced by the algorithm by Darke et al. (2021) 
(where separate columns indicate the beginnings and endingsof period with continuous 
EHR ascertainment) to one where the individual observation is an id-year pair 
and a column indicates whether continuous EHR ascertainment is given for that 
participant for that year. It codes years without ascertainment as not having 
EHR ascertainment.
"""

import pandas as pd

periods = pd.read_csv('output_files/data_period.csv', header = 0)

# set aside the ids with only one period
period_ids = periods.loc[periods['period'] > 1, 'eid'].copy().tolist() # ids with several periods
periods_1 = periods.loc[~periods['eid'].isin(period_ids)].copy() # df with ids with a single period
periods_mul = periods.loc[periods['eid'].isin(period_ids)].copy() # df with ids with several periods

# create year column and identify full years in sample
periods['from'] = pd.to_datetime(periods['from'], format = '%Y-%m-%d', )
periods['to'] = pd.to_datetime(periods['to'], format = '%Y-%m-%d')
periods.loc[:, 'from_y'] = periods.loc[:, 'from'].dt.to_period('Y').astype(str).astype(int) + 1
periods.loc[:, 'to_y'] = periods.loc[:, 'to'].dt.to_period('Y').astype(str).astype(int) - 1

# remove periods that include less than a full calendar year
periods = periods.loc[periods['from_y'] <= periods['to_y']].copy()

# save the remaining periods for error testing later
all_periods = periods.drop(['from', 'to', 'period'], inplace = False, axis = 1).copy()

# create columns for each possible year
year_cols = list(range(min(periods['from_y']), max(periods['from_y']) + 1))
periods = periods.reindex(columns = (['eid', 'period', 'from_y', 'to_y'] + year_cols))

# for each year, check whether it falls within the period of continuous EHR for each participant
for col in year_cols:
    # if the year is withing the prescribing period, add it; otherwise 0
    periods[col] = 0
    bool_list = (col >= periods['from_y']) & (col <= periods['to_y'])
    periods.loc[bool_list, col] = col

# keep only relevant columns and transform to long format (observations are id-year pairs)
periods.index = periods['eid']
periods.drop(['eid', 'period', 'from_y', 'to_y'], inplace = True, axis = 1)
periods = periods.stack().reset_index() 
periods = periods.groupby(['eid', 'level_1'], as_index = False).agg('sum').reset_index()
periods.drop(['index'], inplace = True, axis = 1)
periods.columns = ['eid', 'year', 'year_present']

# remove participants that have opted outs
opt_out = pd.read_csv('participant_opt_out.csv')
opt_out.columns = ['id']
periods['eid'] = periods['eid'].astype(str)
opt_out['id'] = opt_out['id'].astype(str)
periods = periods.loc[~periods['eid'].isin(opt_out.id)] 
periods = periods.reset_index(drop=True)

periods.to_csv('output_files/data_period_long.csv',index=False, header=True)