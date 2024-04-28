# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 17:17:27 2023

@author: jurem



This code takes as inputs: 
    (a) `file` (string): the name of the file of prescriptions from UKBB
    (b) `version` (string): "across" vs. "within" (see manuscript for details)
    (d) `complete_sample` (numeric): a binary value indicating whether sampling 
        should be done on generic names for all common drugs in the sample
        (thus reflecting "general" polypharmacy), or only on anticholinergic 
        drugs in the sample (reflecting "anticholinergic" polypharmacy)
    (d) `j_max` (numeric): the number of data frames created
    (e) `i_max` (numeric): the number of pseudo-scale columns in each data frame
        (thus `j_max` x `i_max` is the total number of created pseudo-scales
    (g) `achb_crit` (numeric): relevant only when `complete_sample == 0`; 
        it sets the minimum number of scales that are required to include 
        an anticholinergic drug for it to be considered for the draw
    (h) `keep_years` (list of two values): the period to which the pseudoscales simulation
       should be limited; the list boundaries are inclusive: e.g., [2015, 2015] will
       keep the year 2015; as opposed to `6_generate_ABS.py` which generates
       an annual burden for each year in the sample, regardless of period of interest,
       it is recommended to restrict "keep_years" to the period of interest
       due to the time it takes to run
    (i) `uniform_min` (numeric): minimum of the uniform distribution from which "drug_n" is drawn
    (x) `uniform_max` (numeric): maximum of the uniform distribution from which "drug_n" is drawn 
    (j) `p_scores` (list of length=5): a list of probabilities with with the scores
        4, 3, 2, 1, and 0.5 are drawn
    (k) `output_folder` (string for path): default path for the exported data frames

The running time is long (several days, depending on the chosen time interval), 
so it prints out the iteration (1 for each to-be-exported pseudo-scale .csv-file) 
and the name of the scale on which it is currently working 
(named 'score_0_0' to 'score_19_49').

It saves to the working directory j_max + 1 files:
    - `j_max .csv` files, each containing the participant id, year, 
      and burden according to i_max different pseudo-scales (1 column for each scale)
    - list `l`, containing the number of drugs in each pseudo-scale, 
      arranged in the order in which the scales were created

It **DOES NOT** return anything.
"""
import numpy as np
import pandas as pd



def pseudo_scales(file = 'meds_de_branded.csv', 
                  version = '', 
                  complete_sample = 1, 
                  j_max = 20, 
                  i_max = 50, 
                  keep_years = [2015, 2015],
                  achb_crit = 1, 
                  uniform_min = 14, 
                  uniform_max = 150, 
                  p_scores = [14/1502, 782/1502, 298/1502, 382/1502, 26/1502], 
                  output_folder = 'output_files/',
                  random_seed = 24):
# (1)

    # import a cleaned prescription file, where an observation is a single prescription
    if isinstance(file, str):
        meds = pd.read_csv(file, header=0, sep="|", dtype = str, encoding = 'cp1252')
    elif isinstance(file, pd.DataFrame):
        meds = file.copy()
    else:
        return(False)
    # keep only drugs with valid administration route
    meds['application'] = meds['application'].astype(int).copy()
    meds = meds.loc[meds['application'] == 0]
    # create column for year
    meds['date'] = pd.to_datetime(meds['date'], format = '%d/%m/%Y')   
    meds['year'] = meds['date'].apply(lambda x: x.year)
    # many years we certainly won't use; so keep only a user-defined interval
    meds = meds[(meds['year'] >= keep_years[0]) & (meds['year'] <= keep_years[1])]
    # create a list of the most common unique generic drug names 
    # (with valid administration routes) in the sample
    drugs_set = list(set(meds.loc[meds['name_generic'].notnull(), 'name_generic']))
    # create data frame to which we will be appending the resulting scales
    subset = meds.loc[meds['name_generic'].notnull(), ['id', 'year', 'name_generic']].copy()
    # create a dummy variable so that unstacking (id as row, year as columns) 
    # does not result in an empty data frame that we cannot stack
    meds['dummy'] = 0
    # transform "meds" so that an observation is an id-year combination for each participant and each year
    meds = meds.drop(['date', 'prescription', 'quantity', 'prescription_old', \
                      'name_generic', 'data_provider', 'application', 'name_brand'], axis=1). \
        groupby(['id','year'], as_index=True).agg('sum').unstack(fill_value=0).stack()
    meds.reset_index(drop = False, inplace = True)
    meds.drop(['dummy'], axis=1, inplace = True)
# (2)
    # reduces "drug_set" to only those drugs that are present in at least 
    # 'achb_crit' number of anticholinergic scales
    if complete_sample == 0:
        drugs_set = pd.read_csv('aas_in_sample.csv')
        # also select only those drugs present in the chosen period
        drugs_set = drugs_set.loc[(drugs_set['n'] >= achb_crit) & \
                                  (drugs_set['drug'].\
                                   isin(set(subset.name_generic))), 'drug'].tolist()
    else:
        pass
# (3)

    backup_subset = subset.copy()
    l = []
    random_seed
    rng = np.random.default_rng(random_seed)
    states_drug_n = []
    states_drug_scores = []
    states_drugs_l = []
    if version == 'across':
        if uniform_max > len(drugs_set):
           uniform_max = len(drugs_set)
           print('\n' + 'WARNING: maximum scale size exceeds size of sample drawn upon. \
                 Reducing the maximum scale size to ' + str(uniform_max) + '.' + '\n')
        else: pass
        for j in range(j_max):
            subset = backup_subset.copy()
            meds = meds[['id', 'year']]
            title = 'pseudo_scale_' + str(j) + '.csv'
            print('Iteration ' + str(j))
            for i in range(i_max):
                col_name = 'score' + '_' + str(j) + '_' + str(i)
                # "drug_n" is chosen from a uniform distribution ranging from 14 to 174, 
                # as these represent the minimum and maximum for the n of previous scales for our sample
                drug_n = round(rng.uniform(low=uniform_min, high=uniform_max))
                states_drug_n.append(rng.bit_generator.state)
                
                # each drug gets assigned a burden score from 0.5 to 4, with 
                #probabilities that correspond to the frequencies of those scores 
                # in anticholinergic scales
                drug_scores = rng.choice(a=[0.5,1,2,3,4], size=drug_n, p=p_scores)
                states_drug_scores.append(rng.bit_generator.state)
                
                drugs_l = rng.choice(a=drugs_set, size=drug_n, replace = False)
                states_drugs_l.append(rng.bit_generator.state)
                
                # the new pseudo-scale is used as a dictionary to score the 
                # drugs in the dataframe; those scores get added as a new column 
                # for each pseudo-scale
                drugs_dict = dict(zip(drugs_l, drug_scores))
                              
                new_col = subset.name_generic.replace(drugs_dict)
                new_col = new_col.rename(col_name)
                new_col.loc[~new_col.isin(list(set(drug_scores)))] = 0
                new_col = pd.to_numeric(new_col) # to numeric (currently "object")
                subset = pd.concat([subset, new_col], axis = 1)
                
                # if a drug was not rated by the scale, it gets the score 1 
                # for non-scale burden
                col_alt = (new_col.isin(drug_scores) == False).apply(int) 
                col_alt = col_alt.rename(col_name + '_alt')
                subset = pd.concat([subset, col_alt], axis = 1)                
                
                l.append(drug_n)
            subset = subset.drop(['name_generic'], axis = 1).\
                groupby(['id','year'], as_index=True).agg('sum').\
                    unstack(fill_value=0).stack()
            subset.reset_index(drop = False, inplace = True)
            subset = subset.loc[(subset['year'] >= keep_years[0]) & \
                                (subset['year'] <= keep_years[1])]
            meds = pd.merge(meds, subset,  how='left', left_on=['id','year'], \
                            right_on = ['id','year'])
            meds.fillna(0, inplace = True)    
            meds.to_csv(output_folder + title, index=False, header=True, sep='|')
    elif version == 'within':
        scales = pd.read_csv('scales_summary.csv')
        for j in range(len(scales)):
            subset = backup_subset.copy()
            meds = meds[['id', 'year']]
            title = 'score_' + scales.loc[j, 'scale'] + '.csv'
            drug_n = scales.loc[j, 'count_any']
            # if "drug_n" is greater than the set of all drugs "drugs_set" 
            # (this can happen when sampling from only among anticholinergics 
            # if the criterion for what is considered anticholinergic is set too high)
            if drug_n > len(drugs_set):
                drug_n = len(drugs_set)
                print('\n' + 'WARNING: maximum scale size exceeds size of sample \
                      drawn upon. Reducing the maximum scale size to ' + str(drug_n) + '.' + '\n')
            else: pass
            num_4 = scales.loc[j, 'count_4']
            num_3 = scales.loc[j, 'count_3']
            num_2 = scales.loc[j, 'count_2']
            num_1 = scales.loc[j, 'count_1']
            num_05 = scales.loc[j, 'count_05']
            drug_scores = np.repeat(np.array([4, 3, 2, 1, 0.5]), \
                                    np.array([num_4, num_3, num_2, num_1, num_05])).tolist()
            print(str(j) + ': ' + scales.loc[j, 'scale'])
            for i in range(i_max):
                col_name = 'score' + '_' + scales.loc[j, 'scale'] + '_' + str(i)
                
                drugs_l = rng.choice(a=drugs_set, size=drug_n, replace = False)
                states_drugs_l.append(rng.bit_generator.state)
                
                drugs_dict = dict(zip(drugs_l, drug_scores))
                
                new_col = subset.name_generic.replace(drugs_dict)
                new_col = new_col.rename(col_name)
                new_col.loc[~new_col.isin(list(set(drug_scores)))] = 0
                new_col = pd.to_numeric(new_col)
                subset = pd.concat([subset, new_col], axis = 1)
                
                col_alt = (new_col.isin(drug_scores) == False).apply(int)
                col_alt = col_alt.rename(col_name + '_alt')
                subset = pd.concat([subset, col_alt], axis = 1)
                l.append(drug_n)
            # transform "subset" to id-year format, merge with "meds", tidy, and export  
            subset = subset.drop(['name_generic'], axis = 1).\
                groupby(['id','year'], \
                        as_index=True).agg('sum').unstack(fill_value=0).stack()
            subset.reset_index(drop = False, inplace = True)
            subset = subset.loc[(subset['year'] >= keep_years[0]) & \
                                (subset['year'] <= keep_years[1])]
            meds = pd.merge(meds, subset,  how='left', left_on=['id','year'], \
                            right_on = ['id','year'])
            meds.fillna(0, inplace = True)    
            meds.to_csv(output_folder + title, index=False, header=True, sep='|')                
    else: 
        raise ValueError('Specified version does not exist. Function aborted.')
# (4)
    pd.DataFrame(l).to_csv(output_folder + 'pseudo_scale_size_' + \
                           version + '.csv', index=False, header=True)
    pd.DataFrame(states_drugs_l).to_csv(output_folder + \
                                        'states_drugs_l.csv', index=False, header=True)
    if version == 'across':
        pd.DataFrame(states_drug_n).to_csv(output_folder + \
                                           'states_drug_n.csv', index=False, header=True)
        pd.DataFrame(states_drug_scores).to_csv(output_folder + \
                                                'states_drug_scores.csv', index=False, header=True)
    else:
        pass