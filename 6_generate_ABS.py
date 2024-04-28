# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 18:26:37 2024

@author: jurem

Assigns anticholinergic potency according to each available anticholinergic scale 
to each prescription.

It saves to disk a data frame in which observations are id-year combinations, 
indicating the yearly anticholinergic burden according to each anticholinergic scale.

"""
import pandas as pd

meds = pd.read_csv('output_files/meds_debranded.csv', \
                   header=0, sep="|", dtype = str, encoding = 'cp1252')
# a list of all anticholinergic drugs (according to any scale) that were found in the sample
achb = pd.read_csv ('output_files/aas_in_sample.csv') 
achb = achb.loc[achb['drug'].isin(set(meds.name_generic)), 'drug'].tolist()
# remove unnecessary column 'n'
achb.drop(['n'], axis = 1, inplace = True)
meds['application'] = meds['application'].astype(int).copy()
meds = meds.drop(['data_provider', 'prescription', 'quantity', 'prescription_old', 'name_brand'], \
                 axis = 1)
# data frame with non-anticholinergic polypharmacy that is later used as covariate
meds_alt = meds.copy()[['id', 'date']]
# repeat the below for each anticholinergic scale
for col in achb.columns[1:]:
    print(col)
    # create a drug-score dictionary and a list of the non-0-rated drugs
    scale_dict = achb.loc[achb[col] != 0, ['drug', col]].set_index('drug').to_dict()[col]
    valid_scores = list(set(scale_dict.values()))
    # insert the scores for the scored drugs
    new_col = meds.name_generic.replace(scale_dict)
    new_col = new_col.rename(col)
    # non-scored drugs are given 0 by default  
    new_col.loc[~new_col.isin(valid_scores)] = 0 
    meds = pd.concat([meds, new_col], axis = 1)
    # invalid administration routes are also scored 0
    meds.loc[meds['application'] != 0, col] = 0 

    # non-anticholinergic polypharmacy column - create by inverting the scores 
    # (0s to previously scored drugs, 1s to non-scored 0-value drugs.)
    col_alt = (meds[col].isin(valid_scores) == False).apply(int)
    col_alt = col_alt.rename(col + '_alt')
    meds_alt = pd.concat([meds_alt, col_alt], axis = 1)

# add measures of polypharmacy and transform to id-year format
meds['date'] = pd.to_datetime(meds['date'], format='%d/%m/%Y')
meds['year'] = meds.date.dt.to_period('Y')
polypharmacy = meds[['id', 'year', 'name_generic']]. \
    groupby(['id', 'year'], as_index = True).\
        agg(drug_number = ('name_generic', 'count'), drug_number_unique = ('name_generic', pd.Series.nunique)). \
    unstack(fill_value = 0).stack().reset_index(drop = False, inplace = False)
meds = meds.drop(['date', 'name_generic', 'application'], axis = 1) \
    .groupby(['id','year'], as_index = True).agg('sum').unstack(fill_value = 0).stack()
meds.reset_index(drop = False, inplace = True)
meds.to_csv('output_files/achb_scales.csv',index=False, header=True, sep='|')

# repeat for meds_alt
meds_alt['date'] = pd.to_datetime(meds_alt['date'], format = '%d/%m/%Y', errors = 'coerce')
meds_alt['year'] = meds_alt.date.dt.to_period('Y')
meds_alt = meds_alt.drop(['date'], axis = 1) \
    .groupby(['id','year'], as_index = True).agg('sum').unstack(fill_value = 0).stack()
meds_alt.reset_index(drop = False, inplace = True)
meds_alt = pd.merge(meds_alt, polypharmacy, left_on=['id','year'], right_on = ['id','year'])
meds_alt.to_csv('output_files/achb_scales_poly.csv',index=False, header=True, sep='|')