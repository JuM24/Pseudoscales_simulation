# -*- coding: utf-8 -*-
"""
Created on Sat Apr 6 18:25:04 2024

@author: jurem

It saves to disk the input data frame with three new columns: 
    (1) 'name_brand' gives the brand name in the prescription (default is NaN)
    (2) 'name_generic' gives the generic name of the drug in the prescription
        (for all most common drugs in the sample (>=5,000 occurrences); all other rows 
        have NaN's in this column)
    (3) 'application' indicating administration route of the drug as defined above.
(1) contains counts for combinations drugs, (2) counts for non-combinations 
drugs, (3) counts for brand names replaced by generic names. (2) and (3) also 
have columns with values for the simple (non-regex) search terms (i.e., 
simple search of drug name).

It also outputs .csv files with some statistics on the numbers of found &
renamed drugs in the sample.

"""

import pandas as pd

# set to "False" to disable output of statistics about numbers of found/replaced
# drug names; runs faster if set to "False"
produce_stats = True


drug_names = pd.read_csv('drug_names_all.csv', header=0)
drug_names[['alt', 'generic']] = drug_names[['alt', 'generic']].\
    apply(lambda x: x.str.strip())
drug_names.drop_duplicates(inplace = True)

meds = pd.read_csv('output_files/meds_cleaned.csv', header=0, sep="|", dtype = str, 
                   encoding = 'cp1252')

meds_backup = meds.copy()

print('1/5: ' + 'Preparing the sample: replacing some compound names...')
   
for (i,j) in list(zip(['oestradiol',' estrogen','sulphate','beclomethasone',
                       'glycopyronium','bendrofluazide','betametasone',
                       'chlorthalidone','cholecalciferol'],
                      ['estradiol','oestrogen','sulfate','beclometasone',
                       'glycopyrronium','bendroflumethiazide','betamethasone',
                       'chlortalidone','colecalciferol'])):
    meds.loc[:, 'prescription'] = meds.loc[:, 'prescription'].\
        str.replace(i, j, regex = False)   
meds.loc[:, 'prescription'] = meds.loc[:, 'prescription'].\
    str.replace(' +', ' ', regex = True)

print('Replaced terms in ' + str(len(meds.loc[meds['prescription'] != \
                                              meds_backup['prescription']].\
                                     copy())) + ' rows.')

# "combinations.csv" is a long-format data frame that contains 6 columns: #
# (1) the alternative name of a combination drug, (2) the combination, 
# (3-5) 1-3 alternative names for the combination, (6) the number of compounds 
# composing the combination.
combos = pd.read_csv('combinations.csv')
combos[['generic', 'name_1', 'name_2', 'name_3']] = \
    combos[['generic', 'name_1', 'name_2', 'name_3']].apply(lambda x: x.str.strip())
combos.drop_duplicates(inplace = True)

# the dataframe is sorted ascendingly, so that combos with more compounds don't 
# get overwritten by combos with fewer compounds
combos = combos.sort_values(by = 'number', ascending = True)
combos.reset_index(drop=True, inplace=True)

# A  
print('\n' + '2/5: ' + 'Searching for combination drugs... Combination drugs left: ')    
meds.loc[:, 'name_generic'] = None
count = len(combos)
print(count)
for i in range(len(combos)):
    row_terms = []
    for col in ['name_1', 'name_2', 'name_3']:
        if isinstance(combos[col][i], str):
            # splits each combination name into individual compounds and creates 
            # a list identifying the number of searched-for compounds each prescription contains
            term_list = combos[col][i].split('/')
            # the names are in the form drug1/drug2 (in case of 2-combo), 
            # so its unequivocal, hence `regex=False`
            row_terms_col = sum([meds['prescription'].\
                                 str.contains(word, regex = False) for word in term_list]) 
            row_terms_col = list((row_terms_col[row_terms_col == len(term_list)]).index)
            row_terms = list(set(row_terms + row_terms_col))
        else: pass
    # it then assigns the generic combination name to the 'name_generic' column 
    # for all rows that contained all searched-for compounds; this is repeated 
    # for each combination and for each alternative name of the combination
    meds.loc[(list(row_terms)), 'name_generic'] = combos['generic'][i]
    if produce_stats == True:
        if len(row_terms) != 0:
            combos.loc[i, 'n_meds'] = len(meds.loc[row_terms])
            combos.loc[i, 'n_ids'] = len(set(meds.loc[row_terms, 'id']))
        else: pass
    else: pass
    count -= 1
    if (count != 0) & (count % 50 == 0):
        print(count)
print('Identified generic terms for combination drugs in ' + \
      str(len(meds.loc[meds['name_generic'].notnull()])) + ' rows.')
meds.to_csv('meds_combos.csv', index=False, header=True, sep = '|')
combos.to_csv('combos_stats.csv', index=False, header=True)

# B    
meds_small = meds.loc[meds['name_generic'].isnull()].copy()
generics_df = drug_names.loc[drug_names['combination'] == 0].\
    drop_duplicates(subset=['generic']).copy().drop(['alt'], axis = 1, inplace = False)
# create a list of generic names of all most common drugs and then search for 
# them in the prescriptions
generics = list(drug_names.loc[drug_names['combination'] == 0].\
                drop_duplicates(subset=['generic']).copy()['generic'])
count = len(generics)
print('\n' + '3/5: ' + 'Searching for non-combination generic names... Generic names left: ')
print(count)
# assign the names to the "name_generic" column
for g in generics:
    meds_small.loc[(meds_small['prescription'].\
                    str.contains("(?:^|[^a-z])" + g + "(?:$|[^a-z])")) & \
                   (meds_small['name_generic'].isnull()), 'name_generic'] = g
    count -= 1
    if (count != 0) & (count % 50 == 0):
        print(count)
    if produce_stats == True:
        generics_df.loc[generics_df['generic'] == g, 'n_meds'] = \
            len(meds_small.loc[meds_small['name_generic'] == g])
        generics_df.loc[generics_df['generic'] == g, 'n_ids'] = \
            len(set(meds_small.loc[meds_small['name_generic'] == g, 'id']))
        generics_df.loc[generics_df['generic'] == g, 'n_meds_simple'] = \
            len(meds_small.loc[meds_small['prescription'].str.contains(g)])
        generics_df.loc[generics_df['generic'] == g, 'n_ids_simple'] = \
            len(set(meds_small.loc[meds_small['prescription'].str.contains(g), 'id']))
assert sum(meds.loc[meds_small.index, 'prescription'] != \
           meds_small.loc[meds_small.index, 'prescription']) == 0    
meds.loc[meds_small.index, 'name_generic'] = \
    meds_small.loc[meds_small.index, 'name_generic'].copy()
print('Identified generic names in ' + \
      str(len(meds_small.loc[meds_small['name_generic'].notnull()])) + ' rows.')
meds.to_csv('meds_generics.csv', index=False, header=True, sep = '|')
generics_df.to_csv('meds_generics_stats.csv', index=False, header=True)

# C
meds['name_brand'] = None
# create a brand-generic dictionary with brands as keys and generic names as values; 
# sort alphabetically by key.
drug_dict = dict(zip(drug_names['alt'], drug_names['generic']))
drug_names_l = list(drug_dict.keys())
drug_names_l.sort(reverse = True)
drug_dict = {i: drug_dict[i] for i in drug_dict}
print('\n' + '4/5: ' + 'Replacing brand names with generic names... Brand names left: ') 
count = len(drug_dict)
print(count)
meds_small = meds.loc[meds['name_generic'].isnull()].copy()
# use regex to replace brand names in 'prescription_generic' with generic terms. 
# The regex contains two non-capture groups:
# 1st group: start of line or not character; 
# 2nd group: end of line or not character. 
# Between the groups is the brand name; this should ensure that all brand names 
# are found if they are embedded within symbols or numbers (e.g., drug concentration) 
# but not within other words, as the latter case could actually be different medicines.
for brand in drug_dict:
    meds_small.loc[(meds_small['prescription'].\
                    str.contains("(?:^|[^a-z])" + \
                                 brand + "(?:$|[^a-z])")) & \
                   (meds_small['name_generic'].isnull()), \
                       ['name_generic', 'name_brand']] = [drug_dict[brand], brand]
    count -= 1
    if (count != 0) & (count % 50 == 0):
        print(count)
    if produce_stats == True:
        drug_names.loc[drug_names['alt'] == brand, 'n_meds'] = \
            len(meds_small.loc[meds_small['name_brand'] == brand])
        drug_names.loc[drug_names['alt'] == brand, 'n_ids'] = \
            len(set(meds_small.loc[meds_small['name_brand'] == brand, 'id']))
        drug_names.loc[drug_names['alt'] == brand, 'n_meds_simple'] = \
            len(meds_small.loc[meds_small['prescription'].str.contains(brand)])
        drug_names.loc[drug_names['alt'] == brand, 'n_ids_simple'] = \
            len(set(meds_small.loc[meds_small['prescription'].str.contains(brand), 'id']))

assert sum(meds.loc[meds_small.index, 'prescription'] != \
           meds_small.loc[meds_small.index, 'prescription']) == 0
meds.loc[meds_small.index, ['name_generic', 'name_brand']] = \
    meds_small.loc[meds_small.index, ['name_generic', 'name_brand']].copy()
print('Replaced terms in ' + \
      str(len(meds_small.loc[meds_small['name_generic'].notnull()])) + ' rows.')
print('\n' + '5/5: ' + 'Finalising the data frame: categorising administration routes.' + '\n')

drug_external_brand = list(set(drug_names.loc[drug_names['external'] == 1, 'alt']))
drug_external_gen = list(set(drug_names.loc[drug_names['external'] == 2, 'generic']))
meds['application'] = 0
meds.loc[meds['name_brand'].isin(drug_external_brand), 'application'] = 1
meds.loc[meds['name_generic'].isin(drug_external_gen), 'application'] = 2
print('Administration routes based on drug names: ')
print(meds['application'].value_counts())
print('\n')


external = ['nebuliser', 'nasal', 'spray', 'gel', 'cream', 'inhalator', 'patch', 
            'plaster', 'lotion', 'shampoo', 'inhaler', 'topical', 'ophthalmic', 
            'otic', 'nose drop', 'nose drp', 'eye drp', 'eye susp', 'ear drp', 
            'oint', ' dro ', 'paste', 'drops', 'eye dro', 'ear dro', 'nose dro', 
            'swab']
for e in external:
    meds.loc[meds['prescription'].str.contains("(?:^|[^a-z])" + e + \
                                               "(?:$|[^a-z])"), 'application'] = 1    
meds.loc[meds['prescription'].str.contains("(?:^|[^a-z])" + 'tab' + \
                                           "(?:$|[^a-z])"), 'application'] = 0
print('Administration routes after additional search for relevant terms: ')
print(meds['application'].value_counts())
meds.to_csv('meds_de_branded.csv', index=False, header=True, sep = '|')
drug_names.to_csv('meds_de_branded_stats.csv', index=False, header=True)
print('\n' + 'Finished.')