# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 12:02:08 2019

@author: jurem

The code below:
    1. Supplements the rows that are missing drug codes with codes from the Read_v2-reader.
    2. For drugs that were classified as anticholinergic by at least one scale 
    (and for other common drugs) substitutes the brand names by generic drug names.
    3. Removes from the sample (a) participants that have opted out of the study, 
    (b) rows without prescriptions, (c) rows without dates.

The columns retained in the exported data frame are:
    - 'id': participant id
    - 'data_provider': 1 = England (Vision), 2 = Scotland, 3 = England (TPP), 4 = Wales
    - 'date': when the prescription was issued
    - 'prescription': full name of the prescription
    - 'quantity': quantity/dose of the prescribed drug, usually in mg
"""

## helper function to remove the additional 0's in some read-codes
def remove_00(code):
    # do not change the 'unknown'-strings, change only those with two 0s at the end
    if (code != 'unknown') & (len(code)==7) & (code[-2:len(code)]=='00'):
        # retain everything but the last two characters of the string
        new_code = code[0:-2] 
        return new_code
    else:
        return code


import pandas as pd
import numpy as np


## Read in the data and prepare it.
meds = pd.read_csv('gp_scripts.csv', sep='\t', encoding='cp1252')
# re-name the columns
meds.columns = ['id', 'data_provider', 'date', 'read_code', 'bnf', 'dmd', 
                'prescription', 'quantity'] 
# create a backup for the original prescription column
meds.loc[:,'prescription_old'] = meds.loc[:,'prescription'].copy() 
# to lower case
meds.loc[:,'prescription'] = meds.loc[:,'prescription'].str.lower()
# read in the codes
codes = pd.read_csv('read_codes.csv', sep=",", dtype = str, encoding = "cp1252")
# re-name the columns
codes.columns = ['code', 'drug', 'status_flag'] 
# change codes to strings
codes['code'] = codes.code.astype(str) 
# remove leading and trailing white spaces from read-codes in the read-code data frame
codes['code'] = codes.loc[:,'code'].apply(str.strip)
# change NA values in read_code column into 'unknown'
meds.loc[meds['read_code'].isna(), 'read_code'] = 'unknown' 
# remove white spaces from read-codes in the prescriptions data frame
meds['read_code'] = meds.loc[:,'read_code'].apply(str.strip) 




# Some read-codes contain two 0s at the end; remove those.
meds['read_code'] = meds['read_code'].apply(remove_00)




## Use the read-code list to supplement the data frame.
# create a dictionary with read-code/drug-name pairs
read_code_dict = (codes.groupby('code')['drug'].apply(lambda x: x.tolist())).to_dict()
# the prescriptions are individual lists; transform them to strings
for key in read_code_dict.keys():
    read_code_dict[key] = ''.join(read_code_dict[key])

# helper code to fill read-codes into the 'meds' dataset
def find_read_code(code):
    try:
        read_code = read_code_dict[code]
    except: # if the code doesn't exist, flag as "unknown"
        read_code = 'unknown'
    return read_code

# create a column for read-code-supplemented information
meds['prescription_read'] = np.nan
# plug each read-code in our sample into the dictionary as a key and create an 
# additional column from the values
meds['prescription_read'] = meds.loc[:,'read_code'].apply(lambda x: find_read_code(x))
#convert to lowercase
meds.loc[:,'prescription_read'] = meds.loc[:,'prescription_read'].str.lower()
# put read-code-supplied drugs into the drug column
meds.loc[meds['prescription'].isna(), 'prescription'] = \
    (meds.loc[meds['prescription'].isna(), 'prescription_read']).copy()
# change 'unknown' in prescription column back to NaN
meds.loc[meds['prescription']=='unknown', 'prescription'] = np.nan



# remove participants that have opted outs
opt_out = pd.read_csv('participant_opt_out.csv')
opt_out.columns = ['id']
meds['id'] = meds['id'].astype(str)
opt_out['id'] = opt_out['id'].astype(str)
meds = meds.loc[~meds['id'].isin(opt_out.id)] 
meds = meds.reset_index(drop=True)
# remove rows with blank prescription column ()
meds = meds.loc[meds['prescription'].notnull()]
# remove rows with blank/invalid date (3)
meds = meds.loc[meds['date'].notnull(), :]
invalid_dates = ["01/01/1901", "02/02/1902", "03/03/1903", "07/07/2037"]
meds = meds.loc[~meds['date'].isin(invalid_dates), :]
#remove unnecessary columns
meds.drop(['read_code','bnf','dmd','prescription_read'], axis=1, inplace=True)
# change all prescriptions to strings
meds['prescription'] = meds.prescription.astype(str)
#remove potential white-space from the front of prescription names
meds['prescription'] = meds.loc[:,'prescription'].apply(str.strip)
#remove all '|' characters, as it will be used as a column separator
meds['prescription'] = meds['prescription'].str.replace('|', ' ', regex = False)

#export to .csv
meds.to_csv('output_files/meds_cleaned.csv', index=False, header=True, sep='|')