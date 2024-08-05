# Pseudoscales simulation code
Code to reproduce emulation of target trial for the effect of the use of hearing aids on dementia in people with hearing loss in UKBB.

## Required files
The below table lists the UK Biobank field IDs that are required for the analyses.

Field ID | Description
----------- | -----
`31` |	Sex
`34`, `52` |	Year and month of birth
`53` |	Date of attending assessment centre
`48` |	Waist circumference
`189` |	Townsend deprivation index at recruitment
`6138` |	Qualifications
`1558` |	Alcohol use frequency
`20116` |	Smoking status
`6164` |	Physical activity
`24003`, `24004`, `24006`, `24016`, `24017`, `24018` | Air pollution
`2020` | Loneliness
`2050` | Depressive mood
`709`, `1031`, `6160` | Social isolation
`399` , `4282`, `20016`, `20018`, `20023` | Results from cognitive testing
`2247`, `2257`, `20019`, `20021`, `2257` |	Hearing impairment variables
`131258 - 131261` |	Hearing loss date
`40000` |	Date of death
`40022`, `41270`, `41271`, `41280`, `41281`, `41234` |	Hospital inpatient data
`130846` | Delirium date
`42018` |	Dementia date
`130890 - 130902` * |	Mood disorder date
`131296 - 131306` * |	Heart disease date
`131438 - 131456` * |	Influenza/pneumonia date
`131484 - 131498` * |	Lower respiratory system disorder date
`131658 - 131670` * |	Liver disease date
`131360 - 131378` * |	Cerebrovascular disease date
`130714` |	Hypertension date
`130814` |	Hyperlipidaemia date
`130874 - 130888` * |	Psychotic disorder date
`131212` | Visual impairment
`131060, 130920` | Sleep disorder
`42006`, `42028`, `130992-131032`, `131042-131050`, `131056`, `131058`, `130992-131010`, `131038-13120`, `131442` * | CNS disorder date
`130690-130748` * | Endocrine disorder date
`130750-130788` * | Nutritional deficiency date
`130798-130834` * | Metabolic disorder date
`40005`, `40006`, `40013` | Cancer
`41259` |	Hospital inpatient records ("hesin.txt")
`42038` | GP prescription records ("gp_registrations.txt")
`42039` | GP prescription records ("gp_scripts.csv")
`42040` | GP clinical event records ("gp_clinical.txt")

*_only even-numbered field IDs_

Additionally, the following files are required:
- "participant_opt_out.csv": a table with one column `id` that contains as observations the UK Biobank participant IDs for participants that have opted out of the study. This list will change over time and researchers with access to UK Biobank data will be regularly informed of additions to the list.


## The coding environment
1. Download the contents of this repository and extract them to the working directory.
2. Install R version 4.3.2 (https://cran.rstudio.com/bin/windows/base/old/4.3.2/) and Rstudio ([https://www.rstudio.com/categories/rstudio-ide/](https://posit.co/products/open-source/rstudio/)), and run the following (choose "activate the project and use the project library" when prompted):
```R
  install.packages('renv')
  renv::restore()
```
- `renv::restore()` might need to be run again to install the correct versions of the required packages.
3. Install Python via the Anaconda distribution (https://www.anaconda.com/download/success) and perform the following steps:
- (a) Launch the Anaconda prompt/terminal and navigate the the working directory.
- (b) Run `conda env create -f environment.yml` to setup the correct Python environment.
- (c) Run `conda activate pseudoscales_env` to activate the Python environment.
- The project can be discontinued and resumed later, but before resuming, always repeat steps (a) and (c).

## Running the code

1. Run the scripts with the prefixes "`0_`" to "`8_`" sequentially. Short descriptions are below (they are also available within each script):
- `0_extract_variables.R`: extract the field IDs from the UK Biobank masterfile and save them as an “.Rds” file.
- `1_` files: the final output (“data_period.csv”) indicates the periods of continuous ascertainment using the algorithm developed by Darke et al. ([![DOI](https://img.shields.io/badge/DOI-10.1093/jamia/ocab260-blue)](https://doi.org/10.1093/jamia/ocab260)).
- `2_registration_periods.py`: further processing of the output frin the previous step to enable downstream analyses. Writes `data_period_long.csv` to folder.
- `3_clean_prescriptions.py`: pre-process the UK Biobank prescriptions file. This includes - where possible - filling in missing drug names and removing invalid rows. Writes `meds_cleaned.csv` to folder.
- `4_find_generics.py`: separates combination drugs into individual compounds, replaces brand drug names with generic names, and identifies invalid routes of administration. It saves to disk the updated data frame (`meds_de_branded.csv`) frame and a file containing some statistics on the numbers of replaced drug names (`meds_de_branded_stats.csv`).
- `5_prepare_covariates.R`: extracts and cleans the covariates and other variables used in downstream analyses and saves them as `covariates.Rds`. Also updates the ABS masterfile (`aas_combined.csv`) to exclude drugs not prescribed in the period of interest and saves the new file (`aas_in_sample.csv`) to disk.
- `6_generate_abs.py`: calculates for each participant their annual cumulative anticholinergic burden according to each ABS; saves the results to disk as `achb_scales.csv`. Also saves the non-anticholinergic burden according to each scale (to be used as covariate) as `achb_scales_poly_csv`.
- `7_simulate_pseudoscales.py`: defines the function to perform the pseudoscales simulation. The details of the arguments are provided in the script. 
- `8_1_run_simulation.py`: executes the function to reproduce the pseudoscales in the manuscript. The code runs `pseudo_scales()` (defined in the step above) 4x: (1) across-sampling for general polypharmacy, (2) across-sampling for anticholinergic polypharmacy, (3) within-sampling for general polypharmacy, and (4) within-sampling for anticholinergic polypharmacy.
- `8_2_modelling.R`: executes the functions `combine_scales()`, `prepare_scales()`, and `outcome_effect_parallel()` for the four sets of simulated pseudoscales from the previous step. The functions have high RAM requirements due to the size of the data. The "core_number" argument for `outcome_effect_parallel` might need to be altered, depending on the machine that the code is being run on.
