# brain_states-analyses
R code to generate results for Gui et al., Diminished engagement of attentive brain states predicts later Autism Spectrum Disorder


The authors do not have permission to share the real data. 

pheno_dataset.csv is a .csv file with one row per participant and one column per phenotypic variable. The selected variables for these analyses were: Outcome group, age in days, sex , Mullen Early Learning Composite at 8 months, VABS Socialization Score at 3 years, VABS Motor Skills Score at 3 years)

NrTrials.csv is a pre-processing record .csv file with one row per participant and one column per condition indicating the number of valid trials per participant

The ~/ERPs_folder/ contains one .txt file per participant per condition (labelled as ID_Condition.txt, i.e., 001_FaceDirect.txt) with the entire time series of pre-processed multichannel ERP data (498 time-stamps as rows, 124 channels as columns).
