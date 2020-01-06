# Pathways Heart Study analytical file repository

## Current anlayses are based on data received on 11/12/2019:

### Cases with tumor characteristics: 
    cases.sas7bdat; 
### Controls: 
    controls.sas7bdat; 

### Baseline risk factors: 
    baseline_bmi.sas7bdat; 
    baseline_bp.sas7bdat; 
    dyslipidemia.sas7bdat; 
    diabetes.sas7bdat; 
    baseline_labs.sas7bdat
    smoking.sas7bdat; 
    menopause.sas7bdat; 
    census.sas7bdat

### CVD outcome: 
    cvd_events.sas7bdat

### Censoring events: 
    censoring.sas7bdat

## There are 2 R codes:
    0_data.R is for data importing, cleaning, and variable definition. Please replace the data file path with your own local path.
    
    1_analysis.R is for data analysis of Aim 1. This code will generate the tables for manuscript.
  






