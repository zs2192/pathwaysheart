# Pathways Heart Study analytical file repository

Project setup:


Current analytical datasets:

Cases data received on 2019-06-27: 
  cases_final_27mar19.sas7bdat; 
  cases_tumor_char_26jun19.sas7bdat

Controls data received on 2019-07-03: 
  controls_group1.sas7bdat; 
  controls_group2.sas7bdat

Risk factor data received on 4/17/2019: 
  bmi.sas7bdat; 
  bp.sas7bdat; 
  dyslipidemia.sas7bdat; 
  diabetes.sas7bdat; 
  smoking.sas7bdat; 
  smoking_6months.sas7bdat

Menopause, parity, census data received on 4/17/2019: 
  menopause.sas7bdat; 
  parity.sas7bdat; 
  census.sas7bdat

CVD outcome data recevied on 6/27/2019: 
  cvd_events.sas7bdat

censoring events received on 6/27/2019: 
  censoring_27jun19.sas7bdat

lab data received on 6/27/2019: 
  labs.sas7bdat



There are 2 R codes:

  0_data.R is for data importing, cleaning, and variable definition. Please replace the data file path with your own local path.
  1_analysis.R is for data analysis of Aim 1. This code will generate the tables for manuscript.
  






