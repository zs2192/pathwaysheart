# Pathways Heart Study Aim 1 data cleaning
# Zaixing Shi, 1/7/2020


library(tidyverse)
library(haven)
library(ggplot2)
library(dplyr)
library(plyr)
library(naniar)
library(reshape2)
library(doParallel)
library(lsr)

# setup parallel computing
nodes <- detectCores()
cl <- makeCluster(nodes)
registerDoParallel(cl)

################################################################################
# Load data

# import cases data received on 2019-11-12
pathpref = 'data/20191112105146/'
cases <- as.data.frame(read_sas(paste0(pathpref,'cases.sas7bdat')))

# import controls data received on 2019-11-12
controls <- as.data.frame(read_sas(paste0(pathpref,'controls.sas7bdat')))

## risk factor data received on 2019-11-12
bmi <- as.data.frame(read_sas(paste0(pathpref,'baseline_bmi.sas7bdat')))
bp <- as.data.frame(read_sas(paste0(pathpref,'baseline_bp.sas7bdat')))
lipid <- as.data.frame(read_sas(paste0(pathpref,'dyslipidemia.sas7bdat')))
diab <- as.data.frame(read_sas(paste0(pathpref,'diabetes.sas7bdat')))
smok <- as.data.frame(read_sas(paste0(pathpref,'smoking.sas7bdat')))
hyper <- as.data.frame(read_sas(paste0(pathpref,'phase_hypertension.sas7bdat')))

## menopause, parity, census data received on 2019-11-12
menop <- as.data.frame(read_sas(paste0(pathpref,'menopause.sas7bdat')))
census <- as.data.frame(read_sas(paste0(pathpref,'census.sas7bdat')))

# import CVD outcome data recevied on 2019-11-12
cvd <- as.data.frame(read_sas(paste0(pathpref,'cvd_events.sas7bdat')))

# import censoring events received on 2019-11-12
censor <- as.data.frame(read_sas(paste0(pathpref,'censoring.sas7bdat')))

# import updated lab data received on 2019-11-12 
labs <- as.data.frame(read_sas(paste0(pathpref,'baseline_labs.sas7bdat')))






################################################################################
# Combine cases and controls data

# check missingness
png('case_missing.png', height=3.5, width = 2.5, unit='in', res=200)
gg_miss_var(cases, show_pct =TRUE) + labs(y = "% missing, Cases")
dev.off()

png('control_missing.png', height=3.5, width = 2.5, unit='in', res=200)
gg_miss_var(controls, show_pct =TRUE) + labs(y = "% missing, Controls")
dev.off()

# change all variables to lower case 
names(cases) <- tolower(names(cases))
names(controls) <- tolower(names(controls))


# create another surgery date var for cases and controls for risk factor data selection
cases$rf_date <- cases$daysto_surg
## fill patients without surg_date with dxdate plus median days between dx to surg
cases$rf_date[which(is.na(cases$rf_date))] <- as.numeric(median(cases$daysto_surg, na.rm=T))

controls$rf_date <- cases$rf_date[match(controls$num7_caseid,cases$num7_studyid)]

# stack cases and controls, matched on age and race/ethnicity
# reorder controls1 variable to be the same as cases
controls$group <- 'Control'
controls$dxage <- controls$index_age
controls[,names(cases)[!(names(cases) %in% names(controls))]] <- NA

cases[,names(controls)[!(names(controls) %in% names(cases))]] <- NA
cases$group <- 'Case'
controls <- controls[,names(cases)]

# combine cases and control data, n=89644
all <- rbind(cases, controls)






################################################################################
# Select labs data
# The following steps are recommended because the original data may include data 
# that fall outside the date range - within 36 mos prior to index date and before 
# surgery


# BMI
names(bmi) <- tolower(names(bmi))
## merge the ref_date to the BMI dataset
bmi <- merge(bmi, all[,c("num7_studyid",'rf_date')], by='num7_studyid')
## select the bmi measures within 3 years prior to index date and before the ref_date
bmi1 <- bmi[which(bmi$daysto_bmi >= -365*3 & bmi$daysto_bmi <= bmi$rf_date),]
## calculate the date absolute difference between index date and each measurement date
#bmi1$datediff <- abs(as.numeric(bmi1$measure_date - bmi1$index_date))
## for each patient, select the measurement with the minimum absolute date difference
#bmi1 <- ddply(bmi1, .(num7_studyid), function(d){
#  d <- d[which(d$datediff==min(d$datediff)),]
#  d
#},.parallel=TRUE)
## clean up BMI to remove non-numeric characters from data
#bmi1$bmi_num <- gsub('<|>|>=','',bmi1$bmi)
#bmi1$bmi_num <- as.numeric(gsub('-.*','',bmi1$bmi_num))
## aggregate multiple values within same day
#bmi2 <- aggregate(x = bmi1$bmi_num, 
#                  by = list(bmi1$num7_studyid,bmi1$datediff), 
#                  mean, na.rm=T)
## rename the final dataset
#names(bmi2) <- c('num7_studyid','bmi_datediff','bmi')




# BP
names(bp) <- tolower(names(bp))
bp <- merge(bp, all[,c("num7_studyid",'rf_date')], by='num7_studyid')
bp1 <- bp[which(bp$daysto_bp >= -365*2 & bp$daysto_bp <= bp$rf_date),]
#bp1$datediff <- abs(as.numeric(bp1$measure_date - bp1$index_date))
#bp1 <- ddply(bp1, .(num7_studyid), function(d){
#  d <- d[which(d$datediff==min(d$datediff)),]
#  d
#},.parallel=TRUE)

## clean up bp
#bp1$hypertension <- factor(bp1$hypertension,
#                           labels = c('Normal','Elevated','Stage 1','Stage 2','Undefined'))

## aggregate multiple values within same day
#systol <- aggregate(x = bp1$systolic, 
#                  by = list(bp1$num7_studyid,bp1$datediff), 
#                  mean, na.rm=T)
#diastol <- aggregate(x = bp1$diastolic, 
#                   by = list(bp1$num7_studyid,bp1$datediff), 
#                   mean, na.rm=T)

#bp2 <- merge(systol, diastol, by=c("Group.1","Group.2"))
#names(bp2) <- c('num7_studyid','bp_datediff','systolic','diastolic')



# Labs
labs$result <- as.numeric(labs$RESULT_C)

# merge surgery date
names(labs) <- tolower(names(labs))
labs <- merge(labs, all[,c("num7_studyid",'rf_date')], by='num7_studyid')
labs1 <- labs[which(labs$daysto_labs >= -365*3 & labs$daysto_labs <= labs$rf_date),]
#labs1$datediff <- abs(as.numeric(labs1$test_dt - labs1$index_date))
#system.time(labs1 <- ddply(labs1, .(num7_studyid,test_type), function(d){
#  d <- d[which(d$datediff==min(d$datediff)),]
#  d
#},.parallel=TRUE))

## Turn lab data from long to wide format, so that each patient has 1 row, and 
## each type of test is in a separate column

labs_w <- dcast(data=labs,num7_studyid~test_type,
                value.var = 'result', mean, na.rm=T)

# save a copy of labs2 data
write_csv(labs_w,'labs wide.csv')





################################################################################
# Clean CVD outcome data

## make long to wide format, keep 1st occurence of each condition only
### condition groups

names(cvd) <- tolower(names(cvd))
# recode cvd group
cvd$cvd_group[cvd$cvd_condition=='Heart Failure'] <- 'Heart Failure'
cvd$cvd_group[cvd$cvd_condition=='Cardiomyopathy'] <- 'Cardiomyopathy'
cvd$cvd_group[cvd$cvd_condition=='Transient Ischemic Attack (TIA)'] <- 'TIA'
cvd$cvd_group[cvd$cvd_condition %in% c('Acute Ischemic stroke',
                                       'Acute myocardial infarction',
                                       'Intracerebral hemorrhage',
                                       'Other cerebrovascular disease',
                                       'Retinal vascular occlusion',
                                       'Subarachnoid hemorrhage')] <- 'Stroke'

#### define true incidence and prevalence
cvd2 <- dcast(data=cvd, num7_studyid ~ cvd_group,
              value.var = 'daysto_cvd_event', min,na.rm=T)

cvd2[,-1] <- lapply(cvd2[,-1], function(x) {
  x[x==Inf] <- NA
  #x <- as.Date(x, origin='1970-01-01')
  x
})

names(cvd2) <- gsub('/| ','_',names(cvd2))
cvd_grp <- paste0(names(cvd2)[-1],'_grp')
names(cvd2)[-1] <- paste0(cvd_grp,'_dt')

cvd2[,cvd_grp] <- lapply(cvd2[,-1], function(x) {
  y <- ifelse(!is.na(x),1,0)
  y[y==1 & x>0] <- 2
  y <- factor(y, levels=0:2, labels = c('No','Prevalent','Incident'))
  y
  })

cvd2list <- lapply(cvd_grp, function(x){
  dates <- eval(parse(text=paste0("dcast(data=cvd2, num7_studyid~",x,", value.var = '",x,"_dt')")))
  prev <- ifelse(cvd2[,x]=='Prevalent',1,0)
  inc <- ifelse(cvd2[,x]=='Incident',1,0)
  sum <- data.frame(cbind(prev,inc,dates))
  sum <- sum[,c('num7_studyid','prev','inc','Incident')]
  names(sum) <- c('num7_studyid',paste0(x,c('_prev','_inc','_incdt')))
  sum
})

#### define recurrence
cvd2r <- dcast(data=cvd[which(cvd$daysto_cvd_event>=0),], 
               num7_studyid~cvd_group, value.var = 'daysto_cvd_event', min,na.rm=T)

cvd2r[,-1] <- lapply(cvd2r[,-1], function(x) {
  x[x==Inf] <- NA
  x
  })

names(cvd2r) <- gsub('/| ','_',names(cvd2r))
cvd_grpr <- paste0(names(cvd2r)[-1],'_grp_rec')
names(cvd2r)[-1] <- paste0(cvd_grpr,'dt')

cvd2r[,cvd_grpr] <- lapply(cvd2r[,-1], function(x) {
  y <- ifelse(!is.na(x),1,0)
  y
})

#### merge all together
cvd2list[[12]] <- cvd2
cvd2list[[13]] <- cvd2r
cvd2_final <- Reduce(function(...) merge(...,by='num7_studyid', all.x=T), cvd2list)




### condition

#### define true incidence and prevalence
cvd3 <- dcast(data=cvd, num7_studyid~cvd_condition, 
              value.var = 'daysto_cvd_event', min,na.rm=T)


names(cvd3) <- gsub('/| ','_',names(cvd3))
names(cvd3) <- gsub('&_|\\(|\\)','',names(cvd3))
names(cvd3)[17] <- "Percutaneous_coronary_intervention_PCI"

cvd3[,-1] <- lapply(cvd3[,-1], function(x) {
  x[x==Inf] <- NA
  x
})

cvd_cond <- names(cvd3)[-1]
names(cvd3)[-1] <- paste0(names(cvd3)[-1],'_dt')
cvd3[,cvd_cond] <- lapply(cvd3[,-1],  function(x) {
  y <- ifelse(!is.na(x),1,0)
  y[y==1 & x>0] <- 2
  y <- factor(y, levels=0:2, labels = c('No','Prevalent','Incident'))
  y
})

cvd3list <- lapply(cvd_cond, function(x){
  dates <- eval(parse(text=paste0("dcast(data=cvd3, num7_studyid~",x,", value.var = '",x,"_dt')")))
  prev <- ifelse(cvd3[,x]=='Prevalent',1,0)
  inc <- ifelse(cvd3[,x]=='Incident',1,0)
  sum <- data.frame(cbind(prev,inc,dates))
  sum <- sum[,c('num7_studyid','prev','inc','Incident')]
  names(sum) <- c('num7_studyid',paste0(x,c('_prev','_inc','_incdt')))
  sum
})

#### define recurrence
cvd3r <- dcast(data=cvd[which(cvd$daysto_cvd_event>=0),], 
               num7_studyid~cvd_condition, value.var = 'daysto_cvd_event', min,na.rm=T)

cvd3r[,-1] <- lapply(cvd3r[,-1], function(x) {
  x[x==Inf] <- NA
  x
})

names(cvd3r) <- gsub('/| ','_',names(cvd3r))
names(cvd3r) <- gsub('&_|\\(|\\)','',names(cvd3r))
names(cvd3r)[17] <- "Percutaneous_coronary_intervention_PCI"
cvd_condr <- paste0(names(cvd3r)[-1],'_rec')
names(cvd3r)[-1] <- paste0(cvd_condr,'dt')

cvd3r[,cvd_condr] <- lapply(cvd3r[,-1], function(x) {
  y <- ifelse(!is.na(x),1,0)
  y
})


cvd3list[[23]] <- cvd3
cvd3list[[24]] <- cvd3r
cvd3_final <- Reduce(function(...) merge(...,by='num7_studyid', all.x=T), cvd3list)

# recode recurrence patients with prevalent 


# export condition list
write_csv(cvd_grp,'cvd_grp.csv')
write_csv(cvd_cond,'cvd_cond.csv')

# save a copy of the cvd data
write_csv(cvd2_final,'cvd2_final.csv')
write_csv(cvd3_final,'cvd3_final.csv')







################################################################################
# Merge all data

## combine all data into a list
datalist <- list(all, 
                 bmi[,c("num7_studyid","bmi","bmi_category_ii","daysto_bmi")],
                 bp[,c("num7_studyid","systolic","diastolic","daysto_bp","hypertension")],
                 labs_w[,c("num7_studyid","GLU_F","GTT75_PRE","HDL","HGBA1C","LDL_CLC_NS","TOT_CHOLES","TRIGL_NS")],
                 diab[,c(2:4)],lipid[,1:3], hyper[,c(2,3,5)],
                 smok[,c(2:3)],menop,
                 census[,c(1:14)],
                 cvd2_final,
                 cvd3_final,
                 censor[,-2])

## make all var names lower case
datalist <- lapply(datalist,function(x){
  names(x) <- tolower(names(x))
  x
})

## merge all data by num7_studyid
a1 <- Reduce(function(...) merge(...,by='num7_studyid', all.x=T), datalist)

## set all NA for cvd outcomes as 0
a1[,tolower(c(cvd_grp,cvd_cond))] <- lapply(a1[,tolower(c(cvd_grp,cvd_cond))], function(x) {
  x[is.na(x)] <- 'No'
  x
})

a1[,grep('inc$|prev$|rec$',names(a1))] <- 
  lapply(a1[,grep('inc$|prev$|rec$',names(a1))], function(x) {
    x[is.na(x)] <- 0
    x
  })

# save a copy of the merged data
write_csv(a1,'a1.csv')







################################################################################
# Create some more variables...

# recode categorical variables
# age groups
a1$agegrp <- cut(a1$dxage, c(0, 40, 50, 60, 70, 101), right=F,
                  labels = c('<40 yo','40-49','50-59','60-69','70+ yo'))

# race
a1$raceethn1 <- factor(a1$raceethn1, 
                        labels=c("WHITE","BLACK","ASIAN","HISPANIC","PI","AI-AN"))

# postive nodes
a1$nodal <- factor(a1$nodal, labels = c('Positive','Negative'))

# tumor markers
a1[,c('er','pr','her2')] <- lapply(a1[,c('er','pr','her2')], function(x)
  factor(x, levels = c(0,1,2,3,8,9),
         labels = c('Not done','Positive','Negative',
                    'Borderline','Ordered, N/A',
                    'Unknown')))

# erpr status
a1$erpr <- factor(a1$erpr, labels = c('ER+/PR+','ER+/PR-','ER-/PR+','ER-/PR-',
                                        'UNKNOWN'))

# more tumor markers and treatment received
a1[,c('tri_neg','chemo_yn','rad_yn','horm_yn')] <- 
  lapply(a1[,c('tri_neg','chemo_yn','rad_yn','horm_yn')], function(x)
    factor(x, levels = c(0,1,9),
           labels = c('No','Yes','Other')))

# AJCC stage
a1$ajcc_stage <- factor(a1$ajcc_stage,
                         labels = c('Stage I','Stage II', 'Stage III','Stage IV'))

# create enrollment length
a1$enr_len <- as.numeric((a1$daysto_enr_end - a1$daysto_enr_start)/30)


# BMI category
a1$bmi1 <- as.numeric(gsub('[<|>|=]','',a1$bmi))
a1$bmicat <- cut(a1$bmi1, c(0,18.5,25,30,35,Inf),right=F,
                 labels = c('Underweight','Normal','Overweight','Obese I','Obese II+'))
a1$bmicat <- as.character(a1$bmicat)
a1$bmicat[a1$raceethn1=='ASIAN' & (a1$bmi1>=23 & a1$bmi1<27.5)] <- 'Overweight'
a1$bmicat[a1$raceethn1=='ASIAN' & (a1$bmi1>=27.5 & a1$bmi1<35)] <- 'Obese I'
a1$bmicat[is.na(a1$bmicat)] <- 'Unknown'
a1$bmicat <- factor(a1$bmicat,
                    levels = c('Underweight','Normal','Overweight','Obese I','Obese II+','Unknown'))

# smoking status
a1$smok <- a1$smoke_status_6m
a1$smok[a1$smok %in% c('Q')] <- 'Quited'
a1$smok[a1$smok %in% c('Y')] <- 'Current smoker'
a1$smok[a1$smok %in% c('N')] <- 'Never smoker'
a1$smok[(a1$smok %in% c(''))|is.na(a1$smok)] <- 'Unknown'

a1$smok <- factor(a1$smok, levels=c('Never smoker','Current smoker','Quited','Unknown'))

# menopausal status
a1$menop <- factor(a1$bl_meno_status)
a1$menop[is.na(a1$menop) & a1$dxage>51] <- '1'
a1$menop[is.na(a1$menop) & a1$dxage<=51] <- '0'

# baseline diabetes                                                 })
a1$diab_bl <- ifelse(a1$daysto_diabetes<=0, 1, 0)
a1$diab_bl[is.na(a1$diab_bl)] <- 0
a1$diab_bl <- factor(a1$diab_bl,levels=c(1,0),labels = c('Yes','No'))

# baseline dyslipidemia
a1$dyslipid_bl <- ifelse(a1$daysto_dyslipidemia<=0, 1, 0)
a1$dyslipid_bl[is.na(a1$dyslipid_bl)] <- 0
a1$dyslipid_bl <- factor(a1$dyslipid_bl,levels=c(1,0),labels = c('Yes','No'))

# baseline hypertension
a1$htn_bl <- ifelse(a1$daysto_phase_htn<=0, 1, 0)
a1$htn_bl[is.na(a1$htn_bl)] <- 0
a1$htn_bl <- factor(a1$htn_bl,levels=c(1,0),labels = c('Yes','No'))

# household income
a1$hhincome1 <- rowSums(a1[,grep('famincome[1-9]$',names(a1), value=T)], na.rm = T)
a1$hhincome2 <- rowSums(a1[,grep('famincome1[0-2]',names(a1), value=T)], na.rm = T)
a1$hhincome3 <- rowSums(a1[,grep('famincome1[3-6]',names(a1), value=T)], na.rm = T)

# education
a1$edu1 <- rowSums(a1[,grep('education[1-3]$',names(a1), value=T)], na.rm = T)
a1$edu2 <- rowSums(a1[,grep('education[4-5]$',names(a1), value=T)], na.rm = T)
a1$edu3 <- a1$education6
a1$edu4 <- rowSums(a1[,grep('education[7-8]$',names(a1), value=T)], na.rm = T)








################################################################################
# Define CVD outcomes


# define true incindence
# combined event of ischemic heart disease, stroke/TIA, cardiomyopathy/heart failure
a1$cvdcombo_grp_prev <- ifelse(a1$ischemic_heart_disease_grp_prev==1 | 
                                a1$stroke_grp_prev==1 | 
                                a1$cardiomyopathy_grp_prev==1|
                                 a1$heart_failure_grp_prev==1, 1, 0)


a1$cvdcombo_grp_inc <- ifelse(a1$ischemic_heart_disease_grp_inc==1 | 
                                a1$stroke_grp_inc==1 |  
                                a1$cardiomyopathy_grp_inc==1|
                                a1$heart_failure_grp_inc==1, 1, 0)


a1$cvdcombo_grp_incdt <- apply(a1[,c("ischemic_heart_disease_grp_incdt",
                                     "stroke_grp_incdt",
                                     "cardiomyopathy_grp_incdt",
                                     "heart_failure_grp_incdt"
                                     )],1,min,na.rm=T)


# define censoring time
a1$censor_dt <- apply(a1[,c("daysto_death","daysto_enr_end","daysto_end_study")],1,min,na.rm=T)


# calculate time of follow up
a1$ischemic_heart_disease_grp_inc_fu <- a1$ischemic_heart_disease_grp_incdt
a1$ischemic_heart_disease_grp_inc_fu[which(a1$ischemic_heart_disease_grp_inc==0)] <- 
  a1$censor_dt[which(a1$ischemic_heart_disease_grp_inc==0)] 

a1$stroke_grp_inc_fu <- a1$stroke_grp_incdt
a1$stroke_grp_inc_fu[which(a1$stroke_grp_inc==0)] <- 
  a1$censor_dt[which(a1$stroke_grp_inc==0)]

a1$cardiomyopathy_grp_inc_fu <- a1$cardiomyopathy_grp_incdt
a1$cardiomyopathy_grp_inc_fu[which(a1$cardiomyopathy_grp_inc==0)] <- 
  a1$censor_dt[which(a1$cardiomyopathy_grp_inc==0)] 

a1$heart_failure_grp_inc_fu <- a1$heart_failure_grp_incdt
a1$heart_failure_grp_inc_fu[which(a1$heart_failure_grp_inc==0)] <- 
  a1$censor_dt[which(a1$heart_failure_grp_inc==0)] 

a1$cvdcombo_grp_inc_fu <- a1$cvdcombo_grp_incdt
a1$cvdcombo_grp_inc_fu[which(a1$cvdcombo_grp_inc==0)] <- 
  a1$censor_dt[which(a1$cvdcombo_grp_inc==0)] 

a1$arrhythmia_grp_inc_fu <- a1$arrhythmia_grp_incdt
a1$arrhythmia_grp_inc_fu[which(a1$arrhythmia_grp_inc==0)] <- 
  a1$censor_dt[which(a1$arrhythmia_grp_inc==0)] 

a1$cardiac_arrest_grp_inc_fu <- a1$cardiac_arrest_grp_incdt
a1$cardiac_arrest_grp_inc_fu[which(a1$cardiac_arrest_grp_inc==0)] <- 
  a1$censor_dt[which(a1$cardiac_arrest_grp_inc==0)] 

a1$carotid_disease_grp_inc_fu <- a1$carotid_disease_grp_incdt
a1$carotid_disease_grp_inc_fu[which(a1$carotid_disease_grp_inc==0)] <- 
  a1$censor_dt[which(a1$carotid_disease_grp_inc==0)] 

a1$myocarditis_pericarditis_grp_inc_fu <- a1$myocarditis_pericarditis_grp_incdt
a1$myocarditis_pericarditis_grp_inc_fu[which(a1$myocarditis_pericarditis_grp_inc==0)] <- 
  a1$censor_dt[which(a1$myocarditis_pericarditis_grp_inc==0)] 

a1$tia_grp_inc_fu <- a1$tia_grp_incdt
a1$tia_grp_inc_fu[which(a1$tia_grp_inc==0)] <- 
  a1$censor_dt[which(a1$tia_grp_inc==0)] 

a1$valvular_disease_grp_inc_fu <- a1$valvular_disease_grp_incdt
a1$valvular_disease_grp_inc_fu[which(a1$valvular_disease_grp_inc==0)] <- 
  a1$censor_dt[which(a1$valvular_disease_grp_inc==0)] 

a1$venous_thromboembolic_disease_grp_inc_fu <- a1$venous_thromboembolic_disease_grp_incdt
a1$venous_thromboembolic_disease_grp_inc_fu[which(a1$venous_thromboembolic_disease_grp_inc==0)] <- 
  a1$censor_dt[which(a1$venous_thromboembolic_disease_grp_inc==0)] 


# define any new onset = true incindence + recurrence
# combined event of ischemic heart disease, stroke/TIA, cardiomyopathy/heart failure
a1$cvdcombo_grp_rec <- ifelse(a1$ischemic_heart_disease_grp_rec==1 | 
                                a1$stroke_grp_rec==1 |  
                                a1$cardiomyopathy_grp_rec==1|
                                a1$heart_failure_grp_rec==1, 1, 0)


a1$cvdcombo_grp_recdt <- apply(a1[,c("ischemic_heart_disease_grp_recdt",
                                     "stroke_grp_recdt",
                                     "cardiomyopathy_grp_recdt",
                                     "heart_failure_grp_recdt"
)],1,min,na.rm=T)



# calculate time of follow up
a1$ischemic_heart_disease_grp_rec_fu <- a1$ischemic_heart_disease_grp_recdt
a1$ischemic_heart_disease_grp_rec_fu[which(a1$ischemic_heart_disease_grp_rec==0)] <- 
  a1$censor_dt[which(a1$ischemic_heart_disease_grp_rec==0)] 

a1$stroke_grp_rec_fu <- a1$stroke_grp_recdt
a1$stroke_grp_rec_fu[which(a1$stroke_grp_rec==0)] <- 
  a1$censor_dt[which(a1$stroke_grp_rec==0)]

a1$cardiomyopathy_grp_rec_fu <- a1$cardiomyopathy_grp_recdt
a1$cardiomyopathy_grp_rec_fu[which(a1$cardiomyopathy_grp_rec==0)] <- 
  a1$censor_dt[which(a1$cardiomyopathy_grp_rec==0)] 

a1$heart_failure_grp_rec_fu <- a1$heart_failure_grp_recdt
a1$heart_failure_grp_rec_fu[which(a1$heart_failure_grp_rec==0)] <- 
  a1$censor_dt[which(a1$heart_failure_grp_rec==0)] 

a1$cvdcombo_grp_rec_fu <- a1$cvdcombo_grp_recdt
a1$cvdcombo_grp_rec_fu[which(a1$cvdcombo_grp_rec==0)] <- 
  a1$censor_dt[which(a1$cvdcombo_grp_rec==0)] 








################################################################################
# Define CVD risk factor outcomes

# define incident diabetes, hypertension and dyslipidemia
a1$cvdrf_diab <- ifelse(a1$daysto_diabetes>0,1,0)
a1$cvdrf_diab[is.na(a1$cvdrf_diab)] <- 0
a1$cvdrf_htn <- ifelse(a1$daysto_phase_htn>0,1,0)
a1$cvdrf_htn[is.na(a1$cvdrf_htn)] <- 0
a1$cvdrf_dyslipid <- ifelse(a1$daysto_dyslipidemia>0,1,0)
a1$cvdrf_dyslipid[is.na(a1$cvdrf_dyslipid)] <- 0


# combined event of diabetes, hypertension or dyslipidemia
a1$cvdrfcombo <- ifelse(a1$cvdrf_diab==1 | 
                                a1$cvdrf_htn==1 | 
                                a1$cvdrf_dyslipid==1, 1, 0)


a1$cvdrfcombo_dt <- apply(a1[,c("daysto_diabetes","daysto_phase_htn",
                                     "daysto_dyslipidemia")],1,function(x){
                                       min(x[x>0], na.rm=T)
                                     })
a1$cvdrfcombo_dt[a1$cvdrfcombo_dt==Inf] <- NA

# define censoring time
a1$censor_dt <- apply(a1[,c("daysto_death","daysto_enr_end","daysto_end_study")],1,min,na.rm=T)


# calculate time of follow up
a1$cvdrf_diab_fu <- a1$daysto_diabetes
a1$cvdrf_diab_fu[which(a1$cvdrf_diab==0)] <- a1$censor_dt[which(a1$cvdrf_diab==0)] 

a1$cvdrf_htn_fu <- a1$daysto_phase_htn
a1$cvdrf_htn_fu[which(a1$cvdrf_htn==0)] <- a1$censor_dt[which(a1$cvdrf_htn==0)] 

a1$cvdrf_dyslipid_fu <- a1$daysto_dyslipidemia
a1$cvdrf_dyslipid_fu[which(a1$cvdrf_dyslipid==0)] <- a1$censor_dt[which(a1$cvdrf_dyslipid==0)] 

a1$cvdrfcombo_fu <- a1$cvdrfcombo_dt
a1$cvdrfcombo_fu[which(a1$cvdrfcombo==0)] <- a1$censor_dt[which(a1$cvdrfcombo==0)] 









################################################################################
# Data process for modeling

# relevel the group var, with control as referent
a1$group1 <- factor(a1$group, levels = c('Control','Case'))

# relevel bmi var
a1$bmicat1 <- factor(a1$bmicat, levels = c('Normal','Underweight','Overweight','Obese I','Obese II+','Unknown'))

# relevel baseline diabetes 
a1$diab_bl1 <- factor(a1$diab_bl, levels = c('No','Yes'))

# relevel dyslipidemia
a1$dyslipid_bl1 <- factor(a1$dyslipid_bl, levels=c('No/Unknown','Yes'))

# cut continuous var into quantiles
# variables to analyze
varlist <- c('enr_len','cops2','bmi1','agegrp','raceethn1',
             'bmicat','smok','menop','diab_bl','htn_bl','dyslipid_bl',
             "systolic" ,"diastolic","glu_f","hdl","hgba1c",
             "ldl_clc_ns","tot_choles","trigl_ns",
             'medhousincome','houspoverty', "edu1","edu2","edu3","edu4")

a1$glu_q <- quantileCut(a1$glu_f,4)
a1[,paste0(varlist[15:25],'_q')] <- lapply(a1[,varlist[15:25]], function(x) {
  x <- factor(quantileCut(x,4),labels = c('Q1','Q2','Q3','Q4'))
  x <- as.character(x)
  x[is.na(x)] <- 'Missing'
  if(length(table(x))==4){
    x <- factor(x, levels = c('Q1','Q2','Q3','Q4'))
  } else {
    x <- factor(x, levels = c('Q1','Q2','Q3','Q4','Missing'))
  }
  x
})

# cut blood pressure
a1$systolic_2 <- cut(a1$systolic,c(0,140,Inf),right=F,labels = c('Normal','Abnormal'))
a1$systolic_2 <- as.character(a1$systolic_2)
a1$systolic_2[is.na(a1$systolic_2)] <- 'Missing'
a1$systolic_2 <- factor(a1$systolic_2,levels = c('Normal','Abnormal','Missing'))

a1$diastolic_2 <- cut(a1$diastolic,c(0,90,Inf),right=F,labels = c('Normal','Abnormal'))
a1$diastolic_2 <- as.character(a1$diastolic_2)
a1$diastolic_2[is.na(a1$diastolic_2)] <- 'Missing'
a1$diastolic_2 <- factor(a1$diastolic_2,levels = c('Normal','Abnormal','Missing'))








################################################################################
# Select patients by excluding prevalent CVD for survival analysis of true incidence

# all sample
a1_ihd <- which(a1$ischemic_heart_disease_grp_prev==0)
a1_hf <- which(a1$heart_failure_grp_prev==0)
a1_cm <- which(a1$cardiomyopathy_grp_prev==0)
a1_stroke <- which(a1$stroke_grp_prev==0)
a1_combo <- which(a1$cvdcombo_grp_prev==0)
a1_arrhythmia <- which(a1$arrhythmia_grp_prev==0)
a1_cardiac <- which(a1$cardiac_arrest_grp_prev==0)
a1_carotid <- which(a1$carotid_disease_grp_prev==0)
a1_myocarditis <- which(a1$myocarditis_pericarditis_grp_prev==0)
a1_tia <- which(a1$tia_grp_prev==0)
a1_valvular <- which(a1$valvular_disease_grp_prev==0)
a1_dvt <- which(a1$venous_thromboembolic_disease_grp_prev==0)
a1_diab <- which(a1$diab_bl=='No')
a1_htn <- which(a1$htn_bl=='No')
a1_dyslipid <- which(a1$dyslipid_bl=='No')
a1_rfcombo <- which(a1$diab_bl=='No' & a1$htn_bl=='No' & a1$dyslipid_bl=='No')

# receive chemo
id_chemo <- a1$num7_studyid[which(a1$chemo_yn=='Yes')]
a1_chemo <- which(a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo)
a1_ihd_chemo <- which(a1$ischemic_heart_disease_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_hf_chemo <- which(a1$heart_failure_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_cm_chemo <- which(a1$cardiomyopathy_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_stroke_chemo <- which(a1$stroke_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_combo_chemo <- which(a1$cvdcombo_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_arrhythmia_chemo <- which(a1$arrhythmia_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_cardiac_chemo <- which(a1$cardiac_arrest_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_carotid_chemo <- which(a1$carotid_disease_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_myocarditis_chemo <- which(a1$myocarditis_pericarditis_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_tia_chemo <- which(a1$tia_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_valvular_chemo <- which(a1$valvular_disease_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_dvt_chemo <- which(a1$venous_thromboembolic_disease_grp_prev==0 & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_diab_chemo <- which(a1$diab_bl=='No' & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_htn_chemo <- which(a1$htn_bl=='No' & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_dyslipid_chemo <- which(a1$dyslipid_bl=='No' & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))
a1_rfcombo_chemo <- which(a1$diab_bl=='No' & a1$htn_bl=='No' & a1$dyslipid_bl=='No' & (a1$num7_studyid %in% id_chemo | a1$num7_caseid %in% id_chemo))

# receive HT
id_horm <- a1$num7_studyid[which(a1$horm_yn=='Yes')]
a1_horm <- which(a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm)
a1_ihd_horm <- which(a1$ischemic_heart_disease_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_hf_horm <- which(a1$heart_failure_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_cm_horm <- which(a1$cardiomyopathy_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_stroke_horm <- which(a1$stroke_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_combo_horm <- which(a1$cvdcombo_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_arrhythmia_horm <- which(a1$arrhythmia_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_cardiac_horm <- which(a1$cardiac_arrest_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_carotid_horm <- which(a1$carotid_disease_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_myocarditis_horm <- which(a1$myocarditis_pericarditis_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_tia_horm <- which(a1$tia_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_valvular_horm <- which(a1$valvular_disease_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_dvt_horm <- which(a1$venous_thromboembolic_disease_grp_prev==0 & (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_diab_horm <- which(a1$diab_bl=='No'& (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_htn_horm <- which(a1$htn_bl=='No'& (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_dyslipid_horm <- which(a1$dyslipid_bl=='No'& (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))
a1_rfcombo_horm <- which(a1$diab_bl=='No' & a1$htn_bl=='No' & a1$dyslipid_bl=='No'& (a1$num7_studyid %in% id_horm | a1$num7_caseid %in% id_horm))

# receive RT
id_rad <- a1$num7_studyid[which(a1$rad_yn=='Yes')]
a1_rad <- which(a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad)
a1_ihd_rad <- which(a1$ischemic_heart_disease_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_hf_rad <- which(a1$heart_failure_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_cm_rad <- which(a1$cardiomyopathy_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_stroke_rad <- which(a1$stroke_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_combo_rad <- which(a1$cvdcombo_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_arrhythmia_rad <- which(a1$arrhythmia_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_cardiac_rad <- which(a1$cardiac_arrest_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_carotid_rad <- which(a1$carotid_disease_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_myocarditis_rad <- which(a1$myocarditis_pericarditis_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_tia_rad <- which(a1$tia_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_valvular_rad <- which(a1$valvular_disease_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_dvt_rad <- which(a1$venous_thromboembolic_disease_grp_prev==0 & (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_diab_rad <- which(a1$diab_bl=='No'& (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_htn_rad <- which(a1$htn_bl=='No'& (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_dyslipid_rad <- which(a1$dyslipid_bl=='No'& (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))
a1_rfcombo_rad <- which(a1$diab_bl=='No' & a1$htn_bl=='No' & a1$dyslipid_bl=='No'& (a1$num7_studyid %in% id_rad | a1$num7_caseid %in% id_rad))

# left side BC, receive RT
id_rad_l <- a1$num7_studyid[which(a1$rad_yn=='Yes' & a1$laterality==2)]
a1_rad_l <- which(a1$num7_studyid %in% id_rad_l | a1$num7_caseid %in% id_rad_l)
a1_ihd_rad_l <- which(a1$ischemic_heart_disease_grp_prev==0 & (a1$num7_studyid %in% id_rad_l | a1$num7_caseid %in% id_rad_l))
a1_hf_rad_l <- which(a1$heart_failure_grp_prev==0 & (a1$num7_studyid %in% id_rad_l | a1$num7_caseid %in% id_rad_l))
a1_cm_rad_l <- which(a1$cardiomyopathy_grp_prev==0 & (a1$num7_studyid %in% id_rad_l | a1$num7_caseid %in% id_rad_l))
a1_stroke_rad_l <- which(a1$stroke_grp_prev==0 & (a1$num7_studyid %in% id_rad_l | a1$num7_caseid %in% id_rad_l))
a1_combo_rad_l <- which(a1$cvdcombo_grp_prev==0 & (a1$num7_studyid %in% id_rad_l | a1$num7_caseid %in% id_rad_l))


# right side BC, receive RT
id_rad_r <- a1$num7_studyid[which(a1$rad_yn=='Yes' & a1$laterality==1)]
a1_rad_r <- which(a1$num7_studyid %in% id_rad_r | a1$num7_caseid %in% id_rad_r)
a1_ihd_rad_r <- which(a1$ischemic_heart_disease_grp_prev==0 & (a1$num7_studyid %in% id_rad_r | a1$num7_caseid %in% id_rad_r))
a1_hf_rad_r <- which(a1$heart_failure_grp_prev==0 & (a1$num7_studyid %in% id_rad_r | a1$num7_caseid %in% id_rad_r))
a1_cm_rad_r <- which(a1$cardiomyopathy_grp_prev==0 & (a1$num7_studyid %in% id_rad_r | a1$num7_caseid %in% id_rad_r))
a1_stroke_rad_r <- which(a1$stroke_grp_prev==0 & (a1$num7_studyid %in% id_rad_r | a1$num7_caseid %in% id_rad_r))
a1_combo_rad_r <- which(a1$cvdcombo_grp_prev==0 & (a1$num7_studyid %in% id_rad_r | a1$num7_caseid %in% id_rad_r))


