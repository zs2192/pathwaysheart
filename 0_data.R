
# Pathways Heart Study - Aim 1 data cleaning
# Zaixing Shi, 7/20/2019

library(tidyverse)
library(haven)
library(ggplot2)
library(dplyr)
library(plyr)
library(naniar)
library(reshape2)
library(doParallel)

# setup parallel computing
nodes <- detectCores()
cl <- makeCluster(nodes)
registerDoParallel(cl)

################################################################################
# Load data

# import cases data received on 2019-06-27
pathpref = 'Q:/HGREENLEE/Data Working/Pathways Heart Study/Aim 1 cases/2019-06-27/'
cases_new <- as.data.frame(read_sas(paste0(pathpref,'cases_final_27mar19.sas7bdat')))
cases_tumor <- as.data.frame(read_sas(paste0(pathpref,'cases_tumor_char_26jun19.sas7bdat')))

# import controls data received on 2019-07-03
pathpref = 'Q:/HGREENLEE/Data Working/Pathways Heart Study/Aim 1 cases/2019-07-03/'
controls1 <- as.data.frame(read_sas(paste0(pathpref,'controls_group1.sas7bdat')))
controls2 <- as.data.frame(read_sas(paste0(pathpref,'controls_group1.sas7bdat')))

## risk factor data received on 4/17/2019
pathpref = 'Q:/HGREENLEE/Data Working/Pathways Heart Study/Aim 1 cases/2019-04-17/'
bmi <- as.data.frame(read_sas(paste0(pathpref,'bmi.sas7bdat')))
bp <- as.data.frame(read_sas(paste0(pathpref,'bp.sas7bdat')))
lipid <- as.data.frame(read_sas(paste0(pathpref,'dyslipidemia.sas7bdat')))
diab <- as.data.frame(read_sas(paste0(pathpref,'diabetes.sas7bdat')))
smok <- as.data.frame(read_sas(paste0(pathpref,'smoking.sas7bdat')))
smok6 <- as.data.frame(read_sas(paste0(pathpref,'smoking_6months.sas7bdat')))

## menopause, parity, census data received on 4/17/2019
menop <- as.data.frame(read_sas(paste0(pathpref,'menopause.sas7bdat')))
parity <- as.data.frame(read_sas(paste0(pathpref,'parity.sas7bdat')))
census <- as.data.frame(read_sas(paste0(pathpref,'census.sas7bdat')))

# import CVD outcome data recevied on 6/27/2019
pathpref = 'Q:/HGREENLEE/Data Working/Pathways Heart Study/Aim 1 cases/2019-06-27/'
cvd <- as.data.frame(read_sas(paste0(pathpref,'cvd_events.sas7bdat')))

# import censoring events received on 6/27/2019
censor <- as.data.frame(read_sas(paste0(pathpref,'censoring_27jun19.sas7bdat')))

# import updated lab data received on 6/27/2019 
labs <- as.data.frame(read_sas(paste0(pathpref,'labs.sas7bdat')))






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
names(cases_new) <- tolower(names(cases_new))
names(cases_tumor) <- tolower(names(cases_tumor))
names(controls1) <- tolower(names(controls1))

# combine case and tumor characteristics data
cases <- merge(cases_new[,c("cvd_studyid","gender","raceethn1","birth_year",
                            "n_controls_group1","n_controls_group2","cops2")], 
               cases_tumor, by="cvd_studyid")

# create another surgery date var for cases and controls for risk factor data selection
cases$rf_date <- cases$surg_date
## fill patients without surg_date with dxdate plus median days between dx to surg
cases$rf_date[which(is.na(cases$rf_date))] <- cases$dxdate[which(is.na(cases$rf_date))]+
                                              as.numeric(median(cases$surg_date-cases$dxdate, na.rm=T))

controls1$rf_date <- cases$rf_date[match(controls1$case_id,cases$cvd_studyid)]

# stack cases and controls 1, matched on age and race/ethnicity
# reorder controls1 variable to be the same as cases
controls1$group <- 'Control'
controls1$dxage <- controls1$index_age
controls1[,names(cases)[!(names(cases) %in% names(controls1))]] <- NA

cases[,names(controls1)[!(names(controls1) %in% names(cases))]] <- NA
cases$group <- 'Case'
controls1 <- controls1[,names(cases)]

# combine cases and control data, n=89644
all <- rbind(cases, controls1)





################################################################################
# Select labs data
# Date range: closest measure within 36 mos prior to index date and before surgery

# impute missing surgery dates
cases$dx2surgdays <- as.numeric(cases$surg_date - cases$index_date) 
cases$surg_date[which(is.na(cases$surg_date) & cases$surg_sum>0)] <- 
  cases$index_date[which(is.na(cases$surg_date) & cases$surg_sum>0)] + 
  median(cases$dx2surgdays,na.rm=T)


# BMI
names(bmi) <- tolower(names(bmi))
bmi <- merge(bmi, all[,c("cvd_studyid","dxdate",'rf_date')], by='cvd_studyid')
bmi1 <- bmi[which(bmi$measure_date >= bmi$index_date-365*3 & bmi$measure_date <= bmi$rf_date),]
bmi1$datediff <- abs(as.numeric(bmi1$measure_date - bmi1$index_date))
bmi1 <- ddply(bmi1, .(cvd_studyid), function(d){
  d <- d[which(d$datediff==min(d$datediff)),]
  d
},.parallel=TRUE)

## clean up BMI
bmi1$bmi_num <- gsub('<|>|>=','',bmi1$bmi)
bmi1$bmi_num <- as.numeric(gsub('-.*','',bmi1$bmi_num))

## aggregate multiple values within same day
bmi2 <- aggregate(x = bmi1$bmi_num, 
                  by = list(bmi1$cvd_studyid,bmi1$datediff), 
                  mean, na.rm=T)
names(bmi2) <- c('cvd_studyid','bmi_datediff','bmi')




# BP
names(bp) <- tolower(names(bp))
bp <- merge(bp, all[,c("cvd_studyid","dxdate",'rf_date')], by='cvd_studyid')
bp1 <- bp[which(bp$measure_date >= bp$index_date-365*2 & bp$measure_date <= bp$rf_date),]
bp1$datediff <- abs(as.numeric(bp1$measure_date - bp1$index_date))
bp1 <- ddply(bp1, .(cvd_studyid), function(d){
  d <- d[which(d$datediff==min(d$datediff)),]
  d
},.parallel=TRUE)

## clean up bp
bp1$hypertension <- factor(bp1$hypertension,
                           labels = c('Normal','Elevated','Stage 1','Stage 2','Undefined'))

## aggregate multiple values within same day
systol <- aggregate(x = bp1$systolic, 
                  by = list(bp1$cvd_studyid,bp1$datediff), 
                  mean, na.rm=T)
diastol <- aggregate(x = bp1$diastolic, 
                   by = list(bp1$cvd_studyid,bp1$datediff), 
                   mean, na.rm=T)

bp2 <- merge(systol, diastol, by=c("Group.1","Group.2"))
names(bp2) <- c('cvd_studyid','bp_datediff','systolic','diastolic')



# Labs
labs$result <- as.numeric(labs$RESULT_C)

# merge surgery date
names(labs) <- tolower(names(labs))
labs <- merge(labs, all[,c("cvd_studyid","dxdate",'rf_date')], by='cvd_studyid')
labs1 <- labs[which(labs$test_dt >= labs$index_date-365*3 & labs$test_dt <= labs$rf_date),]
labs1$datediff <- abs(as.numeric(labs1$test_dt - labs1$index_date))
system.time(labs1 <- ddply(labs1, .(cvd_studyid,test_type), function(d){
  d <- d[which(d$datediff==min(d$datediff)),]
  d
},.parallel=TRUE))

labs2 <- dcast(data=labs1,cvd_studyid+index_date+control_group~test_type,
                value.var = 'result', mean, na.rm=T)

# save a copy of labs2 data
write_csv(labs2,'labs2.csv')

# Smoking - data pending
#smok$datediff <- as.numeric(smok$CONTACT_DATE - smok$index_date)

#smok1 <- smok[which(smok$datediff %in% -365:0),]

#smok1 <- ddply(smok1, .(CVD_studyid), function(d){
#  d <- d[order(d$datediff, decreasing = T),]
#  d[1,]
#})




################################################################################
# Clean CVD outcome data

## consolidate duplicated cvd group
cvd$CVD_condition[cvd$CVD_condition=='Percutaneous transluminal coronary angioplasty status'] <- 'Percutaneous transluminal coronary angioplasty'

## make long to wide format, keep 1st occurence of each condition only
### condition groups
cvd2 <- dcast(data=cvd, CVD_STUDYID+index_date~CVD_group, value.var = 'adate', min,na.rm=T)

cvd2[,-1:-2] <- lapply(cvd2[,-1:-2], function(x) as.Date(x, origin='1970-01-01'))
cvd2[,-1:-2] <- lapply(cvd2[,-1:-2], function(x) as.Date(as.character(x)))

names(cvd2) <- gsub('/| ','_',names(cvd2))
cvd_grp <- paste0(names(cvd2)[-1:-2],'_grp')
names(cvd2)[-1:-2] <- paste0(cvd_grp,'_dt')

cvd2[,cvd_grp] <- lapply(cvd2[,-1:-2], function(x) {
  y <- ifelse(!is.na(x),1,0)
  y[y==1 & x>cvd2$index_date] <- 2
  y <- factor(y, levels=0:2, labels = c('No','Prevalent','Incident'))
  y
  })

cvd2list <- lapply(cvd_grp, function(x){
  dates <- eval(parse(text=paste0("dcast(data=cvd2, CVD_STUDYID~",x,", value.var = '",x,"_dt')")))
  prev <- ifelse(cvd2[,x]=='Prevalent',1,0)
  inc <- ifelse(cvd2[,x]=='Incident',1,0)
  sum <- data.frame(cbind(prev,inc,dates))
  sum <- sum[,c('CVD_STUDYID','prev','inc','Incident')]
  sum$Incident <- as.Date(sum$Incident, origin='1970-01-01')
  names(sum) <- c('CVD_STUDYID',paste0(x,c('_prev','_inc','_incdt')))
  sum
})

cvd2list[[10]] <- cvd2
cvd2_final <- Reduce(function(...) merge(...,by='CVD_STUDYID', all.x=T), cvd2list)




### condition
cvd3 <- dcast(data=cvd, CVD_STUDYID+index_date~CVD_condition, value.var = 'adate', min,na.rm=T)

cvd3[,-1:-2] <- lapply(cvd3[,-1:-2], function(x) as.Date(x, origin='1970-01-01'))
cvd3[,-1:-2] <- lapply(cvd3[,-1:-2], function(x) as.Date(as.character(x)))

names(cvd3) <- gsub('/| ','_',names(cvd3))
names(cvd3) <- gsub('&_|\\(|\\)','',names(cvd3))
names(cvd3)[18] <- "Percutaneous_coronary_intervention_PCI"

cvd_cond <- names(cvd3)[-1:-2]
names(cvd3)[-1:-2] <- paste0(names(cvd3)[-1:-2],'_dt')
cvd3[,cvd_cond] <- lapply(cvd3[,-1:-2],  function(x) {
  y <- ifelse(!is.na(x),1,0)
  y[y==1 & x>cvd2$index_date] <- 2
  y <- factor(y, levels=0:2, labels = c('No','Prevalent','Incident'))
  y
})

cvd3list <- lapply(cvd_cond, function(x){
  dates <- eval(parse(text=paste0("dcast(data=cvd3, CVD_STUDYID~",x,", value.var = '",x,"_dt')")))
  prev <- ifelse(cvd3[,x]=='Prevalent',1,0)
  inc <- ifelse(cvd3[,x]=='Incident',1,0)
  sum <- data.frame(cbind(prev,inc,dates))
  sum <- sum[,c('CVD_STUDYID','prev','inc','Incident')]
  sum$Incident <- as.Date(sum$Incident, origin='1970-01-01')
  names(sum) <- c('CVD_STUDYID',paste0(x,c('_prev','_inc','_incdt')))
  sum
})

cvd3list[[24]] <- cvd3
cvd3_final <- Reduce(function(...) merge(...,by='CVD_STUDYID', all.x=T), cvd3list)



# export condition list
write_csv(cvd_grp,'cvd_grp.csv')
write_csv(cvd_cond,'cvd_cond.csv')

# save a copy of the cvd data
write_csv(cvd2_final,'cvd2_final.csv')
write_csv(cvd3_final,'cvd3_final.csv')





################################################################################
# Merge all data

## combine all data into a list
datalist <- list(all[,-grep("deathdt|enr_start|enr_end",names(all))], 
                 bmi2,bp2,
                 labs2[,c("cvd_studyid","GLU_F","GTT75_PRE","HDL","HGBA1C","LDL_CLC_NS","TOT_CHOLES","TRIGL_NS")],
                 diab[,c(1,4:10)],lipid[,2:4], 
                 smok1[,c(2,5:8)],menop[,c(1,3)],parity[,c(1,5,6)],
                 census[,c(2,5:48)],
                 cvd2_final[,-grep('index_date',names(cvd2_final))],
                 cvd3_final[,-grep('index_date',names(cvd3_final))],
                 censor[,-c(2,7)])

## make all var names lower case
datalist <- lapply(datalist,function(x){
  names(x) <- tolower(names(x))
  x
})

## merge all data by CVD_studyid
a1 <- Reduce(function(...) merge(...,by='cvd_studyid', all.x=T), datalist)

## set all NA for cvd outcomes as 0
a1[,tolower(c(cvd_grp,cvd_cond))] <- lapply(a1[,tolower(c(cvd_grp,cvd_cond))], function(x) {
  x[is.na(x)] <- 'No'
  x
})

a1[,grep('inc$|prev$',names(a1))] <- 
  lapply(a1[,grep('inc$|prev$',names(a1))], function(x) {
    x[is.na(x)] <- 0
    x
  })

# save a copy of the merged data
write_csv(a1,'a1.csv')






################################################################################
# Create some more variables...

# recode categorical variables
# age groups
a1$agegrp <- cut(a1$dxage, c(0, 40, 50, 60, 70, 80, 101), right=F,
                  labels = c('<40 yo','40-49','50-59','60-69','70-79','80+ yo'))

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

# convert dxdate to date
a1$dxdate <- as.Date(a1$dxdate, origin = '1970-01-01')

# enrolled in Pathways
a1$enrolled <- factor(a1$enrolled, levels = c(0,1), labels=c('No','Yes'))

# create index year
a1$index_yr <- format(a1$index_date,"%Y")
a1$index_yr[a1$index_yr %in% c('2005','2006','2007')] <- '2005-2007' 
a1$index_yr[a1$index_yr %in% c('2008','2009','2010')] <- '2008-2010' 
a1$index_yr[a1$index_yr %in% c('2011','2012','2013')] <- '2011-2013' 
a1$index_yr <- factor(a1$index_yr)

# create enrollment length
a1$enr_len <- as.numeric((a1$enr_end - a1$enr_start)/30)


# BMI category
a1$bmicat <- cut(a1$bmi, c(0,18.5,25,30,35,Inf),right=F,
                 labels = c('Underweight','Normal','Overweight','Obese I','Obese II+'))
a1$bmicat <- as.character(a1$bmicat)
a1$bmicat[is.na(a1$bmicat)] <- 'Unknown'
a1$bmicat <- factor(a1$bmicat,
                    levels = c('Underweight','Normal','Overweight','Obese I','Obese II+','Unknown'))

# smoking status
a1$smok1 <- a1$tobacco_use
a1$smok1[a1$smok %in% c('P','Q')] <- 'Quited/past/former'
a1$smok1[a1$smok %in% c('Y')] <- 'Current smoker'
a1$smok1[a1$smok %in% c('N')] <- 'Never smoker'
a1$smok1[a1$smok %in% c('U')] <- 'Unknown'
a1$smok1[is.na(a1$smok)] <- 'Unknown'

a1$smok1 <- factor(a1$smok1, levels=c('Never smoker','Current smoker','Quited/past/former','Unknown'))

# menopausal status
a1$menop <- factor(a1$bl_meno_status)
a1$menop[is.na(a1$menop) & a1$dxage>51] <- '1'
a1$menop[is.na(a1$menop) & a1$dxage<=51] <- '0'

# gravidity
a1$gravid <- cut(as.numeric(a1$ob_gravidity), c(0,1,2,3,Inf),right=F,
                 labels = c('0','1','2','3+'))

a1$parity <- cut(as.numeric(a1$ob_parity), c(0,1,2,3,Inf),right=F,
                 labels = c('0','1','2','3+'))

# diabetes
a1[,c("inpatient_flg","outpatient_flg","pharmacy_flg","a1c_flg",
             "fgrg_flg","a1cfgrg_flg")] <- lapply(a1[,c("inpatient_flg","outpatient_flg","pharmacy_flg","a1c_flg",
                                                        "fgrg_flg","a1cfgrg_flg")],
                                                  function(x){
                                                    x[is.na(x)] <- 2
                                                    x <- factor(x, levels=c(1,0,2),labels = c('Yes','No','Unknown'))
                                                    x
                                                  })
a1$diab <- factor(ifelse(a1$inpatient_flg=='Yes'|a1$outpatient_flg=='Yes'|a1$pharmacy_flg=='Yes'|a1$a1c_flg=='Yes'|a1$fgrg_flg=='Yes'|a1$a1cfgrg_flg=='Yes',
                  'Yes','No'))

# dyslipidemia
a1$dyslipidemia[is.na(a1$dyslipidemia)] <- 0
a1$dyslipidemia <- factor(a1$dyslipidemia,levels=c(1,0),
                          labels = c('Yes','No/Unknown'))

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

# combined event of ischemic heart disease, stroke/TIA, cardiomyopathy/heart failure
a1$cvdcombo_grp_inc <- ifelse(a1$ischemic_heart_disease_grp_inc==1 | 
                            a1$stroke_tia_grp_inc==1 | 
                            a1$cardiomyopathy_heart_failure_grp_inc==1, 1, 0)


a1$cvdcombo_grp_incdt <- as.Date(apply(a1[,c("ischemic_heart_disease_grp_incdt",
                                             "stroke_tia_grp_incdt",
                                             "cardiomyopathy_heart_failure_grp_incdt")],
                                       1,min,na.rm=T))


# define censoring time
a1$censor_dt <- as.Date(apply(a1[,c("death_date","enr_end","end_of_study")],1,min,na.rm=T))


# calculate time of follow up
a1$ischemic_heart_disease_grp_inc_fu <- as.numeric(a1$ischemic_heart_disease_grp_incdt-a1$index_date)
a1$ischemic_heart_disease_grp_inc_fu[which(a1$ischemic_heart_disease_grp_inc==0)] <- 
  as.numeric(a1$censor_dt[which(a1$ischemic_heart_disease_grp_inc==0)] - 
               a1$index_date[which(a1$ischemic_heart_disease_grp_inc==0)])


a1$stroke_tia_grp_inc_fu <- as.numeric(a1$stroke_tia_grp_incdt-a1$index_date)
a1$stroke_tia_grp_inc_fu[which(a1$stroke_tia_grp_inc==0)] <- 
  as.numeric(a1$censor_dt[which(a1$stroke_tia_grp_inc==0)] - 
               a1$index_date[which(a1$stroke_tia_grp_inc==0)])


a1$cardiomyopathy_heart_failure_grp_inc_fu <- as.numeric(a1$cardiomyopathy_heart_failure_grp_incdt-a1$index_date)
a1$cardiomyopathy_heart_failure_grp_inc_fu[which(a1$cardiomyopathy_heart_failure_grp_inc==0)] <- 
  as.numeric(a1$censor_dt[which(a1$cardiomyopathy_heart_failure_grp_inc==0)] - 
               a1$index_date[which(a1$cardiomyopathy_heart_failure_grp_inc==0)])


a1$cvdcombo_grp_inc_fu <- as.numeric(a1$cvdcombo_grp_incdt-a1$index_date)
a1$cvdcombo_grp_inc_fu[which(a1$cvdcombo_grp_inc==0)] <- 
  as.numeric(a1$censor_dt[which(a1$cvdcombo_grp_inc==0)] - 
               a1$index_date[which(a1$cvdcombo_grp_inc==0)])






################################################################################
# Data process for modeling

# relevel the group var, with control as referent
a1$group1 <- factor(a1$group, levels = c('Control','Case'))

# relevel bmi var
a1$bmicat1 <- factor(a1$bmicat, levels = c('Normal','Underweight','Overweight','Obese I','Obese II+','Unknown'))

library(lsr)
# cut continuous var into quantiles
a1$glu_q <- quantileCut(a1$glu_f,4)
a1[,paste0(varlist[22:36],'_q')] <- lapply(a1[,varlist[22:36]], function(x) {
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

# relevel dyslipidemia
a1$dyslipidemia2 <- factor(a1$dyslipidemia, levels=c('No/Unknown','Yes'))

################################################################################
# Create datasets with cases who received treatment and their matched controls

# pull IDs of cases receiving chemo, HT and RT
id_chemo <- cases$cvd_studyid[which(cases$chemo_yn==1)]
id_horm <- cases$cvd_studyid[which(cases$horm_yn==1)]
id_rad <- cases$cvd_studyid[which(cases$rad_yn==1)]

# select subsample based on these ids
a1chemo <- a1[which(a1$cvd_studyid %in% id_chemo | a1$case_id %in% id_chemo),]
a1horm <- a1[which(a1$cvd_studyid %in% id_horm | a1$case_id %in% id_horm),]
a1rad <- a1[which(a1$cvd_studyid %in% id_rad | a1$case_id %in% id_rad),]






