# CVD and BC data cleaning
# Zaixing Shi, 5/7/2019

library(tidyverse)
library(haven)
library(ggplot2)
library(dplyr)
library(plyr)
library(naniar)
library(reshape2)


################################################################################
# Load data

# import data received on 2018-06-20
pathpref='Q:/HGREENLEE/Data Working/Pathways Heart Study/Aim 1 cases/2018-06-20/'
cases <- as.data.frame(read_sas(paste0(pathpref,'Num0001_group1_cases.sas7bdat')))
controls <- as.data.frame(read_sas(paste0(pathpref,'Num0001_group1_controls.sas7bdat')))


# import data received on 2019-04-17
pathpref = 'Q:/HGREENLEE/Data Working/Pathways Heart Study/Aim 1 cases/2019-04-17/'
## updated cases and controls data
cases_new <- as.data.frame(read_sas(paste0(pathpref,'cases.sas7bdat')))
controls1 <- as.data.frame(read_sas(paste0(pathpref,'controls_group1.sas7bdat')))
controls2 <- as.data.frame(read_sas(paste0(pathpref,'controls_group1.sas7bdat')))

## risk factor data
bmi <- as.data.frame(read_sas(paste0(pathpref,'bmi.sas7bdat')))
bp <- as.data.frame(read_sas(paste0(pathpref,'bp.sas7bdat')))
labs <- as.data.frame(read_sas(paste0(pathpref,'labs.sas7bdat')))
lipid <- as.data.frame(read_sas(paste0(pathpref,'dyslipidemia.sas7bdat')))
diab <- as.data.frame(read_sas(paste0(pathpref,'diabetes.sas7bdat')))
smok <- as.data.frame(read_sas(paste0(pathpref,'smoking.sas7bdat')))
smok6 <- as.data.frame(read_sas(paste0(pathpref,'smoking_6months.sas7bdat')))

## menopause, parity, census
menop <- as.data.frame(read_sas(paste0(pathpref,'menopause.sas7bdat')))
parity <- as.data.frame(read_sas(paste0(pathpref,'parity.sas7bdat')))
census <- as.data.frame(read_sas(paste0(pathpref,'census.sas7bdat')))

# import CVD outcome data recevied on 5/1/2019
pathpref = 'Q:/HGREENLEE/Data Working/Pathways Heart Study/CVD outcome/2019-05-01/'
cvd <- as.data.frame(read_sas(paste0(pathpref,'cvd_outcomes_thruDec2018.sas7bdat')))







################################################################################
# Combine cases and controls data

# check missingness
png('case_missing.png', height=3.5, width = 2.5, unit='in', res=200)
gg_miss_var(cases, show_pct =TRUE) + labs(y = "% missing, Cases")
dev.off()

png('control_missing.png', height=3.5, width = 2.5, unit='in', res=200)
gg_miss_var(controls, show_pct =TRUE) + labs(y = "% missing, Controls")
dev.off()


# combine old and new cases data
cases$CVD_studyid <- cases$case_id 
cases <- merge(cases, cases_new, by=names(cases_new)[names(cases_new) %in% names(cases)])

# stack cases and controls 1, matched on age and race/ethnicity

# change all variables to lower case 
names(cases) <- tolower(names(cases))
names(controls1) <- tolower(names(controls1))

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
# Date range: closest measure within 24 mos prior to index date and before surgery

# BMI
bmi$datediff <- as.numeric(bmi$MEASURE_DATE - bmi$index_date)
bmi1 <- bmi[which(bmi$datediff %in% -365:0),]
bmi1 <- ddply(bmi1, .(CVD_studyid), function(d){
  d <- d[which(d$datediff==max(d$datediff)),]
  d
})

## clean up BMI
bmi1$bmi_num <- gsub('<|>|>=','',bmi1$BMI)
bmi1$bmi_num <- as.numeric(gsub('-.*','',bmi1$bmi_num))

## aggregate multiple values within same day
bmi1 <- aggregate(x = bmi1$bmi_num, 
                  by = list(bmi1$CVD_studyid,bmi1$index_date,
                            bmi1$MEASURE_DATE,bmi1$datediff), 
                  mean, na.rm=T)
names(bmi1) <- c('cvd_studyid','index_date','bmi_measure_date','bmi_datediff','bmi')




# BP
bp$datediff <- as.numeric(bp$MEASURE_DATE - bp$index_date)
bp1 <- bp[which(bp$datediff %in% -365:0),]
bp1 <- ddply(bp1, .(CVD_studyid), function(d){
  d <- d[which(d$datediff==max(d$datediff)),]
  d
})

## clean up bp
bp1$hypertension <- factor(bp1$hypertension,
                           labels = c('Normal','Elevated','Stage 1','Stage 2','Undefined'))

## aggregate multiple values within same day
systol <- aggregate(x = bp1$SYSTOLIC, 
                  by = list(bp1$CVD_studyid,bp1$index_date,
                            bp1$MEASURE_DATE,bp1$datediff), 
                  mean, na.rm=T)
diastol <- aggregate(x = bp1$DIASTOLIC, 
                   by = list(bp1$CVD_studyid,bp1$index_date,
                             bp1$MEASURE_DATE,bp1$datediff), 
                   mean, na.rm=T)

bp1 <- merge(systol, diastol, by=c("Group.1","Group.2", "Group.3" ,"Group.4"))
names(bp1) <- c('cvd_studyid','index_date','bp_measure_date','bp_datediff','systolic','diastolic')



# Labs
labs$result <- as.numeric(labs$RESULT_C)
labs$datediff <- as.numeric(labs$test_dt - labs$index_date)

labs1 <- labs[which(labs$datediff %in% -365:0),]

labs2 <- ddply(labs1, .(CVD_studyid,TEST_TYPE), function(d){
  d <- d[which(d$datediff==max(d$datediff)),]
  d
})

labs2 <- dcast(data=labs2,CVD_studyid+index_date+control_group+caco~TEST_TYPE,
                value.var = 'result', mean, na.rm=T)



# Smoking
smok$datediff <- as.numeric(smok$CONTACT_DATE - smok$index_date)

smok1 <- smok[which(smok$datediff %in% -365:0),]

smok1 <- ddply(smok1, .(CVD_studyid), function(d){
  d <- d[order(d$datediff, decreasing = T),]
  d[1,]
})




################################################################################
# Clean CVD outcome data

## keep events after index date

cvd1 <- cvd[which(cvd$adate > cvd$index_date),]


## make long to wide format, keep 1st occurence of each condition only
### condition groups
cvd2 <- dcast(data=cvd1, CVD_STUDYID~group, value.var = 'adate', min,na.rm=T)

cvd2[,-1] <- lapply(cvd2[,-1], function(x) as.Date(x, origin='1970-01-01'))
cvd2[,-1] <- lapply(cvd2[,-1], function(x) as.Date(as.character(x)))

names(cvd2) <- gsub('/| ','_',names(cvd2))
cvd_grp <- names(cvd2)[-1]
names(cvd2)[-1] <- paste0(names(cvd2)[-1],'_dt')
cvd2[,cvd_grp] <- lapply(cvd2[,-1], function(x) ifelse(!is.na(x),1,0))


### condition
cvd3 <- dcast(data=cvd1, CVD_STUDYID~CONDITION, value.var = 'adate', min,na.rm=T)

cvd3[,-1] <- lapply(cvd3[,-1], function(x) as.Date(x, origin='1970-01-01'))
cvd3[,-1] <- lapply(cvd3[,-1], function(x) as.Date(as.character(x)))

names(cvd3) <- gsub('/| ','_',names(cvd3))
cvd_cond <- names(cvd3)[-1]
names(cvd3)[-1] <- paste0(names(cvd3)[-1],'_dt')
cvd3[,cvd_cond] <- lapply(cvd3[,-1], function(x) ifelse(!is.na(x),1,0))



# export condition list
write.csv(cvd_grp,'cvd_grp.csv')
write.csv(cvd_cond,'cvd_cond.csv')


################################################################################
# Merge all data


## combine all data into a list
datalist <- list(all, 
                 bmi1[,c("cvd_studyid","bmi_measure_date", "bmi_datediff", "bmi")],
                 bp1[,c("cvd_studyid","bp_measure_date", "bp_datediff","systolic","diastolic")],
                 labs2[,c("CVD_studyid","GLU_F","GTT75_PRE","HDL","HGBA1C","LDL_CLC_NS","TOT_CHOLES","TRIGL_NS" )],
                 diab[,c(1,4:10)],lipid[,2:4], 
                 smok1[,c(2,5:8)],menop[,c(1,3)],parity[,c(1,5,6)],
                 census[,c(2,5:48)],
                 cvd2)

## make all var names lower case
datalist <- lapply(datalist,function(x){
  names(x) <- tolower(names(x))
  x
})


## merge all data by CVD_studyid
a1 <- Reduce(function(...) merge(...,by='cvd_studyid', all.x=T), datalist)

## set all NA for cvd outcomes as 0
a1[,124:133] <- lapply(a1[,124:133], function(x) {
  x[is.na(x)] <- 0
  x
})





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


# dyslipidemia
a1$dyslipidemia.x[is.na(a1$dyslipidemia.x)] <- 0
a1$dyslipidemia.x <- factor(a1$dyslipidemia.x,levels=c(1,0),
                          labels = c('Yes','No/Unknown'))


# household income
a1$hhincome1 <- rowSums(a1[,c(96:104)], na.rm = T)
a1$hhincome2 <- rowSums(a1[,c(105:107)], na.rm = T)
a1$hhincome3 <- rowSums(a1[,c(108:111)], na.rm = T)

# education
a1$edu1 <- rowSums(a1[,c(70:72)], na.rm = T)
a1$edu2 <- rowSums(a1[,c(73:74)], na.rm = T)
a1$edu3 <- a1$education6
a1$edu4 <- rowSums(a1[,c(76:77)], na.rm = T)


# examine data distribution


par(mfrow=c(4,6))
lapply(5:114, function(x){
  if(class(a1[,x]) %in% c('numeric')){
    plot(density(a1[,x], na.rm=T), main = names(a1)[x])
  } 
  if(class(a1[,x]) %in% c('factor')){
    plot(a1[,x], main = names(a1)[x])
  } 
  if(class(a1[,x]) %in% c('Date')){
    hist(a1[,x], breaks = 'year', main = names(a1)[x])
  } 
})

dev.off()











