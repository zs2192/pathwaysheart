# CVD and BC exploratory data analysis
# Zaixing Shi, 5/7/2019


library(tidyverse)
library(haven)
library(ggplot2)
library(dplyr)
library(plyr)
library(naniar)
library(reshape2)
library(survival)
library(survminer)


################################################################################
# compare characteristics between cases and controls

tab1 <- function(d,x){
  if(class(d[,x]) %in% c('numeric')){
    n <- c(paste0(length(which(!is.na(d[,x]))),' (',
                  round(100*length(which(is.na(d[,x])))/nrow(d)),'%)'),
           tapply(d[,x], d$group, function(y) 
             paste0(length(which(!is.na(y))),' (',
                    round(100*length(which(is.na(y)))/length(y)),'%)')))
    m <- c(mean(d[,x], na.rm=T),
           tapply(d[,x], d$group,mean, na.rm=T))
    sd <- c(sd(d[,x], na.rm=T),
            tapply(d[,x], d$group,sd, na.rm=T))
    p <- t.test(d[,x]~d$group)$p.value
    sum <- c(x,n[1],m[1],sd[1],'',n[2],m[2],sd[2],'',n[3],m[3],sd[3],'',p)
  } 
  if(class(d[,x]) %in% c('factor')){
    n <- cbind(table(d[,x]),table(d[,x], d$group))
    pct <- cbind(prop.table(table(d[,x])),
                 prop.table(table(d[,x], d$group),2))
    p <- chisq.test(table(d[,x], d$group))$p.value
    sum <- cbind(x,'',n[,1],pct[,1],'','',n[,2],pct[,2],'','',n[,3],pct[,3],'',p)
  } 
  sum
}
# variables to analyze
varlist <- c('dxage','enr_len','cops2','bmi','agegrp','raceethn1','index_yr',
             'bmicat','smok1','menop','gravid','parity','diab','dyslipidemia',
             "systolic" ,"diastolic","glu_f","hdl","hgba1c",
             "ldl_clc_ns","tot_choles","trigl_ns",
             'hhincome1','hhincome2',"hhincome3",'medhousincome',
             'houspoverty', "edu1","edu2","edu3","edu4")

# all sample
table1_all <- do.call(rbind,lapply(varlist,function(x) tab1(a1,x)))
table1_all <- data.frame(cbind(row.names(table1_all), table1_all))

# cases receiving chemo
table1_chemo <- do.call(rbind,lapply(varlist,function(x) tab1(a1chemo,x)))
table1_chemo <- data.frame(cbind(row.names(table1_chemo), table1_chemo))

# cases receiving hormonal therapy
table1_horm <- do.call(rbind,lapply(varlist,function(x) tab1(a1horm,x)))
table1_horm <- data.frame(cbind(row.names(table1_horm), table1_horm))

# cases receiving radiation
table1_rad <- do.call(rbind,lapply(varlist,function(x) tab1(a1rad,x)))
table1_rad <- data.frame(cbind(row.names(table1_rad), table1_rad))

# export tables
write_csv(table1_all,'table1_all.csv')
write_csv(table1_chemo,'table1_chemo.csv')
write_csv(table1_horm,'table1_horm.csv')
write_csv(table1_rad,'table1_rad.csv')








################################################################################
# Compare cvd incidence between cases and controls

tab2 <- function(d,x){
  d[,x] <- factor(d[,x],levels=c(0,1),labels = c('No','Yes'))
  n <- cbind(table(d[,x]),table(d[,x], d$group))
  pct <- cbind(prop.table(table(d[,x])),
               prop.table(table(d[,x], d$group),2))
  p <- chisq.test(table(d[,x], d$group))$p.value
  sum <- cbind(x,'',n[,1],pct[,1],'','',n[,2],pct[,2],'','',n[,3],pct[,3],'',p)
  sum
}

cvdvars <- grep('_inc$',names(a1),value = T)

## compare incidence between all cases and controls
table2_cvd <- do.call(rbind,lapply(cvdvars ,function(x) tab2(a1,x)))
table2_cvd <- data.frame(cbind(row.names(table2_cvd), table2_cvd))

## compare incidence between cases receiving chemo and controls
table2_cvd_chemo <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1chemo,x)))
table2_cvd_chemo <- data.frame(cbind(row.names(table2_cvd_chemo), table2_cvd_chemo))

## compare incidence between cases receiving horm and controls
table2_cvd_horm <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1horm,x)))
table2_cvd_horm <- data.frame(cbind(row.names(table2_cvd_horm), table2_cvd_horm))

## compare incidence between cases receiving rad and controls
table2_cvd_rad <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1rad,x)))
table2_cvd_rad <- data.frame(cbind(row.names(table2_cvd_rad), table2_cvd_rad))


write_csv(table2_cvd[table2_cvd$V1=='Yes',], 'table2_cvd.csv')
write_csv(table2_cvd_chemo[table2_cvd_chemo$V1=='Yes',], 'table2_cvd_chemo.csv')
write_csv(table2_cvd_horm[table2_cvd_horm$V1=='Yes',], 'table2_cvd_horm.csv')
write_csv(table2_cvd_rad[table2_cvd_rad$V1=='Yes',], 'table2_cvd_rad.csv')







################################################################################
# Compare cvd prevalence between cases and controls

cvdvars <- grep('_prev$',names(a1),value = T)

## compare prevalent CVD by age group
cvdp_age <- data.frame(t(sapply(cvdvars, function(x) 
  prop.table(table(factor(a1[,x],levels=0:1), a1$agegrp),2)[2,])))
write_csv(cvdp_age, 'cvd_prev_agegrp.csv')

## compare incidence between all cases and controls
table2_prev_cvd <- do.call(rbind,lapply(cvdvars ,function(x) tab2(a1,x)))
table2_prev_cvd <- data.frame(cbind(row.names(table2_prev_cvd), table2_prev_cvd))

## compare incidence between cases receiving chemo and controls
table2_prev_cvd_chemo <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1chemo,x)))
table2_prev_cvd_chemo <- data.frame(cbind(row.names(table2_prev_cvd_chemo), table2_prev_cvd_chemo))

## compare incidence between cases receiving horm and controls
table2_prev_cvd_horm <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1horm,x)))
table2_prev_cvd_horm <- data.frame(cbind(row.names(table2_prev_cvd_horm), table2_prev_cvd_horm))

## compare incidence between cases receiving rad and controls
table2_prev_cvd_rad <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1rad,x)))
table2_prev_cvd_rad <- data.frame(cbind(row.names(table2_prev_cvd_rad), table2_prev_cvd_rad))

# export table
write_csv(table2_prev_cvd[table2_prev_cvd$V1=='Yes',], 'table2_prev_cvd.csv')
write_csv(table2_prev_cvd_chemo[table2_prev_cvd_chemo$V1=='Yes',], 'table2_prev_cvd_chemo.csv')
write_csv(table2_prev_cvd_horm[table2_prev_cvd_horm$V1=='Yes',], 'table2_prev_cvd_horm.csv')
write_csv(table2_prev_cvd_rad[table2_prev_cvd_rad$V1=='Yes',], 'table2_prev_cvd_rad.csv')







################################################################################
# date buffer vs. available lab data

#bmi

bmi_n <- data.frame(t(sapply(-7:-90, function(x){
  ids <- unique(bmi$CVD_studyid[which(bmi$datediff %in% x:0)])
  c(x,prop.table(table(a1$cvd_studyid %in% ids)))
})))

bp_n <- data.frame(t(sapply(-7:-90, function(x){
  ids <- unique(bp$CVD_studyid[which(bp$datediff %in% x:0)])
  c(x,prop.table(table(a1$cvd_studyid %in% ids)))
})))

lab_n <- data.frame(t(sapply(-7:-90, function(x){
  ids <- unique(labs$CVD_studyid[which(labs$datediff %in% x:0)])
  c(x,prop.table(table(a1$cvd_studyid %in% ids)))
})))

sample_n <- rbind(bmi_n,bp_n,lab_n)
sample_n$sample <- factor(rep(c('BMI','Blood pressure','Labs'), each=84),
                          levels=c('BMI','Blood pressure','Labs'))

# plot the relationship

png('date vs cutoff.png',width = 5,height = 3.5,res=400, units = 'in')
ggplot(sample_n)+
  geom_line(aes(x=V1,y=100*TRUE.,color=sample),size=1)+
  scale_x_continuous(name = 'Cutoff days before index date',breaks = c(-7,-30,-60,-90))+
  scale_y_continuous(name = '% have data',breaks = c(0,10,20,30,40,50))+
  scale_color_discrete(name='Data')+
  theme_bw()
dev.off()



################################################################################
# data distribution

par(mfrow=c(2,5))
lapply(c('cops2','bmi',"systolic" ,"diastolic","glu_f","hdl","hgba1c",
             "ldl_clc_ns","tot_choles","trigl_ns"),
       function(x){
         min_value <- min(a1[,x], na.rm=T)
         max_value <- max(a1[,x], na.rm=T)
         plot(density(a1[,x], na.rm=T),
              xlab = paste0('Min=',min_value,', Max=',max_value),
              main=x)
       })

png('histograms.png',width = 10,height = 4,res=400, units = 'in')
par(mfrow=c(2,5))
lapply(c('cops2','bmi',"systolic" ,"diastolic","glu_f","hdl","hgba1c",
         "ldl_clc_ns","tot_choles","trigl_ns"),
       function(x){
         min_value <- min(a1[,x], na.rm=T)
         max_value <- max(a1[,x], na.rm=T)
         hist(a1[,x],breaks = 20,
              xlab = paste0('Min=',min_value,', Max=',max_value),
              main=x)
       })
dev.off()







################################################################################
# Aim 1a - K-M analysis


# for all cases v controls

# ischemic heart disease
surv_object <- Surv(time = floor(a1$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1)
png('isch_all.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease')
dev.off()

# stroke/tia
surv_object <- Surv(time = floor(a1$stroke_tia_grp_inc_fu/30), 
                    event = a1$stroke_tia_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1)
png('stroke_all.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Stroke/TIA')
dev.off()

# cardiomyopathy/heart failure
surv_object <- Surv(time = floor(a1$cardiomyopathy_heart_failure_grp_inc_fu/30), 
                    event = a1$cardiomyopathy_heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1)
png('chf_all.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Cardiomiopathy/heart failure')
dev.off()

# 3 cvd combined
surv_object <- Surv(time = floor(a1$cvdcombo_grp_inc_fu/30), 
                    event = a1$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1)
png('combo_all.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Any of the 3 outcomes')
dev.off()







################################################################################
# Aim 1b - K-M analysis

###
# for chemo cases v controls

# ischemic heart disease
surv_object <- Surv(time = floor(a1chemo$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1chemo$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1chemo)
png('isch_chemo.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1chemo,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()
# stroke/tia
surv_object <- Surv(time = floor(a1chemo$stroke_tia_grp_inc_fu/30), 
                    event = a1chemo$stroke_tia_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1chemo)
png('stroke_chemo.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1chemo,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()
# cardiomyopathy/heart failure
surv_object <- Surv(time = floor(a1chemo$cardiomyopathy_heart_failure_grp_inc_fu/30), 
                    event = a1chemo$cardiomyopathy_heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1chemo)
png('chf_chemo.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1chemo,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()
# 3 cvd combined
surv_object <- Surv(time = floor(a1chemo$cvdcombo_grp_inc_fu/30), 
                    event = a1chemo$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1chemo)
png('combo_chemo.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1chemo,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()


###
# for horm cases v controls

# ischemic heart disease
surv_object <- Surv(time = floor(a1horm$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1horm$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1horm)
png('isch_horm.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1horm,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()
# stroke/tia
surv_object <- Surv(time = floor(a1horm$stroke_tia_grp_inc_fu/30), 
                    event = a1horm$stroke_tia_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1horm)
png('stroke_horm.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1horm,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()
# cardiomyopathy/heart failure
surv_object <- Surv(time = floor(a1horm$cardiomyopathy_heart_failure_grp_inc_fu/30), 
                    event = a1horm$cardiomyopathy_heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1horm)
png('chf_horm.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1horm,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()
# 3 cvd combined
surv_object <- Surv(time = floor(a1horm$cvdcombo_grp_inc_fu/30), 
                    event = a1horm$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1horm)
png('combo_horm.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1horm,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()


###
# for radiation cases v controls

# ischemic heart disease
surv_object <- Surv(time = floor(a1rad$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1rad$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1rad)
png('isch_rad.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1rad,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()
# stroke/tia
surv_object <- Surv(time = floor(a1rad$stroke_tia_grp_inc_fu/30), 
                    event = a1rad$stroke_tia_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1rad)
png('stroke_rad.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1rad,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()
# cardiomyopathy/heart failure
surv_object <- Surv(time = floor(a1rad$cardiomyopathy_heart_failure_grp_inc_fu/30), 
                    event = a1rad$cardiomyopathy_heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1rad)
png('chf_rad.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1rad,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()
# 3 cvd combined
surv_object <- Surv(time = floor(a1rad$cvdcombo_grp_inc_fu/30), 
                    event = a1rad$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1rad)
png('combo_rad.png',width =7,height = 5.5,units = 'in',res=250)
ggsurvplot(fit1, data = a1rad,fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco')
dev.off()




################################################################################
# Aim 1a - Cox model

coxtab <- function(d,x){
  surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
  fit.coxph <- coxph(surv_object ~ group1+cops2+bmicat1+menop+parity+
                       diab+dyslipidemia2+systolic_2+diastolic_2+
                       glu_f_q+hdl_q+hgba1c_q+ldl_clc_ns_q+tot_choles_q+trigl_ns_q+
                       medhousincome_q+houspoverty_q+edu1_q+edu2_q+edu3_q+edu4_q+
                       #medhousincome+houspoverty+edu1+edu2+edu3+edu4+
                       arrhythmia_grp_prev+cardiac_arrest_grp_prev+
                        cardiomyopathy_heart_failure_grp_prev+
                        ischemic_heart_disease_grp_prev+
                        myocarditis_pericarditis_grp_prev+stroke_tia_grp_prev+
                        valvular_disease_grp_prev+venous_thromboembolic_disease_grp_prev, data = d)
  sum <- cbind(summary(fit.coxph)$coef[,c(2,5)],summary(fit.coxph)$conf.int[,3:4])
  sum[,c(1,3,4,2)]
}


# run models
cox_isch <- coxtab(a1,'ischemic_heart_disease_grp_inc')
cox_stroke <- coxtab(a1,'stroke_tia_grp_inc')
cox_chf <- coxtab(a1,'cardiomyopathy_heart_failure_grp_inc')
cox_cvdcombo <- coxtab(a1,'cvdcombo_grp_inc')


# export tables
write_csv(data.frame(cbind(row.names(cox_isch),cox_isch)),'cox_isch.csv')
write_csv(data.frame(cbind(row.names(cox_isch),cox_stroke)),'cox_stroke.csv')
write_csv(data.frame(cbind(row.names(cox_isch),cox_chf)),'cox_chf.csv')
write_csv(data.frame(cbind(row.names(cox_isch),cox_cvdcombo)),'cox_cvdcombo.csv')





################################################################################
# Aim 1b - Cox model

# chemo
cox_isch_chemo <- data.frame(coxtab(a1chemo,'ischemic_heart_disease_grp_inc'))
cox_stroke_chemo <- data.frame(coxtab(a1chemo,'stroke_tia_grp_inc'))
cox_chf_chemo <- data.frame(coxtab(a1chemo,'cardiomyopathy_heart_failure_grp_inc'))
cox_cvdcombo_chemo <- data.frame(coxtab(a1chemo,'cvdcombo_grp_inc'))

# horm
cox_isch_horm <- data.frame(coxtab(a1horm,'ischemic_heart_disease_grp_inc'))
cox_stroke_horm <- data.frame(coxtab(a1horm,'stroke_tia_grp_inc'))
cox_chf_horm <- data.frame(coxtab(a1horm,'cardiomyopathy_heart_failure_grp_inc'))
cox_cvdcombo_horm <- data.frame(coxtab(a1horm,'cvdcombo_grp_inc'))

# rad
cox_isch_rad <- data.frame(coxtab(a1rad,'ischemic_heart_disease_grp_inc'))
cox_stroke_rad <- data.frame(coxtab(a1rad,'stroke_tia_grp_inc'))
cox_chf_rad <- data.frame(coxtab(a1rad,'cardiomyopathy_heart_failure_grp_inc'))
cox_cvdcombo_rad <- data.frame(coxtab(a1rad,'cvdcombo_grp_inc'))


# export tables
write_csv(cbind(row.names(cox_isch),cox_isch_chemo),'cox_isch_chemo.csv')
write_csv(cbind(row.names(cox_isch),cox_stroke_chemo),'cox_stroke_chemo.csv')
write_csv(cbind(row.names(cox_isch),cox_chf_chemo),'cox_chf_chemo.csv')
write_csv(cbind(row.names(cox_isch),cox_cvdcombo_chemo),'cox_cvdcombo_chemo.csv')
write_csv(cox_isch_horm,'cox_isch_horm.csv')
write_csv(cox_stroke_horm,'cox_stroke_horm.csv')
write_csv(cox_chf_horm,'cox_chf_horm.csv')
write_csv(cox_cvdcombo_horm,'cox_cvdcombo_horm.csv')
write_csv(cox_isch_rad,'cox_isch_rad.csv')
write_csv(cox_stroke_rad,'cox_stroke_rad.csv')
write_csv(cox_chf_rad,'cox_chf_rad.csv')
write_csv(cox_cvdcombo_rad,'cox_cvdcombo_rad.csv')







################################################################################
# compare cases enrolled and not enrolled in Pathways


tab2 <- function(x){
  d <- a1[which(a1$group=='Case'),]
  if(class(d[,x]) %in% c('numeric')){
    n <- c(paste0(length(which(!is.na(d[,x]))),' (',
                  round(100*length(which(is.na(d[,x])))/nrow(d)),'%)'),
           tapply(d[,x], d$enrolled, function(y) 
             paste0(length(which(!is.na(y))),' (',
                    round(100*length(which(is.na(y)))/length(y)),'%)')))
    m <- c(mean(d[,x], na.rm=T),
           tapply(d[,x], d$enrolled,mean, na.rm=T))
    sd <- c(sd(d[,x], na.rm=T),
            tapply(d[,x], d$enrolled,sd, na.rm=T))
    p <- t.test(d[,x]~d$enrolled)$p.value
    sum <- c(x,n[1],m[1],sd[1],'',n[2],m[2],sd[2],'',n[3],m[3],sd[3],'',p)
  } 
  if(class(d[,x]) %in% c('factor')){
    n <- cbind(table(d[,x]),table(d[,x], d$enrolled))
    pct <- cbind(prop.table(table(d[,x])),
                 prop.table(table(d[,x], d$enrolled),2))
    p <- chisq.test(n[which(n[,1]>0),-1])$p.value
    sum <- cbind(x,'',n[,1],pct[,1],'','',n[,2],pct[,2],'','',n[,3],pct[,3],'',p)
  } 
  sum
}


table2 <- do.call(rbind,lapply(c('dxage','agegrp','raceethn1','index_yr','enr_len',
                                 "ajcc_stage", "nodal" ,"er", "pr", "erpr","her2","tri_neg",
                                 "chemo_yn","rad_yn", "horm_yn",'cops2',
                                 'bmi','bmicat','smok1','menop','gravid','parity',
                                 "inpatient_flg","outpatient_flg","pharmacy_flg","a1c_flg",
                                 "fgrg_flg","a1cfgrg_flg",'dyslipidemia.x',
                                 "systolic" ,"diastolic","glu_f","hdl","hgba1c",
                                 "ldl_clc_ns","tot_choles","trigl_ns",
                                 'hhincome1','hhincome2',"hhincome3",'medhousincome',
                                 'houspoverty', "edu1","edu2","edu3","edu4"),tab2))
tabel2 <- cbind(row.names(table2), table2)






write.csv(table1, 'table 1.csv')
write.csv(table2, 'table 2.csv')




