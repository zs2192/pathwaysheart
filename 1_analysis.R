# CVD and BC exploratory data analysis
# Zaixing Shi, 12/11/2019


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
varlist <- c('enr_len','cops2','bmi1','agegrp','raceethn1',
             'bmicat','smok','menop','diab','dyslipidemia',
             "systolic" ,"diastolic","glu_f","hdl","hgba1c",
             "ldl_clc_ns","tot_choles","trigl_ns",
             'medhousincome','houspoverty', "edu1","edu2","edu3","edu4")

# all sample
table1_all <- do.call(rbind,lapply(varlist,function(x) tab1(a1,x)))
table1_all <- data.frame(cbind(row.names(table1_all), table1_all))

# cases receiving chemo
table1_chemo <- do.call(rbind,lapply(varlist,function(x) tab1(a1[a1_chemo,],x)))
table1_chemo <- data.frame(cbind(row.names(table1_chemo), table1_chemo))

# cases receiving hormonal therapy
table1_horm <- do.call(rbind,lapply(varlist,function(x) tab1(a1[a1_horm,],x)))
table1_horm <- data.frame(cbind(row.names(table1_horm), table1_horm))

# cases receiving radiation
table1_rad <- do.call(rbind,lapply(varlist,function(x) tab1(a1[a1_rad,],x)))
table1_rad <- data.frame(cbind(row.names(table1_rad), table1_rad))

# export tables
write_csv(table1_all,'table1_all.csv')
write_csv(table1_chemo,'table1_chemo.csv')
write_csv(table1_horm,'table1_horm.csv')
write_csv(table1_rad,'table1_rad.csv')







################################################################################
# Compare cvd prevalence between cases and controls

# get cvd prevalence variables
cvdvars <- grep('_prev$',names(a1),value = T)
cvdvars <- cvdvars[c(1:3,5,4,6:11,c(16,17,18,21,19,13:15,20,24,26,27,28,23,12,22,25,29:33))]
# function to create
tab2 <- function(d,x){
  d[,x] <- factor(d[,x],levels=c(0,1),labels = c('No','Yes'))
  n <- cbind(table(d[,x]),table(d[,x], d$group))
  pct <- cbind(prop.table(table(d[,x])),
               prop.table(table(d[,x], d$group),2))
  p <- chisq.test(table(d[,x], d$group))$p.value
  sum <- cbind(x,'',n[,1],pct[,1],'',n[,2],pct[,2],'',n[,3],pct[,3],'',p)
  sum
}

## compare prevalent CVD by age group
cvdp_age <- data.frame(t(sapply(cvdvars, function(x) 
  prop.table(table(factor(a1[,x],levels=0:1), a1$agegrp),2)[2,])))
write_csv(cvdp_age, 'cvd_prev_agegrp.csv')

## compare incidence between all cases and controls
table2_prev_cvd <- do.call(rbind,lapply(cvdvars ,function(x) tab2(a1,x)))
table2_prev_cvd <- data.frame(cbind(row.names(table2_prev_cvd), table2_prev_cvd))

## compare incidence between cases receiving chemo and controls
table2_prev_cvd_chemo <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_chemo,],x)))
table2_prev_cvd_chemo <- data.frame(cbind(row.names(table2_prev_cvd_chemo), table2_prev_cvd_chemo))

## compare incidence between cases receiving horm and controls
table2_prev_cvd_horm <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_horm,],x)))
table2_prev_cvd_horm <- data.frame(cbind(row.names(table2_prev_cvd_horm), table2_prev_cvd_horm))

## compare incidence between cases receiving rad and controls
table2_prev_cvd_rad <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_rad,],x)))
table2_prev_cvd_rad <- data.frame(cbind(row.names(table2_prev_cvd_rad), table2_prev_cvd_rad))

# export table
write_csv(table2_prev_cvd[table2_prev_cvd$V1=='Yes',], 'table2_prev_cvd.csv')
write_csv(table2_prev_cvd_chemo[table2_prev_cvd_chemo$V1=='Yes',], 'table2_prev_cvd_chemo.csv')
write_csv(table2_prev_cvd_horm[table2_prev_cvd_horm$V1=='Yes',], 'table2_prev_cvd_horm.csv')
write_csv(table2_prev_cvd_rad[table2_prev_cvd_rad$V1=='Yes',], 'table2_prev_cvd_rad.csv')









################################################################################
# Compare cvd frequencies between cases and controls

## True incidence
cvdvars <- grep('_inc$',names(a1),value = T)
cvdvars <- cvdvars[c(1:3,5,4,6:11,c(16,17,18,21,19,13:15,20,24,26,27,28,23,12,22,25,29:33))]


## compare incidence between all cases and controls
table2_cvd <- do.call(rbind,lapply(cvdvars ,function(x) tab2(a1,x)))
table2_cvd <- data.frame(cbind(row.names(table2_cvd), table2_cvd))

## compare incidence between cases receiving chemo and controls
table2_cvd_chemo <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_chemo,],x)))
table2_cvd_chemo <- data.frame(cbind(row.names(table2_cvd_chemo), table2_cvd_chemo))

## compare incidence between cases receiving horm and controls
table2_cvd_horm <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_horm,],x)))
table2_cvd_horm <- data.frame(cbind(row.names(table2_cvd_horm), table2_cvd_horm))

## compare incidence between cases receiving rad and controls
table2_cvd_rad <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_rad,],x)))
table2_cvd_rad <- data.frame(cbind(row.names(table2_cvd_rad), table2_cvd_rad))


write_csv(table2_cvd[table2_cvd$V1=='Yes',], 'table2_cvd.csv')
write_csv(table2_cvd_chemo[table2_cvd_chemo$V1=='Yes',], 'table2_cvd_chemo.csv')
write_csv(table2_cvd_horm[table2_cvd_horm$V1=='Yes',], 'table2_cvd_horm.csv')
write_csv(table2_cvd_rad[table2_cvd_rad$V1=='Yes',], 'table2_cvd_rad.csv')



## Any new onset = true incidence + recurrence
cvdvars <- grep('_rec$',names(a1),value = T)
cvdvars <- cvdvars[c(1:3,5,4,6:11,c(16,17,18,21,19,13:15,20,24,26,27,28,23,12,22,25,29:33))]

## compare incidence between all cases and controls
table2_rec_cvd <- do.call(rbind,lapply(cvdvars ,function(x) tab2(a1,x)))
table2_rec_cvd <- data.frame(cbind(row.names(table2_rec_cvd), table2_rec_cvd))

## compare incidence between cases receiving chemo and controls
table2_rec_cvd_chemo <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_chemo,],x)))
table2_rec_cvd_chemo <- data.frame(cbind(row.names(table2_rec_cvd_chemo), table2_rec_cvd_chemo))

## compare incidence between cases receiving horm and controls
table2_rec_cvd_horm <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_horm,],x)))
table2_rec_cvd_horm <- data.frame(cbind(row.names(table2_rec_cvd_horm), table2_rec_cvd_horm))

## compare incidence between cases receiving rad and controls
table2_rec_cvd_rad <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_rad,],x)))
table2_rec_cvd_rad <- data.frame(cbind(row.names(table2_rec_cvd_rad), table2_rec_cvd_rad))


write_csv(table2_rec_cvd[table2_rec_cvd$V1=='Yes',], 'table2_rec_cvd.csv')
write_csv(table2_rec_cvd_chemo[table2_rec_cvd_chemo$V1=='Yes',], 'table2_rec_cvd_chemo.csv')
write_csv(table2_rec_cvd_horm[table2_rec_cvd_horm$V1=='Yes',], 'table2_rec_cvd_horm.csv')
write_csv(table2_rec_cvd_rad[table2_rec_cvd_rad$V1=='Yes',], 'table2_rec_cvd_rad.csv')







################################################################################
# Aim 1a - K-M analysis figures


# for all cases v controls

# true incidence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd,],fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(a1[a1_stroke,]$stroke_grp_inc_fu/30), 
                    event = a1[a1_stroke,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke,],fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Stroke/TIA',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(a1[a1_hf,]$heart_failure_grp_inc_fu/30), 
                    event = a1[a1_hf,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_hf,],fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# 3 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo,]$cvdcombo_grp_inc_fu/30), 
                    event = a1[a1_combo,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_combo,],fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Any of the 3 outcomes',xlab='Time/months')

fig1 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 2, risk.table.height = 0.2)
ggsave("fig1.png", fig1,width = 13.5, height = 12)


# any new onset=true incidence+recurrence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1$ischemic_heart_disease_grp_rec_fu/30), 
                    event = a1$ischemic_heart_disease_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1)
splots[[1]] <- ggsurvplot(fit1, data = a1,fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke/tia
surv_object <- Surv(time = floor(a1$stroke_tia_grp_rec_fu/30), 
                    event = a1$stroke_tia_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1)
splots[[2]] <- ggsurvplot(fit1, data = a1,fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke/TIA',xlab='Time/months')

# cardiomyopathy/heart failure
surv_object <- Surv(time = floor(a1$cardiomyopathy_heart_failure_grp_rec_fu/30), 
                    event = a1$cardiomyopathy_heart_failure_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1)
splots[[3]] <- ggsurvplot(fit1, data = a1,fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomiopathy/heart failure',xlab='Time/months')

# 3 cvd combined
surv_object <- Surv(time = floor(a1$cvdcombo_grp_rec_fu/30), 
                    event = a1$cvdcombo_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1)
splots[[4]] <- ggsurvplot(fit1, data = a1,fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 3 outcomes',xlab='Time/months')

fig2 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 2, risk.table.height = 0.2)
ggsave("fig2.png", fig2,width = 13.5, height = 12)







################################################################################
# Aim 1b - K-M analysis figures

###
# for chemo cases v controls

# true incidence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd_chemo,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd_chemo,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_chemo,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke/tia
surv_object <- Surv(time = floor(a1[a1_stroke_chemo,]$stroke_tia_grp_inc_fu/30), 
                    event = a1[a1_stroke_chemo,]$stroke_tia_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_chemo,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke/TIA',xlab='Time/months')

# cardiomyopathy/heart failure
surv_object <- Surv(time = floor(a1[a1_chf_chemo,]$cardiomyopathy_heart_failure_grp_inc_fu/30), 
                    event = a1[a1_chf_chemo,]$cardiomyopathy_heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_chf_chemo,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_chf_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomiopathy/heart failure',xlab='Time/months')

# 3 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo_chemo,]$cvdcombo_grp_inc_fu/30), 
                    event = a1[a1_combo_chemo,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo_chemo,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_combo_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 3 outcomes',xlab='Time/months')

fig3 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 2, risk.table.height = 0.2)
ggsave("fig3.png", fig3,width = 13.5, height = 12)


# any new onset=true incidence+recurrence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_chemo,]$ischemic_heart_disease_grp_rec_fu/30), 
                    event = a1[a1_chemo,]$ischemic_heart_disease_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_chemo,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke/tia
surv_object <- Surv(time = floor(a1[a1_chemo,]$stroke_tia_grp_rec_fu/30), 
                    event = a1[a1_chemo,]$stroke_tia_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_chemo,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke/TIA',xlab='Time/months')

# cardiomyopathy/heart failure
surv_object <- Surv(time = floor(a1[a1_chemo,]$cardiomyopathy_heart_failure_grp_rec_fu/30), 
                    event = a1[a1_chemo,]$cardiomyopathy_heart_failure_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_chemo,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomiopathy/heart failure',xlab='Time/months')

# 3 cvd combined
surv_object <- Surv(time = floor(a1[a1_chemo,]$cvdcombo_grp_rec_fu/30), 
                    event = a1[a1_chemo,]$cvdcombo_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_chemo,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 3 outcomes',xlab='Time/months')

fig4 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 2, risk.table.height = 0.2)
ggsave("fig4.png", fig4,width = 13.5, height = 12)





###
# for hormonal therapy cases v controls

# true incidence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd_horm,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd_horm,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_horm,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke/tia
surv_object <- Surv(time = floor(a1[a1_stroke_horm,]$stroke_tia_grp_inc_fu/30), 
                    event = a1[a1_stroke_horm,]$stroke_tia_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_horm,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke/TIA',xlab='Time/months')

# cardiomyopathy/heart failure
surv_object <- Surv(time = floor(a1[a1_chf_horm,]$cardiomyopathy_heart_failure_grp_inc_fu/30), 
                    event = a1[a1_chf_horm,]$cardiomyopathy_heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_chf_horm,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_chf_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomiopathy/heart failure',xlab='Time/months')

# 3 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo_horm,]$cvdcombo_grp_inc_fu/30), 
                    event = a1[a1_combo_horm,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo_horm,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_combo_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 3 outcomes',xlab='Time/months')

fig5 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 2, risk.table.height = 0.2)
ggsave("fig5.png", fig5,width = 13.5, height = 12)


# any new onset=true incidence+recurrence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_horm,]$ischemic_heart_disease_grp_rec_fu/30), 
                    event = a1[a1_horm,]$ischemic_heart_disease_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_horm,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke/tia
surv_object <- Surv(time = floor(a1[a1_horm,]$stroke_tia_grp_rec_fu/30), 
                    event = a1[a1_horm,]$stroke_tia_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_horm,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke/TIA',xlab='Time/months')

# cardiomyopathy/heart failure
surv_object <- Surv(time = floor(a1[a1_horm,]$cardiomyopathy_heart_failure_grp_rec_fu/30), 
                    event = a1[a1_horm,]$cardiomyopathy_heart_failure_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_horm,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomiopathy/heart failure',xlab='Time/months')

# 3 cvd combined
surv_object <- Surv(time = floor(a1[a1_horm,]$cvdcombo_grp_rec_fu/30), 
                    event = a1[a1_horm,]$cvdcombo_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_horm,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 3 outcomes',xlab='Time/months')

fig6 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 2, risk.table.height = 0.2)
ggsave("fig6.png", fig6,width = 13.5, height = 12)





###
# for radiation cases v controls

# true incidence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd_rad,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd_rad,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_rad,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke/tia
surv_object <- Surv(time = floor(a1[a1_stroke_rad,]$stroke_tia_grp_inc_fu/30), 
                    event = a1[a1_stroke_rad,]$stroke_tia_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_rad,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke/TIA',xlab='Time/months')

# cardiomyopathy/heart failure
surv_object <- Surv(time = floor(a1[a1_chf_rad,]$cardiomyopathy_heart_failure_grp_inc_fu/30), 
                    event = a1[a1_chf_rad,]$cardiomyopathy_heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_chf_rad,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_chf_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomiopathy/heart failure',xlab='Time/months')

# 3 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo_rad,]$cvdcombo_grp_inc_fu/30), 
                    event = a1[a1_combo_rad,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo_rad,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_combo_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 3 outcomes',xlab='Time/months')

fig7 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 2, risk.table.height = 0.2)
ggsave("fig7.png", fig7,width = 13.5, height = 12)


# any new onset=true incidence+recurrence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_rad,]$ischemic_heart_disease_grp_rec_fu/30), 
                    event = a1[a1_rad,]$ischemic_heart_disease_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_rad,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke/tia
surv_object <- Surv(time = floor(a1[a1_rad,]$stroke_tia_grp_rec_fu/30), 
                    event = a1[a1_rad,]$stroke_tia_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_rad,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke/TIA',xlab='Time/months')

# cardiomyopathy/heart failure
surv_object <- Surv(time = floor(a1[a1_rad,]$cardiomyopathy_heart_failure_grp_rec_fu/30), 
                    event = a1[a1_rad,]$cardiomyopathy_heart_failure_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_rad,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomiopathy/heart failure',xlab='Time/months')

# 3 cvd combined
surv_object <- Surv(time = floor(a1[a1_rad,]$cvdcombo_grp_rec_fu/30), 
                    event = a1[a1_rad,]$cvdcombo_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_rad,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 3 outcomes',xlab='Time/months')

fig8 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 2, risk.table.height = 0.2)
ggsave("fig8.png", fig8,width = 13.5, height = 12)









################################################################################
# Aim 1a - Cox model

# function to summary cox model results
coxtab <- function(d,x,covar){
  covars <- paste(covar, collapse = '+')
  surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
  fit.coxph1 <- coxph(surv_object ~ group1, data = d)
  fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ group1+bmicat1+menop+parity+
                        diab+systolic_2+diastolic_2+glu_f_q+medhousincome_q+edu1_q+",
                                        covars,", data = d)"))))
  sum1 <- c(summary(fit.coxph1)$coef[c(2,5)],summary(fit.coxph1)$conf.int[3:4])
  sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
  sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
  sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
  sum4$var <- factor(sum4$var, levels=row.names(sum3))
  sum4
}

# list of prevalent cvd
prevcvd <-c("arrhythmia_grp_prev","cardiomyopathy_heart_failure_grp_prev",
            "ischemic_heart_disease_grp_prev","myocarditis_pericarditis_grp_prev",    
            "stroke_tia_grp_prev","valvular_disease_grp_prev","venous_thromboembolic_disease_grp_prev")

###
# run models in CASES and CONTROL 1
# Pure incidence
cox_isch <- coxtab(a1[a1_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[-3])
cox_stroke <- coxtab(a1[a1_stroke,],'stroke_tia_grp_inc',prevcvd[-5])
cox_chf <- coxtab(a1[a1_chf,],'cardiomyopathy_heart_failure_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-5)])
cox1 <- Reduce(function(...) merge(...,by='var', all=T),
               list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox1[,-1] <- lapply(cox1[,-1], as.character)
cox1 <- rbind(cox1[1,],'',cox1[2,],'',cox1[3:7,],cox1[8,],'',cox1[9:11,],cox1[12,],'',
              cox1[13:14,],'',cox1[15:16,],'',cox1[17:20,],'',cox1[21:24,],'',
              cox1[25:27,],'',cox1[28:34,])

cox1 <- cbind(cox1[,1:5],'',cox1[,6:9],'',cox1[,10:13],'',cox1[,14:17])

# Any new onset
cox_isch <- coxtab(a1,'ischemic_heart_disease_grp_rec',prevcvd)
cox_stroke <- coxtab(a1,'stroke_tia_grp_rec',prevcvd)
cox_chf <- coxtab(a1,'cardiomyopathy_heart_failure_grp_rec',prevcvd)
cox_cvdcombo <- coxtab(a1,'cvdcombo_grp_rec',prevcvd)
cox1rec <- Reduce(function(...) merge(...,by='var', all=T),
               list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox1rec$var <- factor(cox1rec$var, levels=levels(cox_isch$var))
cox1rec[,-1] <- lapply(cox1rec[,-1], as.character)
cox1rec <- cox1rec[order(cox1rec$var),]
cox1rec <- rbind(cox1rec[1,],'',cox1rec[2,],'',cox1rec[3:7,],cox1rec[8,],'',cox1rec[9:11,],cox1rec[12,],'',
              cox1rec[13:14,],'',cox1rec[15:16,],'',cox1rec[17:20,],'',cox1rec[21:24,],'',
              cox1rec[25:27,],'',cox1rec[28:34,])
cox1rec <- cbind(cox1rec[,1:5],'',cox1rec[,6:9],'',cox1rec[,10:13],'',cox1rec[,14:17])


###
# run models in CASES and CONTROL 2

# Pure incidence
a2_cox_isch <- coxtab(a2[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[-3])
a2_cox_stroke <- coxtab(a2[a2_stroke,],'stroke_tia_grp_inc',prevcvd[-5])
a2_cox_chf <- coxtab(a2[a2_chf,],'cardiomyopathy_heart_failure_grp_inc',prevcvd[-2])
a2_cox_cvdcombo <- coxtab(a2[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-5)])
a2_cox1 <- Reduce(function(...) merge(...,by='var', all=T),
               list(a2_cox_isch, a2_cox_stroke, a2_cox_chf, a2_cox_cvdcombo))
a2_cox1[,-1] <- lapply(a2_cox1[,-1], as.character)
a2_cox1 <- rbind(a2_cox1[1,],'',a2_cox1[2,],'',a2_cox1[3:7,],a2_cox1[8,],'',a2_cox1[9:11,],a2_cox1[12,],'',
              a2_cox1[13:14,],'',a2_cox1[15:16,],'',a2_cox1[17:20,],'',a2_cox1[21:24,],'',
              a2_cox1[25:27,],'',a2_cox1[28:34,])
a2_cox1 <- cbind(a2_cox1[,1:5],'',a2_cox1[,6:9],'',a2_cox1[,10:13],'',a2_cox1[,14:17])

# Any new onset
a2_cox_isch <- coxtab(a2,'ischemic_heart_disease_grp_rec',prevcvd)
a2_cox_stroke <- coxtab(a2,'stroke_tia_grp_rec',prevcvd)
a2_cox_chf <- coxtab(a2,'cardiomyopathy_heart_failure_grp_rec',prevcvd)
a2_cox_cvdcombo <- coxtab(a2,'cvdcombo_grp_rec',prevcvd)
a2_cox1rec <- Reduce(function(...) merge(...,by='var', all=T),
                  list(a2_cox_isch, a2_cox_stroke, a2_cox_chf, a2_cox_cvdcombo))
a2_cox1rec[,-1] <- lapply(a2_cox1rec[,-1], as.character)
a2_cox1rec <- a2_cox1rec[order(a2_cox1rec$var),]
a2_cox1rec <- rbind(a2_cox1rec[1,],'',a2_cox1rec[2,],'',a2_cox1rec[3:7,],a2_cox1rec[8,],'',a2_cox1rec[9:11,],a2_cox1rec[12,],'',
                 a2_cox1rec[13:14,],'',a2_cox1rec[15:16,],'',a2_cox1rec[17:20,],'',a2_cox1rec[21:24,],'',
                 a2_cox1rec[25:27,],'',a2_cox1rec[28:34,])
a2_cox1rec <- cbind(a2_cox1rec[,1:5],'',a2_cox1rec[,6:9],'',a2_cox1rec[,10:13],'',a2_cox1rec[,14:17])



# export tables as csv file

write_csv(cox1,'cox_all_incid.csv')
write_csv(cox1rec,'cox_all_newonset.csv')
write_csv(a2_cox1,'cox_cntrl2_incid.csv')
write_csv(a2_cox1rec,'cox_cntrl2_newonset.csv')







################################################################################
# Aim 1b - Cox model

###
# chemo cases
# Pure incidence
cox_isch <- coxtab(a1[a1_ihd_chemo,],'ischemic_heart_disease_grp_inc',prevcvd[-3])
cox_stroke <- coxtab(a1[a1_stroke_chemo,],'stroke_tia_grp_inc',prevcvd[-5])
cox_chf <- coxtab(a1[a1_chf_chemo,],'cardiomyopathy_heart_failure_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_chemo,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-5)])
cox2 <- Reduce(function(...) merge(...,by='var', all=T),
               list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox2[,-1] <- lapply(cox2[,-1], as.character)
cox2 <- rbind(cox2[1,],'',cox2[2,],'',cox2[3:7,],cox2[8,],'',cox2[9:11,],cox2[12,],'',
              cox2[13:14,],'',cox2[15:16,],'',cox2[17:20,],'',cox2[21:24,],'',
              cox2[25:27,],'',cox2[28:34,])
cox2 <- cbind(cox2[,1:5],'',cox2[,6:9],'',cox2[,10:13],'',cox2[,14:17])

# Any new onset
cox_isch <- coxtab(a1[a1_chemo,],'ischemic_heart_disease_grp_rec',prevcvd)
cox_stroke <- coxtab(a1[a1_chemo,],'stroke_tia_grp_rec',prevcvd)
cox_chf <- coxtab(a1[a1_chemo,],'cardiomyopathy_heart_failure_grp_rec',prevcvd)
cox_cvdcombo <- coxtab(a1[a1_chemo,],'cvdcombo_grp_rec',prevcvd)
cox2rec <- Reduce(function(...) merge(...,by='var', all=T),
                  list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox2rec$var <- factor(cox2rec$var, levels=levels(cox_isch$var))
cox2rec[,-1] <- lapply(cox2rec[,-1], as.character)
cox2rec <- cox2rec[order(cox2rec$var),]
cox2rec <- rbind(cox2rec[1,],'',cox2rec[2,],'',cox2rec[3:7,],cox2rec[8,],'',cox2rec[9:11,],cox2rec[12,],'',
                 cox2rec[13:14,],'',cox2rec[15:16,],'',cox2rec[17:20,],'',cox2rec[21:24,],'',
                 cox2rec[25:27,],'',cox2rec[28:34,])
cox2rec <- cbind(cox2rec[,1:5],'',cox2rec[,6:9],'',cox2rec[,10:13],'',cox2rec[,14:17])



# horm
# Pure incidence
cox_isch <- coxtab(a1[a1_ihd_horm,],'ischemic_heart_disease_grp_inc',prevcvd[-3])
cox_stroke <- coxtab(a1[a1_stroke_horm,],'stroke_tia_grp_inc',prevcvd[-5])
cox_chf <- coxtab(a1[a1_chf_horm,],'cardiomyopathy_heart_failure_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_horm,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-5)])
cox3 <- Reduce(function(...) merge(...,by='var', all=T),
               list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox3[,-1] <- lapply(cox3[,-1], as.character)
cox3 <- rbind(cox3[1,],'',cox3[2,],'',cox3[3:7,],cox3[8,],'',cox3[9:11,],cox3[12,],'',
              cox3[13:14,],'',cox3[15:16,],'',cox3[17:20,],'',cox3[21:24,],'',
              cox3[25:27,],'',cox3[28:34,])
cox3 <- cbind(cox3[,1:5],'',cox3[,6:9],'',cox3[,10:13],'',cox3[,14:17])

# Any new onset
cox_isch <- coxtab(a1[a1_horm,],'ischemic_heart_disease_grp_rec',prevcvd)
cox_stroke <- coxtab(a1[a1_horm,],'stroke_tia_grp_rec',prevcvd)
cox_chf <- coxtab(a1[a1_horm,],'cardiomyopathy_heart_failure_grp_rec',prevcvd)
cox_cvdcombo <- coxtab(a1[a1_horm,],'cvdcombo_grp_rec',prevcvd)
cox3rec <- Reduce(function(...) merge(...,by='var', all=T),
                  list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox3rec$var <- factor(cox3rec$var, levels=levels(cox_isch$var))
cox3rec[,-1] <- lapply(cox3rec[,-1], as.character)
cox3rec <- cox3rec[order(cox3rec$var),]
cox3rec <- rbind(cox3rec[1,],'',cox3rec[2,],'',cox3rec[3:7,],cox3rec[8,],'',cox3rec[9:11,],cox3rec[12,],'',
                 cox3rec[13:14,],'',cox3rec[15:16,],'',cox3rec[17:20,],'',cox3rec[21:24,],'',
                 cox3rec[25:27,],'',cox3rec[28:34,])
cox3rec <- cbind(cox3rec[,1:5],'',cox3rec[,6:9],'',cox3rec[,10:13],'',cox3rec[,14:17])


# rad
# Pure incidence
cox_isch <- coxtab(a1[a1_ihd_rad,],'ischemic_heart_disease_grp_inc',prevcvd[-3])
cox_stroke <- coxtab(a1[a1_stroke_rad,],'stroke_tia_grp_inc',prevcvd[-5])
cox_chf <- coxtab(a1[a1_chf_rad,],'cardiomyopathy_heart_failure_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_rad,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-5)])
cox4 <- Reduce(function(...) merge(...,by='var', all=T),
               list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox4[,-1] <- lapply(cox4[,-1], as.character)
cox4 <- rbind(cox4[1,],'',cox4[2,],'',cox4[3:7,],cox4[8,],'',cox4[9:11,],cox4[12,],'',
              cox4[13:14,],'',cox4[15:16,],'',cox4[17:20,],'',cox4[21:24,],'',
              cox4[25:27,],'',cox4[28:34,])
cox4 <- cbind(cox4[,1:5],'',cox4[,6:9],'',cox4[,10:13],'',cox4[,14:17])

# Any new onset
cox_isch <- coxtab(a1[a1_rad,],'ischemic_heart_disease_grp_rec',prevcvd)
cox_stroke <- coxtab(a1[a1_rad,],'stroke_tia_grp_rec',prevcvd)
cox_chf <- coxtab(a1[a1_rad,],'cardiomyopathy_heart_failure_grp_rec',prevcvd)
cox_cvdcombo <- coxtab(a1[a1_rad,],'cvdcombo_grp_rec',prevcvd)
cox4rec <- Reduce(function(...) merge(...,by='var', all=T),
                  list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox4rec$var <- factor(cox4rec$var, levels=levels(cox_isch$var))
cox4rec[,-1] <- lapply(cox4rec[,-1], as.character)
cox4rec <- cox4rec[order(cox4rec$var),]
cox4rec <- rbind(cox4rec[1,],'',cox4rec[2,],'',cox4rec[3:7,],cox4rec[8,],'',cox4rec[9:11,],cox4rec[12,],'',
                 cox4rec[13:14,],'',cox4rec[15:16,],'',cox4rec[17:20,],'',cox4rec[21:24,],'',
                 cox4rec[25:27,],'',cox4rec[28:34,])
cox4rec <- cbind(cox4rec[,1:5],'',cox4rec[,6:9],'',cox4rec[,10:13],'',cox4rec[,14:17])


# left side tumor, rad
# Pure incidence
cox_isch <- coxtab(a1[a1_ihd_rad_l,],'ischemic_heart_disease_grp_inc',prevcvd[-3])
cox_stroke <- coxtab(a1[a1_stroke_rad_l,],'stroke_tia_grp_inc',prevcvd[-5])
cox_chf <- coxtab(a1[a1_chf_rad_l,],'cardiomyopathy_heart_failure_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_rad_l,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-5)])
cox5 <- Reduce(function(...) merge(...,by='var', all=T),
               list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox5[,-1] <- lapply(cox5[,-1], as.character)
cox5 <- rbind(cox5[1,],'',cox5[2,],'',cox5[3:7,],cox5[8,],'',cox5[9:11,],cox5[12,],'',
              cox5[13:14,],'',cox5[15:16,],'',cox5[17:20,],'',cox5[21:24,],'',
              cox5[25:27,],'',cox5[28:34,])
cox5 <- cbind(cox5[,1:5],'',cox5[,6:9],'',cox5[,10:13],'',cox5[,14:17])

# Any new onset
cox_isch <- coxtab(a1[a1_rad_l,],'ischemic_heart_disease_grp_rec',prevcvd)
cox_stroke <- coxtab(a1[a1_rad_l,],'stroke_tia_grp_rec',prevcvd)
cox_chf <- coxtab(a1[a1_rad_l,],'cardiomyopathy_heart_failure_grp_rec',prevcvd)
cox_cvdcombo <- coxtab(a1[a1_rad_l,],'cvdcombo_grp_rec',prevcvd)
cox5rec <- Reduce(function(...) merge(...,by='var', all=T),
                  list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox5rec$var <- factor(cox5rec$var, levels=levels(cox_isch$var))
cox5rec[,-1] <- lapply(cox5rec[,-1], as.character)
cox5rec <- cox5rec[order(cox5rec$var),]
cox5rec <- rbind(cox5rec[1,],'',cox5rec[2,],'',cox5rec[3:7,],cox5rec[8,],'',cox5rec[9:11,],cox5rec[12,],'',
                 cox5rec[13:14,],'',cox5rec[15:16,],'',cox5rec[17:20,],'',cox5rec[21:24,],'',
                 cox5rec[25:27,],'',cox5rec[28:34,])
cox5rec <- cbind(cox5rec[,1:5],'',cox5rec[,6:9],'',cox5rec[,10:13],'',cox5rec[,14:17])


# right side tumor, rad
# Pure incidence
cox_isch <- coxtab(a1[a1_ihd_rad_r,],'ischemic_heart_disease_grp_inc',prevcvd[-3])
cox_stroke <- coxtab(a1[a1_stroke_rad_r,],'stroke_tia_grp_inc',prevcvd[-5])
cox_chf <- coxtab(a1[a1_chf_rad_r,],'cardiomyopathy_heart_failure_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_rad_r,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-5)])
cox6 <- Reduce(function(...) merge(...,by='var', all=T),
               list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox6[,-1] <- lapply(cox6[,-1], as.character)
cox6 <- rbind(cox6[1,],'',cox6[2,],'',cox6[3:7,],cox6[8,],'',cox6[9:11,],cox6[12,],'',
              cox6[13:14,],'',cox6[15:16,],'',cox6[17:20,],'',cox6[21:24,],'',
              cox6[25:27,],'',cox6[28:34,])
cox6 <- cbind(cox6[,1:5],'',cox6[,6:9],'',cox6[,10:13],'',cox6[,14:17])

# Any new onset
cox_isch <- coxtab(a1[a1_rad_r,],'ischemic_heart_disease_grp_rec',prevcvd)
cox_stroke <- coxtab(a1[a1_rad_r,],'stroke_tia_grp_rec',prevcvd)
cox_chf <- coxtab(a1[a1_rad_r,],'cardiomyopathy_heart_failure_grp_rec',prevcvd)
cox_cvdcombo <- coxtab(a1[a1_rad_r,],'cvdcombo_grp_rec',prevcvd)
cox6rec <- Reduce(function(...) merge(...,by='var', all=T),
                  list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
cox6rec$var <- factor(cox6rec$var, levels=levels(cox_isch$var))
cox6rec[,-1] <- lapply(cox6rec[,-1], as.character)
cox6rec <- cox6rec[order(cox6rec$var),]
cox6rec <- rbind(cox6rec[1,],'',cox6rec[2,],'',cox6rec[3:7,],cox6rec[8,],'',cox6rec[9:11,],cox6rec[12,],'',
                 cox6rec[13:14,],'',cox6rec[15:16,],'',cox6rec[17:20,],'',cox6rec[21:24,],'',
                 cox6rec[25:27,],'',cox6rec[28:34,])
cox6rec <- cbind(cox6rec[,1:5],'',cox6rec[,6:9],'',cox6rec[,10:13],'',cox6rec[,14:17])


# export tables
write_csv(cox2,'cox_chemo_incid.csv')
write_csv(cox2rec,'cox_chemo_newonset.csv')
write_csv(cox3,'cox_horm_incid.csv')
write_csv(cox3rec,'cox_horm_newonset.csv')
write_csv(cox4,'cox_rad_incid.csv')
write_csv(cox4rec,'cox_rad_newonset.csv')
write_csv(cox5,'cox_rad_left_incid.csv')
write_csv(cox5rec,'cox_rad_left_newonset.csv')
write_csv(cox6,'cox_rad_right_incid.csv')
write_csv(cox6rec,'cox_rad_right_newonset.csv')







################################################################################
# Sample size for each analysis

tab1ns <- data.frame(rbind(c(nrow(a1),unlist(table(a1$group))),
                c(nrow(a1[a1_chemo,]),unlist(table(a1[a1_chemo,]$group))),
                c(nrow(a1[a1_horm,]),unlist(table(a1[a1_horm,]$group))),
                c(nrow(a1[a1_rad,]),unlist(table(a1[a1_rad,]$group)))))
names(tab1ns) <- c('all','cases','controls')
tab1ns$sample <- c('all','chemo','hormonal','radiation')


coxns <- data.frame(rbind(c(nrow(a1),sapply(list(a1_ihd,a1_stroke,a1_chf,a1_combo),length)),
               c(nrow(a2),sapply(list(a2_ihd,a2_stroke,a2_chf,a2_combo),length)),
               c(length(a1_chemo),sapply(list(a1_ihd_chemo,a1_stroke_chemo,a1_chf_chemo,a1_combo_chemo),length)),
               c(length(a1_horm),sapply(list(a1_ihd_horm,a1_stroke_horm,a1_chf_horm,a1_combo_horm),length)),
               c(length(a1_rad),sapply(list(a1_ihd_rad,a1_stroke_rad,a1_chf_rad,a1_combo_rad),length)),
               c(length(a1_rad_l),sapply(list(a1_ihd_rad_l,a1_stroke_rad_l,a1_chf_rad_l,a1_combo_rad_l),length)),
               c(length(a1_rad_r),sapply(list(a1_ihd_rad_r,a1_stroke_rad_r,a1_chf_rad_r,a1_combo_rad_r),length))))
names(coxns) <- c('all','no IHD','no stroke tia','no CHF','no any 3 cvd')
coxns$sample <- c('all','all control2','chemo','hormonal','radiation','left tumor radiation','right tumor radiation')

# export sample sizes
write_csv(tab1ns, 'sample_size_tab1.csv')
write_csv(coxns, 'sample_size_cox_models.csv')




################################################################################
#
#   Lab data quality check - not for main analysis
#
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







