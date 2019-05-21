# CVD and BC exploratory data analysis
# Zaixing Shi, 5/7/2019


library(tidyverse)
library(haven)
library(ggplot2)
library(dplyr)
library(plyr)
library(naniar)
library(reshape2)


################################################################################
# compare characteristics between cases and controls

tab1 <- function(x){
  if(class(a1[,x]) %in% c('numeric')){
    n <- c(paste0(length(which(!is.na(a1[,x]))),' (',
                  round(100*length(which(is.na(a1[,x])))/nrow(a1)),'%)'),
           tapply(a1[,x], a1$group, function(y) 
             paste0(length(which(!is.na(y))),' (',
                    round(100*length(which(is.na(y)))/length(y)),'%)')))
    m <- c(mean(a1[,x], na.rm=T),
           tapply(a1[,x], a1$group,mean, na.rm=T))
    sd <- c(sd(a1[,x], na.rm=T),
            tapply(a1[,x], a1$group,sd, na.rm=T))
    p <- t.test(a1[,x]~a1$group)$p.value
    sum <- c(x,n[1],m[1],sd[1],'',n[2],m[2],sd[2],'',n[3],m[3],sd[3],'',p)
  } 
  if(class(a1[,x]) %in% c('factor')){
    n <- cbind(table(a1[,x]),table(a1[,x], a1$group))
    pct <- cbind(prop.table(table(a1[,x])),
                 prop.table(table(a1[,x], a1$group),2))
    p <- chisq.test(table(a1[,x], a1$group))$p.value
    sum <- cbind(x,'',n[,1],pct[,1],'','',n[,2],pct[,2],'','',n[,3],pct[,3],'',p)
  } 
  sum
}


table1 <- do.call(rbind,lapply(c('dxage','agegrp','raceethn1','index_yr','enr_len','cops2',
                                 'bmi','bmicat','smok1','menop','gravid','parity',
                                 "inpatient_flg","outpatient_flg","pharmacy_flg","a1c_flg",
                                 "fgrg_flg","a1cfgrg_flg",'dyslipidemia.x',
                                 "systolic" ,"diastolic","glu_f","hdl","hgba1c",
                                 "ldl_clc_ns","tot_choles","trigl_ns",
                                 'hhincome1','hhincome2',"hhincome3",'medhousincome',
                                 'houspoverty', "edu1","edu2","edu3","edu4"
                                 ),tab1))

tabel1 <- cbind(row.names(table1), table1)









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





################################################################################
# Compare cvd events between cases and controls

cvdvars <-names(a1[,c(124:133,158:181)])

## create binary vars for prevalent cvd
a1[,paste0(cvdvars,'_p')] <- lapply(a1[,cvdvars], function(x){
  y <- factor(ifelse(x=='Prevalent','Yes','No'), levels=c('Yes','No'))
  y
})


## create binary vars for incident cvd
a1[,paste0(cvdvars,'_i')] <- lapply(a1[,cvdvars], function(x){
  y <- factor(ifelse(x=='Incident','Yes','No'), levels=c('Yes','No'))
  y
})


## compare between all cases and controls

table1_cvd <- do.call(rbind,lapply(c(paste0(cvdvars,'_p'),
                                 paste0(cvdvars,'_i')),tab1))

table1_cvd <- data.frame(cbind(row.names(table1_cvd), table1_cvd))


## compare between Pathways cases and non-PW cases

table1_cvd_pw <- do.call(rbind,lapply(c(paste0(cvdvars,'_p'),
                                     paste0(cvdvars,'_i')),tab2))

table1_cvd_pw <- data.frame(cbind(row.names(table1_cvd_pw), table1_cvd_pw))




write.csv(table1_cvd[table1_cvd$V1=='Yes',], 'table 1 cvd.csv')
write.csv(table1_cvd_pw[table1_cvd_pw$V1=='Yes',], 'table 1 cvd PW.csv')




## compare prevalent CVD by age group
cvdp_age <- data.frame(t(sapply(c(paste0(cvdvars,'_p')), function(x) prop.table(table(a1[,x], a1$agegrp),2)[1,])))
write.csv(cvdp_age, 'cvd prev by age.csv')






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















