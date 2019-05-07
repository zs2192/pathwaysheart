# Check on chemo infusion data
# 12/20/2018

library(haven)
library(tidyverse)
library(gridExtra)

sample <- read_sas('data/heart_chemodata_sample_raw.sas7bdat')


# some data cleaning
sample$cycle <- gsub('[^0-9\\.]','',sample$TX_Cycle)
sample$newstudyid <- paste0('PW',sample$studyid)
sample$dose <- as.numeric(sample$drug_recd_amt_admd)


# plot drug infusion history
p1 <- ggplot(sample[which(grepl('DOX',sample$medication_name_dspd)),])+
  geom_point(aes(x=tx_date, y=newstudyid, size=dose))+
  ggtitle('Doxorubicin')+
  theme_bw()

p2 <- ggplot(sample[which(grepl('CARBOPLATIN',sample$medication_name_dspd)),])+
  geom_point(aes(x=tx_date, y=newstudyid, size=dose))+
  ggtitle('Carboplatin')+
  theme_bw()

p3 <- ggplot(sample[which(grepl('CYCLOPHOSPHAMIDE',sample$medication_name_dspd)),])+
  geom_point(aes(x=tx_date, y=newstudyid, size=dose))+
  ggtitle('Cyclophosphamide')+
  theme_bw()

p4 <- ggplot(sample[which(grepl('TAXEL|TAXOL',sample$medication_name_dspd)),])+
  geom_point(aes(x=tx_date, y=newstudyid, size=dose))+
  ggtitle('Taxane')+
  theme_bw()

p5 <- ggplot(sample[which(grepl('STUDY|INVESTIGATIONAL',sample$medication_name_dspd)),])+
  geom_point(aes(x=tx_date, y=newstudyid, size=dose))+
  ggtitle('Study drugs')+
  theme_bw()

png("chemo_plot.png", width=15, height=5, unit='in',res=200)

grid.arrange(p1, p2, p3, p4, p5, nrow = 2)

dev.off()






# table of cycles
lapply(c('DOX','CARBOPLATIN','CYCLOPHOSPHAMIDE','TAXEL|TAXOL','STUDY|INVESTIGATIONAL'),
       function(x){
         with(sample[which(grepl(x,sample$medication_name_dspd)),],plot(TX_DAY_UNIQ_ID, studyid))
       })
