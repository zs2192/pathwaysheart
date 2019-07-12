# Text extraction radiation report
# Zaixing Shi, 2019-07-02

library(stringr)


t1 <- read.csv('C:/Users/Zaixing/OneDrive - Fred Hutchinson Cancer Research Center/0. Pathways/CVD project/Data/Rad_tx_notes_type1.csv',
               stringsAsFactors=F)

# remove line breaks
t1[,1] <- str_replace_all(t1[,1], "[\r\n]" , "")

# separate line by "?"
tlist <- strsplit(t1[,1],"\\?")

# single out radiation dose to date 
rad_dose <- sapply(tlist, function(x){
  grep('Radiation dose|^Dose:|^DOSE:',unlist(x),value = T)
})

# single out radiation dose to date 
rad_dose <- sapply(tlist, function(x){
  length(grep('Gy',unlist(x),value = T))
})


# extract