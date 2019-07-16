# Text extraction radiation report
# Zaixing Shi, 2019-07-02

library(stringr)

# load radiationd sample report
t1 <- read.csv('Q:/HGREENLEE/Data Working/Pathways Heart Study/Radiation data/Rad_tx_notes_type1.csv',
               stringsAsFactors=F)
t2 <- read.csv('Q:/HGREENLEE/Data Working/Pathways Heart Study/Radiation data/Rad_tx_notes_type2.csv',
               stringsAsFactors=F)


# Preprocess type 1 report
## remove line breaks
t1[,1] <- str_replace_all(t1[,1], "[\r\n]" , "")

## separate line by "?"
tlist <- strsplit(t1[,1],"\\?")

## single out radiation dose to date 
rad_dose <- sapply(tlist, function(x){
  grep('^.*?Dose|DOSE|dose.*?Gy.*?$',unlist(x),value = T)
})


# extract data from type 1 report

## cumulative dose
sapply(rad_dose, function(x){
  x <- unlist(x)
  cum_dose <- str_extract(string = x, pattern = "date[:]\\s+[0-9]+\\s+cGy|[0-9]+\\s+cGy\\s+out")
  cum_dose <- str_extract(string=cum_dose, pattern="[0-9]+")
  #total_dose[!is.na(total_dose)]
})
       
## total planned dose
unlist(sapply(rad_dose, function(x){
  x <- unlist(x)
  total_dose <- str_extract(string = x, pattern = "of\\s*[A-Za-z]*\\s*[0-9]*\\s*cGy|of\\s*[0-9]*\\s*cGy")
  total_dose <- str_extract(string=total_dose, pattern="[0-9]+")
  total_dose <- total_dose[!is.na(total_dose)]
  if(length(total_dose)==0){
    NA
  }
  else {
    total_dose
  }
}))
