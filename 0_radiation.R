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
cum_dose <- unlist(sapply(rad_dose, function(x){
  x <- unlist(x)
  d <- str_extract(string = x, pattern = "date[:]\\s*[0-9]+\\s*cGy|[0-9]+\\s*cGy\\s+out")
  d <- str_extract(string=d, pattern="[0-9]+")
  d <- d[!is.na(d)]
  if(length(d)==0){
    NA
  }
  else {
    d
  }
}))
       
## total planned dose
total_dose <- unlist(sapply(rad_dose, function(x){
  x <- unlist(x)
  d <- str_extract(string = x, pattern = "of\\s*[A-Za-z]*\\s*[0-9]*\\s*cGy|of\\s*[0-9]*\\s*cGy")
  d <- str_extract(string=d, pattern="[0-9]+")
  d <- d[!is.na(d)]
  if(length(d)==0){
    NA
  }
  else {
    d
  }
}))

## number of fractions or dose per fraction
fraction_dose <- unlist(sapply(rad_dose, function(x){
  x <- unlist(x)
  d <- str_extract(string = x, pattern = "[0-9]*\\s*([A-Za-z ]*|[A-Za-z ]*\\/)fraction")
  d <- str_extract(string=d[grepl('cGy',d)], pattern="[0-9]+")
  #d <- str_extract(string=d[!grepl('cGy',d)], pattern="[0-9]+")
  d <- d[!is.na(d)]
  if(length(d)==0){
    NA
  }
  else {
    d
  }
}))

fraction_num <- unlist(sapply(rad_dose, function(x){
  x <- unlist(x)
  d <- str_extract(string = x, pattern = "[0-9]*\\s*([A-Za-z ]*|[A-Za-z ]*\\/)fraction")
  d <- str_extract(string=d[!grepl('cGy',d)], pattern="[0-9]+")
  d <- unique(d[!is.na(d)])[1]
  if(length(d)==0){
    NA
  }
  else {
    d
  }
}))

## location
location <- unlist(sapply(rad_dose, function(x){
  x <- unlist(x)
  d <- str_extract(string = x, pattern = "\\sto[A-Za-z ]*(breast|chest wall|brain)")
  d <- gsub('(to the\\s*)|(to\\s*)','',d)
  d <- d[!is.na(d)]
  if(length(d)==0){
    NA
  }
  else {
    d
  }
}))

## start date
start_date <- unlist(sapply(rad_dose, function(x){
  x <- unlist(x)
  d <- str_extract(string = x, pattern = "(started|Started|began|Began)\\s*[A-Za-z ]*[0-9]*\\/[0-9]*\\/[0-9]*")
  d <- str_extract(string = d, pattern = "[0-9]*\\/[0-9]*\\/[0-9]*")
  d <- d[!is.na(d)]
  if(length(d)==0){
    NA
  }
  else {
    d
  }
}))

## end date
end_date <- unlist(sapply(rad_dose, function(x){
  x <- unlist(x)
  d <- str_extract(string = x, pattern = "(complete|anticipate|fini)\\s*[A-Za-z ]*[0-9]*\\/[0-9]*\\/[0-9]*")
  d <- str_extract(string = d, pattern = "[0-9]*\\/[0-9]*\\/[0-9]*")
  d <- d[!is.na(d)]
  if(length(d)==0){
    NA
  }
  else {
    d
  }
}))


# combine extracted data with original data
t1_new <- cbind(t1, cum_dose, total_dose, fraction_dose, fraction_num, location, start_date, end_date)

write.csv(t1_new, 'Q:/HGREENLEE/Data Working/Pathways Heart Study/Radiation data/t1_extract.csv')

