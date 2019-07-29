# Text extraction radiation report
# Zaixing Shi, 2019-07-27

library(stringr)

################################################################################
# load radiationd sample report
t1 <- read.csv('Q:/HGREENLEE/Data Working/Pathways Heart Study/Radiation data/Rad_tx_notes_type1.csv',
               stringsAsFactors=F)
t2 <- read.csv('Q:/HGREENLEE/Data Working/Pathways Heart Study/Radiation data/Rad_tx_notes_type2.csv',
               stringsAsFactors=F)




################################################################################
# extract data from type 1 report


t1_new <-  lapply(1:nrow(t1), function(x){
  ## remove line breaks
  txt <- str_replace_all(t1[x,1], "[\r\n]" , "")
  
  ## separate line by "?"
  txtlist <- strsplit(txt,"\\?")
  
  ## single out radiation dose to date 
  txtlist <- grep('^.*?Dose|DOSE|dose.*?Gy.*?$',unlist(txtlist),value = T)
  
  # extract data
  lapply(txtlist, function(y){
    t <- y
    # cum dose
    cd <- str_extract(string = y, pattern = "(dose|Dose|DOSE)\\s*[A-Za-z ]*[:]\\s*\\d+\\s*cGy|\\d+\\s*cGy\\s+out")
    cd <- str_extract(string=cd, pattern="\\d+")
    # total dose
    td <- str_extract(string = y, pattern = "of\\s*[A-Za-z ]*\\s*\\d*\\s*cGy|of\\s*\\d*\\s*cGy")
    td <- str_extract(string=td, pattern="\\d+")
    # fraction 
    f <- str_extract(string = y, pattern = "\\d*\\s*([A-Za-z ]*|[A-Za-z ]*\\/)fraction")

    if(grepl('cGy',f)) {
      # dose per frac 
      df <- str_extract(string=f, pattern="\\d+")
      # num of frac 
      nf <- NA
    } else {
      df <- NA
      nf <- str_extract(string=f, pattern="\\d+")
    }
    # location
    l <- str_extract(string = y, pattern = "([S,s]ite\\s*.*?[:]|\\sto)[A-Za-z ]*(breast|chest wall|brain|supraclavicular|axillary)")
    l <- gsub('(to the\\s*)|(to\\s*)','',l)
    # technique
    tec <- str_extract_all(string = y, pattern = "(\\d+\\s*MV[A-Za-z]*\\s*[A-Za-z ]*photon|photon|\\d+\\s*(MeV|MEV)\\s*electron|electron|d+\\s*MEV\\s*.*beam|beam|subsegmental\\s*integrated\\s*fields|,\\susing\\s*.*technique)")
    tec <- paste(unique(tec[[1]]),collapse = ', ')
    # start date
    sd <- str_extract(string = y, pattern = "(started|Started|began|Began)\\s*[A-Za-z ]*\\d*\\/\\d*\\/\\d*")
    sd <- str_extract(string = sd, pattern = "\\d*\\/\\d*\\/\\d*")
    # end date
    ed <- str_extract(string = y, pattern = "(complete|anticipate|fini)\\s*[A-Za-z ]*\\d*\\/\\d*\\/\\d*")
    ed <- str_extract(string = ed, pattern = "\\d*\\/\\d*\\/\\d*")

    # summary
    sum <- c(txt,t,cd,td,nf,df,l,tec,sd,ed)
    sum
  })
})

t1_new <- data.frame(do.call(rbind,lapply(t1_new, function(x) do.call(rbind,x))))

names(t1_new) <- c('report','treatment','cummlative_dose','total_dose',
                   'fraction_num','fraction_dose','location','technique',
                   'start_date','end_date')



write.csv(t1_new, 'Q:/HGREENLEE/Data Working/Pathways Heart Study/Radiation data/t1_extract_20190728.csv')

################################################################################
# extract data from type 2 report


t2_new <- lapply(1:nrow(t2),function(x){
  ## remove line breaks
  txt <- str_replace_all(t2[x,1], "[\r\n]" , "")
  
  ## separate line by "? ? "
  t2list <- strsplit(txt,"(\\?\\s){2}")
  
  ## select treatment summary section
  t2rad_dose <- grep('(Treatment)\\s*(Summary)',unlist(t2list),value = T)
  
  ## separate line by "?"
  t2rad_dose <- strsplit(t2rad_dose,"\\?")[[1]]
  
  # further restrict to lines containing technical details
  t2rad_dose2 <- t2rad_dose[grep("(.*[0-9]+X)|(.*[0-9]+E)",t2rad_dose)]

  # clean up spaces and commas
  t2rad_dose2 <- trimws(gsub("\\s\\/\\s",'\\/',t2rad_dose2))
  t2rad_dose2 <- trimws(gsub(",",'',t2rad_dose2))
  
  # course of RT
  c <- grep("Course[:]\\s*", t2rad_dose, value=TRUE)
  c <- str_extract(string = c, pattern = "\\d")
  
  sum <- lapply(t2rad_dose2, function(y){
    # Site and energy
    d <- str_extract_all(string = y, pattern = "(.*\\d+X)|(.*\\d+E)")
    s <- trimws(gsub("(\\d+X)|(\\d+E)|\\/|\\\\",'',d))
    e <- paste(unlist(str_extract_all(string = d, pattern = "(\\d+X)|(\\d+E)")),collapse = '/')
    
    # technique
    tec <- str_extract_all(string = y, pattern = "(\\d+X|\\d+E)\\s*.*?[0-9.]{3}")
    tec <- trimws(gsub("(\\d+X)|(\\d+E)|\\d{3}|\\/|\\\\",'',tec))
    
    # fraction
    f <- str_extract(string =y,'^(.*?)\\d{1,2}\\/\\d{1,2}')# number of fraction
    f <- str_extract(string =f,'(\\d{1,4}|\\d{1,4}\\@\\d+%)\\s*\\d{1,2}\\/\\d{1,2}')
    f <- unlist(strsplit(f,'\\s|\\/'))
    f <- f[!is.na(f)]
    if(length(f)==2){
      frac = c(gsub(paste0(f[2],'$'),'',f[1]),f[2],f[2])
    } else {
      if(length(f)==0){
        frac=NA
      }
      else {
        frac =unlist(f)
      }
    }
    
    # start and end dates
    dt <- unlist(str_extract_all(string = y, pattern = "\\d{1,2}\\/\\d{1,2}\\/(1|0)\\d{1}|\\d{1,2}\\/\\d{1,2}\\/2\\d{3}"))
    if(length(dt)==0){
      NA
    }
    else {
      if(nchar(dt[1])==8){ as.Date(dt,'%m/%d/%y')} else {
        as.Date(dt,'%m/%d/%Y')
      }
    }
    c(txt, y, c,s,e,tec,frac,dt)
  })
  
  sum
})


t2_new <- data.frame(do.call(rbind,lapply(t2_new, function(x) do.call(rbind,x))))
names(t2_new) <- c('report','treatment','course','locatin','energy','technique',
                   'fraction_dose','fraction_num','fraction_planned',
                   'start_date','end_date')
t2_new$fraction_dose <- as.numeric(as.character(t2_new$fraction_dose)) 
t2_new$fraction_num <- as.numeric(as.character(t2_new$fraction_num))
t2_new$total_dose <- t2_new$fraction_dose*t2_new$fraction_num

t2_new <- t2_new[,c('report','treatment','course','locatin','energy','technique',
                    'fraction_dose','fraction_num','fraction_planned',
                    'total_dose','start_date','end_date')]

write.csv(t2_new, 'Q:/HGREENLEE/Data Working/Pathways Heart Study/Radiation data/t2_extract_20190728.csv')



