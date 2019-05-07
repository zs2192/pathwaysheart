# create treatment - cvd link

tox <- read.csv('C:/Users/zshi/Dropbox/0. Pathways/CVD project/Reviews/Cardio-Oncology slides/Table 1.csv',
                stringsAsFactors = F)

toxpair <- data.frame(do.call(rbind,sapply(1:nrow(tox), function(x) {
  cvd = trimws(strsplit(tox[x,2], ',')[[1]])
  treat = rep(tox[x, 1], length(cvd))
  pair = cbind(treat, cvd)
  pair
})))

toxpair$cvd <- tolower(as.character(toxpair$cvd))

toxpair2 <- aggregate(toxpair$treat, list(toxpair$cvd), paste, collapse="; ")
toxpair2$x <- gsub('\n','',toxpair2$x)

write.csv(toxpair2, 'CVD and treatments.csv')
