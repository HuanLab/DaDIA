# This script is to construct a msp spectral library from a prepared csv file.
# csv file: last two columns: MS2

library('metaMS')

# read in csv file
setwd("C:/Users/User/Desktop/")
file <- read.csv('template.csv',stringsAsFactors = FALSE)

uni.spectra <- unique(file[,1])

spectra <- list()
extra.info <- file[1,1:(ncol(file)-2)]
for(i in 1:length(uni.spectra)){
  ms2 <- data.frame(file[file[,1]==uni.spectra[i], (ncol(file)-1):ncol(file) ])
  spectra[[i]] <- ms2
  extra.info[i,] <- file[file[,1]==uni.spectra[i],1:(ncol(file)-2)][1,]
}

db <- construct.msp(spectra, extra.info)

write.msp(db, file = 'spectral_library.msp')
