library('metaMS')
library(doParallel)
library(foreach)
library(dplyr)
inputName <- "MSMS-Public-Neg-VS15.msp"
outputNameMSP <- "convertedLibraryNeg.msp"
outputNameRDS <- "convertedLibraryNeg.Rds"
ion.mode <- "N"
directory <- "C:/Users/User/Desktop/SAM DONT TOUCH DONT DELETE"

# Calculate the number of cores
no_cores <- detectCores() - 1
print("Using cores:")
print(no_cores)
# Initiate cluster
registerDoParallel(no_cores)

setwd(directory)
db.name <- inputName
oldDB <- read.delim(db.name,stringsAsFactors = FALSE, header = FALSE)
oldindex <- 1
dbIndex <- 1
dbList <- list()
for(i in 2:nrow(oldDB)){
  if(grepl('NAME: ', oldDB[i,], fixed = TRUE)){
    dbList[dbIndex] <- as.data.frame(oldDB[oldindex:(i-1),])
    dbIndex <- dbIndex + 1
    oldindex <- i
  }
}
r <- foreach(n = 1:length(dbList)) %dopar% {
  outDB <- data.frame(matrix(ncol = 1, nrow = 0))
  tmp <- dbList[[n]]
  if(substr(tmp[12], 1, 3) == "Num"){
    outDB <- rbind(outDB, 
                   paste("Name:", sub(".*NAME: ", "", sub("; PlaSMA.*", "", tmp[1])), sep = " "),
                   paste("InChIKey:", sub(".*INCHIKEY: ", "", tmp[7]), sep = " "),
                   paste("Precursor_type:", sub(".*PRECURSORTYPE: ", "", tmp[3]), sep = " "),
                   paste("PrecursorMZ:", sub(".*PRECURSORMZ: ", "", tmp[2]), sep = " "),
                   paste("Ion_mode:", ion.mode, sep = " "),
                   paste("Collision_energy:", "", sep = " "),
                   paste("Formula:", sub(".*FORMULA: ", "", tmp[5]), sep = " "),
                   paste("Ontology:", sub(".*ONTOLOGY: ", "", tmp[10]), sep = " "),
                   paste("Smiles:", sub(".*SMILES: ", "", tmp[6]), sep = " "),
                   paste("RetentionTime:", sub(".*RETENTIONTIME: ", "", tmp[8]), sep = " "),
                   tmp[9],
                   paste("Comment:", sub(".*COMMENT: ", "", tmp[11]), sep = " "),
                   tmp[12])
    if(as.numeric(sub(".*Num Peaks: ", "", tmp[12])) != 0){
      for(j in seq.int(13L, length(tmp)-1, 2L)){
        outDB <- rbind(outDB, paste(tmp[j], tmp[j+1],  sep = " "))
      }
      outDB <- rbind(outDB, NA)
      colnames(outDB) <- "info"
    }else if(as.numeric(sub(".*Num Peaks: ", "", tmp[12])) == 0){
      outDB <- NULL
    }
  }else if(substr(tmp[13], 1, 3) == "Num"){
    outDB <- rbind(outDB, 
                   paste("Name:", sub(".*NAME: ", "", sub("; PlaSMA.*", "", tmp[1])), sep = " "),
                   paste("InChIKey:", sub(".*INCHIKEY: ", "", tmp[6]), sep = " "),
                   paste("Precursor_type:", sub(".*PRECURSORTYPE: ", "", tmp[3]), sep = " "),
                   paste("PrecursorMZ:", sub(".*PRECURSORMZ: ", "", tmp[2]), sep = " "),
                   paste("Ion_mode:", ion.mode, sep = " "),
                   paste("Collision_energy:", "", sep = " "),
                   paste("Formula:", sub(".*FORMULA: ", "", tmp[4]), sep = " "),
                   tmp[5],
                   paste("Smiles:", sub(".*SMILES: ", "", tmp[7]), sep = " "),
                   paste("RetentionTime:", sub(".*RETENTIONTIME: ", "", tmp[8]), sep = " "),
                   tmp[9],
                   tmp[12],
                   tmp[13])
    if(as.numeric(sub(".*Num Peaks: ", "", tmp[13])) != 0){
      for(j in seq.int(14L, length(tmp)-1, 2L)){
        outDB <- rbind(outDB, paste(tmp[j], tmp[j+1],  sep = " "))
      }
      outDB <- rbind(outDB, NA)
      colnames(outDB) <- "info"
    }else if(as.numeric(sub(".*Num Peaks: ", "", tmp[13])) == 0){
      outDB <- NULL
    }
  }else if(substr(tmp[15], 1, 3) == "Num"){
    outDB <- rbind(outDB, 
                   paste("Name:", sub(".*NAME: ", "", sub("; PlaSMA.*", "", tmp[1])), sep = " "),
                   paste("InChIKey:", sub(".*INCHIKEY: ", "", tmp[7]), sep = " "),
                   paste("Precursor_type:", sub(".*PRECURSORTYPE: ", "", tmp[3]), sep = " "),
                   paste("PrecursorMZ:", sub(".*PRECURSORMZ: ", "", tmp[2]), sep = " "),
                   paste("Ion_mode:", ion.mode, sep = " "),
                   paste("Collision_energy:", "", sep = " "),
                   paste("Formula:", sub(".*FORMULA: ", "", tmp[5]), sep = " "),
                   paste("Ontology:", sub(".*ONTOLOGY: ", "", tmp[13]), sep = " "),
                   paste("Smiles:", sub(".*SMILES: ", "", tmp[6]), sep = " "),
                   paste("RetentionTime:", sub(".*RETENTIONTIME: ", "", tmp[11]), sep = " "),
                   tmp[12],
                   paste("Comment:", sub(".*COMMENT: ", "", tmp[14]), sep = " "),
                   tmp[15])
    if(as.numeric(sub(".*Num Peaks: ", "", tmp[15])) != 0){
      for(j in seq.int(16L, length(tmp)-1, 2L)){
        outDB <- rbind(outDB, paste(tmp[j], tmp[j+1],  sep = " "))
      }
      outDB <- rbind(outDB, NA)
      colnames(outDB) <- "info"
    }else if(as.numeric(sub(".*Num Peaks: ", "", tmp[15])) == 0){
      outDB <- NULL
    }
  }
  return(outDB)
}
output <- bind_rows(r)
write.table(output, file = outputNameMSP, row.names = FALSE, col.names = FALSE, 
            quote = FALSE, na = "")

#clean up the cluster
stopImplicitCluster()

db.name <- outputNameMSP
database <- read.msp(db.name)
saveRDS(database, file = outputNameRDS)

database <- read.msp(outputNameMSP, only.org = FALSE,
                     org.set = c('C','H','N','O','P','S','F','Cl','Br','I'), noNumbers = NULL)
