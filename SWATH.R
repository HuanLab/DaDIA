###############################################################
#This is the main script to perform DIA SWATH feature extraction and annotation
#Sam Shen, 2020-07-28
#Copyright @ University of British Columbia
###############################################################
#Parameters for feature extraction
mass.tol <- 10 #mz tolerance in ppm
mass.const.tol <- 0 #mz tolerance in constant value: used in MS2 spectra matching
rt.tol <- 60 #rt tolerance in seconds
num.samples <- 3 ######IMPORTANT######## <------- enter how many DIA samples here #####
perform.MS2.extraction <- TRUE # whether to perform MS2 extraction and feature annotation
###############################################################
DIA.directory <- "C:/Users/User/Desktop/"
###############################################################
# Database search (dot product)
ms1.tol <- 0.01
ms2.tol <- 0.02
###############################################################
library(xcms)
library(MSnbase)
library(dplyr)
if(num.samples == 1){
  setwd(DIA.directory)
  swath_file <- list.files(pattern = ".mzXML")
  swath_data <- readMSData(swath_file, mode = "onDisk")
  swath_data <- filterEmptySpectra(swath_data)
  cwp <- CentWaveParam(ppm=10,
                       peakwidth=c(5,60),
                       mzdiff = 0.01,
                       snthresh = 6,
                       integrate = 1,
                       prefilter = c(3,100),
                       noise = 100)
  swath_data <- findChromPeaks(swath_data, param = cwp) #DIA MS1 spectra
  swath_data_filtered <- filterMsLevel(swath_data, msLevel = 1L)
  xsetSWATH <- as(swath_data_filtered, 'xcmsSet')
  
  DIAtable <- as.data.frame(DIAtable@peaks)
  colnames(DIAtable)[9] <- "intMax"
  DIAtable <- DIAtable[order(DIAtable[,1]),]
  row.names(DIAtable) <- 1:nrow(DIAtable)
  write.csv(DIAtable, file = "DIAtable.csv")
  
} else if(num.samples > 1){
  #DDA guided DIA SWATH Extraction (multi-sample)--------------------------------------------------------------
  setwd(DIA.directory)
  swath_file <- list.files(pattern = ".mzXML")
  swath_data <- readMSData(swath_file, mode = "onDisk")
  swath_data <- filterEmptySpectra(swath_data)
  cwp <- CentWaveParam(ppm=10,
                       peakwidth=c(5,60),
                       mzdiff = 0.01,
                       snthresh = 6,
                       integrate = 1,
                       prefilter = c(3,100),
                       noise = 100)
  swath_data <- findChromPeaks(swath_data, param = cwp) #DIA MS1 spectra
  swath_data_filtered <- filterMsLevel(swath_data, msLevel = 1L)
  xsetSWATH <- as(swath_data_filtered, 'xcmsSet')
  
  #ALIGNMENT  
  xsetSWATH@peaks <- xsetSWATH@peaks[order(xsetSWATH@peaks[,11]),]
  DIAtable <- as.data.frame(xsetSWATH@peaks)
  for(n in (1:length(swath_file))){
    sampleOutput <- DIAtable[DIAtable$sample == n, ]
    sampleOutput <- sampleOutput[order(sampleOutput[,1]),]
    colnames(sampleOutput)[9] <- "intMax"
    row.names(sampleOutput) <- 1:nrow(sampleOutput)
    write.csv(sampleOutput, file = paste(n,"DIAtable.csv",sep = "_"))
  }
  xsetSWATH <- group(xsetSWATH, bw = 5, minfrac = 0.5, mzwid = 0.015, minsamp = 1, max = 50)
  xsetSWATH <- retcor(xsetSWATH, method = "obiwarp", profStep = 1)
  xsetSWATH <- group(xsetSWATH, bw = 5, minfrac = 0.5, mzwid = 0.015, minsamp = 1, max = 50)
  xsetSWATH <- fillPeaks(xsetSWATH)#You must make this work as gapfilling is a standard module in any software these days
  XCMt <- data.frame(xsetSWATH@groups)
  xcmI <- groupval(xsetSWATH, value = "maxo")
  featureTable <- cbind(XCMt$mzmed, XCMt$rtmed, XCMt$rtmin, XCMt$rtmax, xcmI)
  colnames(featureTable)[1:4] <- c("mz", "RT", "RTmin", "RTmax")
  featureTable <- featureTable[order(featureTable[,1]),]
  featureTable <- cbind(featureTable, 1:nrow(featureTable))
  colnames(featureTable)[ncol(featureTable)] <- "ID"
  featureTable <- as.data.frame(featureTable)
  #Output
  write.csv(featureTable, file = "alignedDIAtable.csv")
}

##Functions----------------------------------------------------------------------------------
############
# Dot product function
dp.score <- function(x,y){
  if(nrow(x)==0 | nrow(y)==0){return(0)}
  x[,2] <- 100*x[,2]/max(x[,2])
  y[,2] <- 100*y[,2]/max(y[,2])
  alignment <- data.frame(matrix(nrow=nrow(x), ncol=3))
  alignment[,1:2] <- x[,1:2]
  y1 <- y  ##in case one row in y can be selected multiple times
  for(i in 1:nrow(x)){
    mass.diff <- abs(y1[,1] - x[i,1])
    if(min(mass.diff) <= ms2.tol){
      alignment[i,3] <- y1[mass.diff==min(mass.diff),2][1]
      y1[mass.diff==min(mass.diff),1][1] <- NA   # after matched, NA assigned
      y1 <- y1[complete.cases(y1),]
      if(is.null(nrow(y1)) ==TRUE) break
      if(nrow(y1)==0) break
    }
  }
  alignment <- alignment[complete.cases(alignment),]
  if(nrow(alignment)==0){score <- 0}
  if(nrow(alignment)>0){
    #dot product calculation
    AB <- sum(alignment[,2]*alignment[,3])
    A <- sum(x[,2]^2)
    B <- sum(y[,2]^2)
    dp.score <- AB/sqrt(A*B)
    score <- as.numeric(dp.score)
  }
  match_No <- nrow(alignment)
  return  <- c(score,match_No)
  return(return)
}  

#------------------------------------------------------------------------------------------
if(perform.MS2.extraction == TRUE){
  if(num.samples == 1){
    #DDA & DIA Extraction single sample--------------------------------------------------------------
    setwd(DIA.directory)
    swath_setting_file <- list.files(pattern = ".txt") #headers cannot contain spaces
    swath_setting <- read.table(swath_setting_file, sep = "" , header = T , nrows = 100,
                                na.strings ="", stringsAsFactors= F) 
    swath_setting$Targetmz <- with(swath_setting, (Minmz+Maxmz) / 2)
    swath_setting$offset <- with(swath_setting, (Maxmz-Minmz) / 2)
    swath_setting[1,] <- NA
    
    tempisomz <- rep(swath_setting[,5],length.out=length(swath_data))
    fData(swath_data)$isolationWindowTargetMZ <- tempisomz
    tempoffset <- rep(swath_setting[,6],length.out=length(swath_data))
    fData(swath_data)$isolationWindowLowerOffset <- tempoffset
    fData(swath_data)$isolationWindowUpperOffset <- tempoffset
    
    #MS2 spectra
    DIAtable <- cbind(DIAtable, FALSE)
    colnames(DIAtable)[ncol(DIAtable)] <- "MS2_match"
    MS2_Spectra_Table <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(MS2_Spectra_Table) <- c("ID", "PrecursorMZ", "MS2mz", "MS2int", "PeaksCount", "Source")
    cwp <- CentWaveParam(snthresh = 3, noise = 10, ppm = 10,
                         peakwidth = c(5,60))
    #performs a peak detection, separately for all spectra belonging to the same isolation 
    ##window and adds them to the chromPeaks() matrix of the result object
    swath_data <- findChromPeaksIsolationWindow(swath_data, param = cwp) #takes some time
    chromPeakData(swath_data) #lists identified peaks (both MS1 and MS2) 
    table(chromPeakData(swath_data)$isolationWindow) #count the number of chromatographic peaks identified within each isolation window
    swath_spectra <- reconstructChromPeakSpectra(swath_data, minCor = 0.2)
    
    tmpDIAtable <- as.data.frame(chromPeaks(swath_data, msLevel = 1L))
    mz.diff <- tmpDIAtable[, "mz"] * mass.tol / 1e6 
    tmpDIAtable[, "mzmin"] <- tmpDIAtable[, "mz"] - mz.diff
    tmpDIAtable[, "mzmax"] <- tmpDIAtable[, "mz"] + mz.diff
    tmpDIAtable[, "rtmin"] <- tmpDIAtable[, "rt"] - rt.tol
    tmpDIAtable[, "rtmax"] <- tmpDIAtable[, "rt"] + rt.tol
    for (i in 1:nrow(DIAtable)) {
      finalSpectra <- swath_spectra@listData[[i]]
      if(finalSpectra@peaksCount > 0){
        tmpMatch <- which((DIAtable$mz > tmpDIAtable$mzmin[i]) & 
                            (DIAtable$mz < tmpDIAtable$mzmax[i]) & 
                            (DIAtable$rt > tmpDIAtable$rtmin[i]) & 
                            (DIAtable$rt < tmpDIAtable$rtmax[i]) &
                            (DIAtable$MS2_match == FALSE))
        if(length(tmpMatch) > 0){
          for(j in 1:length(tmpMatch)){
            DIAtable$MS2_match[tmpMatch[j]] <- TRUE
            MS2_Spectra_Table[nrow(MS2_Spectra_Table) + 1,] = list(tmpMatch[j], 
                                                                   finalSpectra@precursorMz,
                                                                   paste(finalSpectra@mz, sep = ",", collapse = ","),
                                                                   paste(finalSpectra@intensity, sep = ",", collapse = ","),
                                                                   finalSpectra@peaksCount,
                                                                   "SWATH")
          }
        }
      }
    }
    
  } else if(num.samples > 1){
    #DDA & DIA Extraction multi sample--------------------------------------------------------------
    setwd(DIA.directory)
    swath_setting_file <- list.files(pattern = ".txt") #headers cannot contain spaces
    swath_setting <- read.table(swath_setting_file, sep = "" , header = T , nrows = 100,
                                na.strings ="", stringsAsFactors= F) 
    swath_setting$Targetmz <- with(swath_setting, (Minmz+Maxmz) / 2)
    swath_setting$offset <- with(swath_setting, (Maxmz-Minmz) / 2)
    swath_setting[1,] <- NA
    
    tempisomz <- rep(swath_setting[,5],length.out=sum(swath_data@featureData@data$fileIdx == 1))
    for(a in 2:length(swath_data@phenoData@data$sampleNames)){
      tisomz <- rep(swath_setting[,5],length.out=sum(swath_data@featureData@data$fileIdx == a))
      tempisomz <- c(tempisomz, tisomz)
    }
    fData(swath_data)$isolationWindowTargetMZ <- tempisomz
    tempoffset <- rep(swath_setting[,6],length.out=sum(swath_data@featureData@data$fileIdx == 1))
    for(b in 2:length(swath_data@phenoData@data$sampleNames)){
      toffset <- rep(swath_setting[,6],length.out=sum(swath_data@featureData@data$fileIdx == b))
      tempoffset <- c(tempoffset, toffset)
    }
    fData(swath_data)$isolationWindowLowerOffset <- tempoffset
    fData(swath_data)$isolationWindowUpperOffset <- tempoffset
    
    #MS2 spectra
    featureTable <- cbind(featureTable, FALSE)
    colnames(featureTable)[ncol(featureTable)] <- "MS2_match"
    MS2_Spectra_Table <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(MS2_Spectra_Table) <- c("ID", "PrecursorMZ", "MS2mz", "MS2int", "PeaksCount", "Source")
    cwp <- CentWaveParam(snthresh = 3, noise = 10, ppm = 10,
                         peakwidth = c(5,60))
    #performs a peak detection, separately for all spectra belonging to the same isolation 
    ##window and adds them to the chromPeaks() matrix of the result object
    swath_data <- findChromPeaksIsolationWindow(swath_data, param = cwp) #takes some time
    chromPeakData(swath_data) #lists identified peaks (both MS1 and MS2) 
    table(chromPeakData(swath_data)$isolationWindow) #count the number of chromatographic peaks identified within each isolation window
    swath_spectra <- reconstructChromPeakSpectra(swath_data, minCor = 0.2)
    
    tmpDIAtable <- as.data.frame(chromPeaks(swath_data, msLevel = 1L))
    mz.diff <- tmpDIAtable[, "mz"] * mass.tol / 1e6 
    tmpDIAtable[, "mzmin"] <- tmpDIAtable[, "mz"] - mz.diff
    tmpDIAtable[, "mzmax"] <- tmpDIAtable[, "mz"] + mz.diff
    tmpDIAtable[, "rtmin"] <- tmpDIAtable[, "rt"] - rt.tol
    tmpDIAtable[, "rtmax"] <- tmpDIAtable[, "rt"] + rt.tol
    for (i in 1:nrow(DIAtable)) {
      finalSpectra <- swath_spectra@listData[[i]]
      if(finalSpectra@peaksCount > 0){
        tmpMatch <- which((featureTable$mz > tmpDIAtable$mzmin[i]) & 
                            (featureTable$mz < tmpDIAtable$mzmax[i]) & 
                            (featureTable$RT > tmpDIAtable$rtmin[i]) & 
                            (featureTable$RT < tmpDIAtable$rtmax[i]) &
                            (featureTable$MS2_match == FALSE))
        if(length(tmpMatch) > 0){
          for(j in 1:length(tmpMatch)){
            featureTable$MS2_match[tmpMatch[j]] <- TRUE
            MS2_Spectra_Table[nrow(MS2_Spectra_Table) + 1,] = list(tmpMatch[j], 
                                                                   finalSpectra@precursorMz,
                                                                   paste(finalSpectra@mz, sep = ",", collapse = ","),
                                                                   paste(finalSpectra@intensity, sep = ",", collapse = ","),
                                                                   finalSpectra@peaksCount,
                                                                   "SWATH")
          }
        }
      }
    }
    
  }
}

#------------------------------------------------------------------------------------------
if(perform.MS2.extraction == TRUE){
  #############
  # load msp database
  library(CAMERA)
  library('metaMS')
  setwd(DIA.directory)
  #db.name <- 'MoNA-export-MassBank.msp'
  #database <- read.msp(db.name, only.org = FALSE,
  #                      org.set = c('C','H','N','O','P','S','F','Cl','Br','I'), noNumbers = NULL)
  database<-readRDS("database.Rds")
  MS2_Spectra_Table <- cbind(MS2_Spectra_Table, 0)
  colnames(MS2_Spectra_Table)[ncol(MS2_Spectra_Table)] <- "Annotation"
  MS2_Spectra_Table <- cbind(MS2_Spectra_Table, 0)
  colnames(MS2_Spectra_Table)[ncol(MS2_Spectra_Table)] <- "DPscore"
  
  if(num.samples == 1){
    #Metabolite Annotation single  sample--------------------------------------------------------------
    DIAtable <- cbind(DIAtable, 0)
    colnames(DIAtable)[ncol(DIAtable)] <- "Annotation"
    DIAtable <- cbind(DIAtable, 0)
    colnames(DIAtable)[ncol(DIAtable)] <- "DPscore"
    # df: premass, ms2.Q, add a new column'feature.identity'
    for(x in 1:nrow(MS2_Spectra_Table)){
      premass.Q <- MS2_Spectra_Table[x, 2]     ###query precursor ion mass
      ms2.Q <- data.frame(m.z = strsplit(MS2_Spectra_Table[x, 3], ",")[[1]],
                          int = strsplit(MS2_Spectra_Table[x, 4], ",")[[1]])  ###query ms2 input, ncol = 2, m.z & int
      ms2.Q$m.z <- as.numeric(as.character(ms2.Q$m.z))
      ms2.Q$int <- as.numeric(as.character(ms2.Q$int))
      
      output <- data.frame(matrix(ncol=3))
      colnames(output) <- c('std.name','DP.score','match_No')
      h <- 1
      for(i in 1:length(database)){
        if(is.null(database[[i]]$PrecursorMZ)==TRUE) next # no precursor mass
        
        premass.L <- database[[i]]$PrecursorMZ # database precursor
        if(abs(premass.L-premass.Q) > ms1.tol) next # precursor filter
        
        ms2.L <- as.data.frame(database[[i]]$pspectrum) # database spectrum
        name.L <- database[[i]]$Name
        
        output[h,1] <- name.L
        output[h,2] <- dp.score(ms2.Q,ms2.L)[1]
        output[h,3] <- dp.score(ms2.Q,ms2.L)[2]
        
        h <- h + 1 
      }
      output <- output[complete.cases(output),]
      
      #record dot product
      if(nrow(output > 0)){
        output <- output[order(-output[,2]),]
        dot.product <- output[1,2]
        MS2_Spectra_Table[x,8] <- dot.product
        DIAtable[MS2_Spectra_Table$ID[x], 13] <- dot.product
      }
      
      # Dp score threshold, Dp score >= 0.7 , match_No >= 6 (used in GNPS identification)
      output <- output[output[,2] >= 0.7,]
      output <- output[output[,3] >= 6,]
      
      if(nrow(output)==0) {feature.identity <- 'unknown'}
      if(nrow(output)> 0) {
        output <- output[order(-output[,2]),] # sort by scores
        feature.identity <- output[1,1] # Rank 1, std name
      } 
      MS2_Spectra_Table[x,7] <- feature.identity
      DIAtable[MS2_Spectra_Table$ID[x], 12] <- feature.identity
    }
    xsa<-xsAnnotate(xsetSWATH)    
    anF <- groupFWHM(xsa, perfwhm = 0.6)
    anI <- findIsotopes(anF, mzabs = 0.01)
    anIC <- groupCorr(anI, cor_eic_th = 0.75)
    anFA <- findAdducts(anIC, polarity="positive")
    peaklist <- getPeaklist(anFA)
    peaklist <- peaklist[order(peaklist$mz),]
    DIAtable <- cbind(DIAtable, peaklist$isotopes)
    colnames(DIAtable)[ncol(DIAtable)] <- "Isotopes"
    DIAtable <- cbind(DIAtable, peaklist$adduct)
    colnames(DIAtable)[ncol(DIAtable)] <- "Adduct"
    DIAtable <- cbind(DIAtable, as.numeric(peaklist$pcgroup))
    colnames(DIAtable)[ncol(DIAtable)] <- "pcgroup"
    write.csv(DIAtable, file = "annotated_output.csv")
    
  } else if(num.samples > 1){
    #Metabolite Annotation multi  sample--------------------------------------------------------------
    featureTable <- cbind(featureTable, 0)
    colnames(featureTable)[ncol(featureTable)] <- "Annotation"
    featureTable <- cbind(featureTable, 0)
    colnames(featureTable)[ncol(featureTable)] <- "DPscore"
    # df: premass, ms2.Q, add a new column'feature.identity'
    for(x in 1:nrow(MS2_Spectra_Table)){
      premass.Q <- MS2_Spectra_Table[x, 2]     ###query precursor ion mass
      ms2.Q <- data.frame(m.z = strsplit(MS2_Spectra_Table[x, 3], ",")[[1]],
                          int = strsplit(MS2_Spectra_Table[x, 4], ",")[[1]])  ###query ms2 input, ncol = 2, m.z & int
      ms2.Q$m.z <- as.numeric(as.character(ms2.Q$m.z))
      ms2.Q$int <- as.numeric(as.character(ms2.Q$int))
      
      output <- data.frame(matrix(ncol=3))
      colnames(output) <- c('std.name','DP.score','match_No')
      h <- 1
      for(i in 1:length(database)){
        if(is.null(database[[i]]$PrecursorMZ)==TRUE) next # no precursor mass
        
        premass.L <- database[[i]]$PrecursorMZ # database precursor
        if(abs(premass.L-premass.Q) > ms1.tol) next # precursor filter
        
        ms2.L <- as.data.frame(database[[i]]$pspectrum) # database spectrum
        name.L <- database[[i]]$Name
        
        output[h,1] <- name.L
        output[h,2] <- dp.score(ms2.Q,ms2.L)[1]
        output[h,3] <- dp.score(ms2.Q,ms2.L)[2]
        
        h <- h + 1 
      }
      output <- output[complete.cases(output),]
      
      #record dot product
      if(nrow(output > 0)){
        output <- output[order(-output[,2]),]
        dot.product <- output[1,2]
        MS2_Spectra_Table[x,8] <- dot.product
        featureTable[MS2_Spectra_Table$ID[x], length(dda_file) + 8] <- dot.product
      }
      
      # Dp score threshold, Dp score >= 0.7 , match_No >= 6 (used in GNPS identification)
      output <- output[output[,2] >= 0.7,]
      output <- output[output[,3] >= 6,]
      
      if(nrow(output)==0) {feature.identity <- 'unknown'}
      if(nrow(output)> 0) {
        output <- output[order(-output[,2]),] # sort by scores
        feature.identity <- output[1,1] # Rank 1, std name
      } 
      MS2_Spectra_Table[x,7] <- feature.identity
      featureTable[MS2_Spectra_Table$ID[x], length(swath_file) + 7] <- feature.identity
    }
    xsa<-xsAnnotate(xsetSWATH)    
    anF <- groupFWHM(xsa, perfwhm = 0.6)
    anI <- findIsotopes(anF, mzabs = 0.01)
    anIC <- groupCorr(anI, cor_eic_th = 0.75)
    anFA <- findAdducts(anIC, polarity="positive")
    peaklist <- getPeaklist(anFA)
    peaklist <- peaklist[order(peaklist$mz),]
    featureTable <- cbind(featureTable, peaklist$isotopes)
    colnames(featureTable)[ncol(featureTable)] <- "Isotopes"
    featureTable <- cbind(featureTable, peaklist$adduct)
    colnames(featureTable)[ncol(featureTable)] <- "Adduct"
    featureTable <- cbind(featureTable, as.numeric(peaklist$pcgroup))
    colnames(featureTable)[ncol(featureTable)] <- "pcgroup"
    write.csv(featureTable, file = "annotated_output.csv")
  }
}
