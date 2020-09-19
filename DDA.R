###############################################################
#This is the main script to perform DDA feature extraction and annotation
#Sam Shen 2020-08-06
#Copyright @ University of British Columbia
###############################################################
#Parameters for feature extraction
mass.tol <- 10 #mz tolerance in ppm
mass.const.tol <- 0 #mz tolerance in constant value: used in MS2 spectra matching
rt.tol <- 60 #rt tolerance in seconds
num.samples <- 14 ######IMPORTANT########
plot.DDA <- FALSE
plot.DDA.mztol <- 0.1
plot.DDA.rttol <- 30
perform.MS2.extraction <- TRUE # whether to perform MS2 extraction and feature annotation
###############################################################
DDA.directory <- "C:/Users/User/Desktop/SAM DONT TOUCH DONT DELETE/leukemiaApplication/DDA_ONLY"
###############################################################
# Database search (dot product)
ms1.tol <- 0.01
ms2.tol <- 0.02
###############################################################
library(xcms)
library(MSnbase)
library(dplyr)
if(num.samples == 1){
  setwd(DDA.directory)
  dda_file <- list.files(pattern = ".mzXML")
  dda_data <- readMSData(dda_file, mode = "onDisk")
  cwp <- CentWaveParam(ppm=10,
                       peakwidth=c(5,60),
                       mzdiff = 0.01,
                       snthresh = 6,
                       integrate = 1,
                       prefilter = c(3,100),
                       noise = 100)
  dda_data <- findChromPeaks(dda_data, param = cwp) #DDA MS1 spectra
  dda_data_filtered <- filterMsLevel(dda_data, msLevel = 1L)
  xsetDDA <- as(dda_data_filtered, 'xcmsSet')
  
  DDAtable <- as.data.frame(xsetDDA@peaks)
  colnames(DDAtable)[9] <- "intMax"
  DDAtable <- DDAtable[order(DDAtable[,1]),]
  row.names(DDAtable) <- 1:nrow(DDAtable)
  write.csv(DDAtable, file = "DDAtable.csv")
  
} else if(num.samples > 1){
  #DDA guided DIA SWATH Extraction (multi-sample)--------------------------------------------------------------
  setwd(DDA.directory)
  dda_file <- list.files(pattern = ".mzXML")
  dda_data <- readMSData(dda_file, mode = "onDisk")
  cwp <- CentWaveParam(ppm=10,
                       peakwidth=c(5,60),
                       mzdiff = 0.01,
                       snthresh = 6,
                       integrate = 1,
                       prefilter = c(3,100),
                       noise = 100)
  dda_data <- findChromPeaks(dda_data, param = cwp) #DDA MS1 spectra
  dda_data_filtered <- filterMsLevel(dda_data, msLevel = 1L)
  xsetDDA <- as(dda_data_filtered, 'xcmsSet')
  
  #ALIGNMENT  
  xsetDDA@peaks <- xsetDDA@peaks[order(xsetDDA@peaks[,11]),]
  DDAtable <- as.data.frame(xsetDDA@peaks)
  for(n in (1:length(dda_file))){
    sampleOutput <- DDAtable[DDAtable$sample == n, ]
    sampleOutput <- sampleOutput[order(sampleOutput[,1]),]
    colnames(sampleOutput)[9] <- "intMax"
    row.names(sampleOutput) <- 1:nrow(sampleOutput)
    write.csv(sampleOutput, file = paste(n,"DDAtable.csv",sep = "_"))
  }
  xsetDDA <- group(xsetDDA, bw = 5, minfrac = 0.5, mzwid = 0.015, minsamp = 1, max = 100)
  xsetDDA <- retcor(xsetDDA, method = "obiwarp", profStep = 1)
  xsetDDA <- group(xsetDDA, bw = 5, minfrac = 0.5, mzwid = 0.015, minsamp = 1, max = 100)
  xsetDDA <- fillPeaks(xsetDDA)#You must make this work as gapfilling is a standard module in any software these days
  XCMt <- data.frame(xsetDDA@groups)
  xcmI <- groupval(xsetDDA, value = "maxo")
  featureTable <- cbind(XCMt$mzmed, XCMt$rtmed, XCMt$rtmin, XCMt$rtmax, xcmI)
  colnames(featureTable)[1:4] <- c("mz", "RT", "RTmin", "RTmax")
  featureTable <- featureTable[order(featureTable[,1]),]
  featureTable <- cbind(featureTable, 1:nrow(featureTable))
  colnames(featureTable)[ncol(featureTable)] <- "ID"
  featureTable <- as.data.frame(featureTable)
  #Output
  write.csv(featureTable, file = "alignedDDAtable.csv")
  
  if(plot.DDA){
    setwd(DDA.directory)
    dir.create("DDA_EIC")
    setwd("DDA_EIC")
    plot.matrix <- featureTable[,5:(ncol(featureTable)-1)]
    xrawList <- list()
    for(n in 1:length(dda_file)){
      xrawList[n] <- xcmsRaw(filepaths(xsetDDA)[n],profstep=0)
    }
    for(k in 1:nrow(plot.matrix)){
      rt.lower.limit <- featureTable$RT[k] - plot.DDA.rttol
      rt.upper.limit <- featureTable$RT[k] + plot.DDA.rttol
      mass.lower.limit <- featureTable$mz[k] - plot.DDA.mztol
      mass.upper.limit <- featureTable$mz[k] + plot.DDA.mztol
      maxIndex <- as.numeric(which.max(plot.matrix[k,]))
      png(file = paste0(featureTable$mz[k], "_", featureTable$RT[k],".png"), width = 480, height = 480)
      eic <- plotEIC(xrawList[[maxIndex]], mzrange = c(mass.lower.limit, mass.upper.limit), 
                     rtrange = c(rt.lower.limit,rt.upper.limit))
      dev.off()
    }
  }
  setwd(DDA.directory)
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

matchMS2 <- function(x, featuretable, msLevel = 2L, expandRt = 0, expandMz = 0, ppm = 0) {
  return.type <- "MSpectra"
  pks <- featuretable
  if (ppm != 0){
    mz.diff <- pks[, "mz"] * ppm / 1e6 
  }
  if (expandMz != 0 || length(mz.diff) > 1) {
    pks[, "mzmin"] <- pks[, "mz"] - expandMz - mz.diff
    pks[, "mzmax"] <- pks[, "mz"] + expandMz + mz.diff
  }
  if (expandRt != 0) {
    pks[, "rtmin"] <- pks[, "rt"] - expandRt
    pks[, "rtmax"] <- pks[, "rt"] + expandRt
  }
  
  peak_ids <- rownames(pks)
  fromFile <- 1L
  sps <- spectra(x)
  pmz <- precursorMz(x)
  rtm <- rtime(x)
  res <- vector(mode = "list", nrow(pks))
  for (i in 1:nrow(pks)) {
    if (is.na(pks[i, "mz"]))
      next
    idx <- which(pmz >= pks[i, "mzmin"] & pmz <= pks[i, "mzmax"] &
                   rtm >= pks[i, "rtmin"] & rtm <= pks[i, "rtmax"])
    if (length(idx)) {
      res[[i]] <- lapply(sps[idx], function(z) {
        z@fromFile = fromFile
        z
      })
    }
  }
  names(res) <- peak_ids
  return(res)
}

matchMS2multi <- function(dda_sample, alignedDaDIA, sample.num, msLevel = 2L, expandRt = 0, expandMz = 0, ppm = 0) {
  return.type <- "MSpectra"
  pks <- alignedDaDIA
  pks <- cbind(pks, 0)
  colnames(pks)[ncol(pks)] <- "mzmin"
  pks <- cbind(pks, 0)
  colnames(pks)[ncol(pks)] <- "mzmax"
  if (ppm != 0){
    mz.diff <- pks[, "mz"] * ppm / 1e6 
  }
  if (expandMz != 0 || length(mz.diff) > 1) {
    pks[, "mzmin"] <- pks[, "mz"] - expandMz - mz.diff
    pks[, "mzmax"] <- pks[, "mz"] + expandMz + mz.diff
  }
  if (expandRt != 0) {
    pks[, "RTmin"] <- pks[, "RT"] - expandRt
    pks[, "RTmax"] <- pks[, "RT"] + expandRt
  }
  
  peak_ids <- as.vector(pks[,ncol(pks)-2])
  fromFile <- as.integer(sample.num)
  sps <- spectra(dda_sample)
  pmz <- precursorMz(dda_sample)
  rtm <- rtime(dda_sample)
  res <- vector(mode = "list", nrow(pks))
  for (i in 1:nrow(pks)) {
    if (is.na(pks[i, "mz"]))
      next
    idx <- which(pmz >= pks[i, "mzmin"] & pmz <= pks[i, "mzmax"] &
                   rtm >= pks[i, "RTmin"] & rtm <= pks[i, "RTmax"])
    if (length(idx)) {
      res[[i]] <- lapply(sps[idx], function(z) {
        z@fromFile = fromFile
        z
      })
    }
  }
  names(res) <- peak_ids
  return(res)
}

#------------------------------------------------------------------------------------------
if(perform.MS2.extraction == TRUE){
  if(num.samples == 1){
    #DDA & DIA Extraction single sample--------------------------------------------------------------
    #MS2 spectra
    dda_spectra <- matchMS2(dda_data, DDAtable, expandRt = rt.tol, expandMz = mass.const.tol, ppm = mass.tol)
    DDAtable <- cbind(DDAtable, FALSE)
    colnames(DDAtable)[ncol(DDAtable)] <- "MS2_match"
    MS2_Spectra_Table <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(MS2_Spectra_Table) <- c("ID", "PrecursorMZ", "MS2mz", "MS2int", "PeaksCount", "Source")
    for (i in 1:nrow(DDAtable)) {
      if(!is.null(dda_spectra[[i]])){
        tmpSpectra <- dda_spectra[[i]]
        for (j in 1:length(tmpSpectra)){
          if(tmpSpectra[[j]]@peaksCount == 0){
            tmpSpectra[[j]] <- NA
          }
        }
        tmpSpectra <- tmpSpectra[is.na(tmpSpectra)==FALSE]
        if(length(tmpSpectra) > 0){
          currInt = tmpSpectra[[1]]@precursorIntensity
          currIdx = 1
          for(k in 1:length(tmpSpectra)){
            if(tmpSpectra[[k]]@precursorIntensity > currInt){
              currIdx = k
              currInt = tmpSpectra[[k]]@precursorIntensity
            }
          }
          finalSpectra = tmpSpectra[[currIdx]]
          DDAtable$MS2_match[i] <- TRUE
          MS2_Spectra_Table[nrow(MS2_Spectra_Table) + 1,] = list(i, 
                                                                 finalSpectra@precursorMz,
                                                                 paste(finalSpectra@mz, sep = ",", collapse = ","),
                                                                 paste(finalSpectra@intensity, sep = ",", collapse = ","),
                                                                 finalSpectra@peaksCount,
                                                                 "DDA")
          #maybe another list to keep track of ID with Spectrum2 objects?
        }
      }
    }
    
  } else if(num.samples > 1){
    #DDA & DIA Extraction multi sample--------------------------------------------------------------
    #MS2 spectra
    dda_spectra <- matchMS2multi(dda_data, featureTable, 3, expandRt = rt.tol, expandMz = mass.const.tol, ppm = mass.tol)
    featureTable <- cbind(featureTable, FALSE)
    colnames(featureTable)[ncol(featureTable)] <- "MS2_match"
    MS2_Spectra_Table <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(MS2_Spectra_Table) <- c("ID", "PrecursorMZ", "MS2mz", "MS2int", "PeaksCount", "Source")
    for (i in 1:nrow(featureTable)) {
      if(!is.null(dda_spectra[[i]])){
        tmpSpectra <- dda_spectra[[i]]
        for (j in 1:length(tmpSpectra)){
          if(tmpSpectra[[j]]@peaksCount == 0){
            tmpSpectra[[j]] <- NA
          }
        }
        tmpSpectra <- tmpSpectra[is.na(tmpSpectra)==FALSE]
        if(length(tmpSpectra) > 0){
          currInt = tmpSpectra[[1]]@precursorIntensity
          currIdx = 1
          for(k in 1:length(tmpSpectra)){
            if(tmpSpectra[[k]]@precursorIntensity > currInt){
              currIdx = k
              currInt = tmpSpectra[[k]]@precursorIntensity
            }
          }
          finalSpectra = tmpSpectra[[currIdx]]
          featureTable$MS2_match[i] <- TRUE
          MS2_Spectra_Table[nrow(MS2_Spectra_Table) + 1,] = list(i, 
                                                                 finalSpectra@precursorMz,
                                                                 paste(finalSpectra@mz, sep = ",", collapse = ","),
                                                                 paste(finalSpectra@intensity, sep = ",", collapse = ","),
                                                                 finalSpectra@peaksCount,
                                                                 "DDA")
          #maybe another list to keep track of ID with Spectrum2 objects?
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
  setwd(DDA.directory)
  #db.name <- 'MoNA-export-MassBank.msp'
  #database <- read.msp(db.name, only.org = FALSE,
  #                      org.set = c('C','H','N','O','P','S','F','Cl','Br','I'), noNumbers = NULL)
  database<-readRDS("convertedLibraryNeg.Rds")
  MS2_Spectra_Table <- cbind(MS2_Spectra_Table, 0)
  colnames(MS2_Spectra_Table)[ncol(MS2_Spectra_Table)] <- "Annotation"
  MS2_Spectra_Table <- cbind(MS2_Spectra_Table, 0)
  colnames(MS2_Spectra_Table)[ncol(MS2_Spectra_Table)] <- "DPscore"
  
  if(num.samples == 1){
    #Metabolite Annotation single  sample--------------------------------------------------------------
    DDAtable <- cbind(DDAtable, 0)
    colnames(DDAtable)[ncol(DDAtable)] <- "Annotation"
    DDAtable <- cbind(DDAtable, 0)
    colnames(DDAtable)[ncol(DDAtable)] <- "DPscore"
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
        featureTable[MS2_Spectra_Table$ID[x], 13] <- dot.product
      }
      
      # Dp score threshold, Dp score >= 0.7 , match_No >= 6 (used in GNPS identification)
      output <- output[output[,2] >= 0.7,]
      output <- output[output[,3] >= 1,]
      
      if(nrow(output)==0) {feature.identity <- 'unknown'}
      if(nrow(output)> 0) {
        output <- output[order(-output[,2]),] # sort by scores
        feature.identity <- output[1,1] # Rank 1, std name
      } 
      MS2_Spectra_Table[x,7] <- feature.identity
      DDAtable[MS2_Spectra_Table$ID[x], 12] <- feature.identity
    }
    xsa<-xsAnnotate(xsetDDA)    
    anF <- groupFWHM(xsa, perfwhm = 0.6)
    anI <- findIsotopes(anF, mzabs = 0.01)
    anIC <- groupCorr(anI, cor_eic_th = 0.75)
    anFA <- findAdducts(anIC, polarity="positive")
    peaklist <- getPeaklist(anFA)
    peaklist <- peaklist[order(peaklist$mz),]
    DDAtable <- cbind(DDAtable, peaklist$isotopes)
    colnames(DDAtable)[ncol(DDAtable)] <- "Isotopes"
    DDAtable <- cbind(DDAtable, peaklist$adduct)
    colnames(DDAtable)[ncol(DDAtable)] <- "Adduct"
    DDAtable <- cbind(DDAtable, as.numeric(peaklist$pcgroup))
    colnames(DDAtable)[ncol(DDAtable)] <- "pcgroup"
    write.csv(DDAtable, file = "annotated_output.csv")
    
  } else if(num.samples > 1){
    #Metabolite Annotation multi sample--------------------------------------------------------------
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
      output <- output[output[,3] >= 1,]
      
      if(nrow(output)==0) {feature.identity <- 'unknown'}
      if(nrow(output)> 0) {
        output <- output[order(-output[,2]),] # sort by scores
        feature.identity <- output[1,1] # Rank 1, std name
      } 
      MS2_Spectra_Table[x,7] <- feature.identity
      featureTable[MS2_Spectra_Table$ID[x], length(dda_file) + 7] <- feature.identity
    }
    # xsa<-xsAnnotate(xsetDDA)    
    # anF <- groupFWHM(xsa, perfwhm = 0.6)
    # anI <- findIsotopes(anF, mzabs = 0.01)
    # anIC <- groupCorr(anI, cor_eic_th = 0.75)
    # anFA <- findAdducts(anIC, polarity="positive")
    # peaklist <- getPeaklist(anFA)
    # peaklist <- peaklist[order(peaklist$mz),]
    # featureTable <- cbind(featureTable, peaklist$isotopes)
    # colnames(featureTable)[ncol(featureTable)] <- "Isotopes"
    # featureTable <- cbind(featureTable, peaklist$adduct)
    # colnames(featureTable)[ncol(featureTable)] <- "Adduct"
    # featureTable <- cbind(featureTable, as.numeric(peaklist$pcgroup))
    # colnames(featureTable)[ncol(featureTable)] <- "pcgroup"
    write.csv(featureTable, file = "annotated_output.csv")
  }
}
