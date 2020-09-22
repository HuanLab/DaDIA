###############################################################
#This is the main script to perform DDA aided DIA feature extraction and annotation
#Tao Huan, Sam Shen, Jian Guo 2020-07-28
#Copyright @ University of British Columbia
###############################################################
print("Loading required packages ...")
library(xcms) 
library(MSnbase)
library(dplyr)
library(doParallel)
library(foreach)
library(BiocGenerics)
library(S4Vectors)
library(ProtGenerics)
print("Finished loading packages")
###############################################################
#Part 1: Parameters for feature extraction
DDA.directory <- "C:/Users/User/Desktop/DDA"
DIA.directory <- "C:/Users/User/Desktop/DIA"
cwpDDA <- CentWaveParam(ppm=10,
                        peakwidth=c(5,60),
                        mzdiff = 0.01,
                        snthresh = 6,
                        integrate = 1,
                        prefilter = c(3,100),
                        noise = 100) #XCMS parameters for DDA feature extraction
cwpDIA <- CentWaveParam(ppm=10,
                        peakwidth=c(5,60),
                        mzdiff = 0.01,
                        snthresh = 6,
                        integrate = 1,
                        prefilter = c(3,100),
                        noise = 100) #XCMS parameters for DIA feature extraction
mass.tol <- 10 #mz tolerance in ppm: used in feature dereplication and MS2 matching
mass.const.tol <- 0.05 #mz tolerance in constant value: used in feature rescue
rt.tol <- 60 #rt tolerance in seconds
num.samples <- 11 #enter how many DIA samples here
plot.DaDIA <- TRUE #plot DaDIA features
plot.DaDIA.mztol <- 0.5 #DaDIA feature plotting mz window width
plot.DaDIA.rttol <- 30 #DaDIA feature plotting rt window width
#Parameters for alignment
bw <- 5
minfrac <- 0.5 #retention time tolerace in minutes 
mzwid <- 0.025 #mass tolerance
max <- 100
quantitative.method <- "maxo"
# "maxo" = peak height
# "into" = peak area
###############################################################
#Part 2: Parameters for database search (dot product)
feature.annotation <- TRUE #annotate DaDIA features
db.name <- "convertedLibrarPos.Rds" #annotation library name
RDS <- TRUE
# "TRUE" = database is in RDS format
# "FALSE" = database is in MSP format
ms1.tol <- 0.01 #dot product calculation ms1 tolerance
ms2.tol <- 0.02 #dot product calculation ms2 tolerance
dot.product.threshold <- 0.1 #dot product annotation threshold
match.number.threshold <- 1 #annotation match number threshold
adduct_isotope.annotation <- TRUE #perform CAMERA annotation
export.mgf <- TRUE #export individual MS2 spectra as .mgf
combine.mgf <- TRUE #combine all exported .mgf files
###############################################################
DIA.unique <- 1 #do not change
DDA.aid <- 2 #do not change
###############################################################

# Calculate the number of cores
no_cores <- detectCores() - 1
print("Using cores:")
print(no_cores)
# Initiate cluster
registerDoParallel(no_cores)

start_time <- Sys.time()
if(num.samples == 1){
  setwd(DDA.directory)
  dda_file <- list.files(pattern = ".mzXML")
  dda_data <- readMSData(dda_file, mode = "onDisk")
  dda_data <- findChromPeaks(dda_data, param = cwpDDA) #DDA feature extraction
  dda_data_filtered <- filterMsLevel(dda_data, msLevel = 1L)
  xsetDDA <- as(dda_data_filtered, 'xcmsSet')
  
  setwd(DIA.directory)
  swath_file <- list.files(pattern = ".mzXML")
  swath_data <- readMSData(swath_file, mode = "onDisk")
  swath_data <- filterEmptySpectra(swath_data)
  swath_data <- findChromPeaks(swath_data, param = cwpDIA) #DIA feature extraction
  swath_data_filtered <- filterMsLevel(swath_data, msLevel = 1L)
  xsetSWATH <- as(swath_data_filtered, 'xcmsSet')
  
  #DDA guided DIA SWATH Extraction
  DIAtable <- as.data.frame(xsetSWATH@peaks) #DIA features
  DDAtable <- as.data.frame(xsetDDA@peaks) #DDA features
  colnames(DIAtable)[ncol(DIAtable)] <- "DDA_DIA"
  DIAtable$DDA_DIA <- DIA.unique
  colnames(DDAtable)[ncol(DDAtable)] <- "present_in_DIA"
  DDAtable$present_in_DIA <- FALSE 
  
  #label DDA feature present in DIA as TRUE
  for(i in 1:nrow(DDAtable)){
    mass.lower.limit <- DDAtable$mz[i] * (1 - mass.tol * 1e-6)
    mass.upper.limit <- DDAtable$mz[i] * (1 + mass.tol * 1e-6)
    rt.lower.limit <- DDAtable$rt[i] - rt.tol
    rt.upper.limit <- DDAtable$rt[i] + rt.tol
    short.list <- DIAtable[DIAtable$mz >= mass.lower.limit & DIAtable$mz <= mass.upper.limit,]
    short.list <- short.list[short.list$rt >= rt.lower.limit & short.list$rt <= rt.upper.limit,]
    if(nrow (short.list) > 0){
      DDAtable$present_in_DIA[i] <- TRUE
    }
  }
  #Rescue features in DIA guided by DDA
  xrawSWATH <- xcmsRaw(swath_file, profstep=0)
  uniqueDDAtable <- DDAtable[DDAtable$present_in_DIA == FALSE, ]
  is.inDIA <- logical(length = nrow(uniqueDDAtable)) #vector with length of row numbers of unique DDA table
  inDIA.matrix <- data.frame(matrix(nrow = nrow(uniqueDDAtable), ncol = ncol(DIAtable)))
  colnames(inDIA.matrix) <- colnames(DIAtable)
  for (j in 1:nrow(uniqueDDAtable)){
    mass.lower.limit <- uniqueDDAtable$mz[j] - mass.const.tol
    mass.upper.limit <- uniqueDDAtable$mz[j] + mass.const.tol
    rt.lower.limit <- uniqueDDAtable$rt[j] - rt.tol
    rt.upper.limit <- uniqueDDAtable$rt[j] + rt.tol
    # filter the features out of the retention time range
    if(rt.lower.limit > tail(xrawSWATH@scantime, n=1) | rt.upper.limit > tail(xrawSWATH@scantime, n=1)) next  
    if(rt.lower.limit < xrawSWATH@scantime[1]+1){
      rt.lower.limit <- xrawSWATH@scantime[1]+1
    }
    if(rt.lower.limit < 1){
      rt.lower.limit <- 1
    }
    if(rt.upper.limit > tail(xrawSWATH@scantime, n=1)){
      rt.upper.limit <- tail(xrawSWATH@scantime, n=1) -1
    }
    # filter the features out of the m/z range
    if(mass.lower.limit < xrawSWATH@mzrange[1]) next()
    if(mass.upper.limit > xrawSWATH@mzrange[2]) next()
    mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
    RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
    eeic <- getEIC(xrawSWATH, mzrange=mzRange, rtrange=RTRange) #extracted EIC object
    
    eic.matrix <- as.data.frame(eeic@eic[["xcmsRaw"]][[1]][,"intensity"])
    eic.matrix <- eic.matrix[is.na(eic.matrix[,1])==FALSE,]
    peak.int <- max(eic.matrix) #find the max intensity in the EIC
    
    if(is.na(peak.int)) next
    is.inDIA[j] <- TRUE
    #Put rescued features in DaDIA featureTable
    inDIA.matrix[j,1]  <- uniqueDDAtable$mz[j]
    inDIA.matrix[j,2]  <- uniqueDDAtable$mz[j]
    inDIA.matrix[j,3]  <- uniqueDDAtable$mz[j]
    inDIA.matrix[j,4]  <- uniqueDDAtable$rt[j]
    inDIA.matrix[j,5]  <- uniqueDDAtable$rt[j]
    inDIA.matrix[j,6]  <- uniqueDDAtable$rt[j]
    inDIA.matrix[j,9]  <- peak.int
    inDIA.matrix[j,11] <- DDA.aid
  }
  inDIA.matrix <- inDIA.matrix[is.na(inDIA.matrix$mz)==FALSE,]
  xsetSWATH@peaks <- rbind(xsetSWATH@peaks, as.matrix(inDIA.matrix))
  DaDIAtable <- rbind(DIAtable, inDIA.matrix)
  colnames(DaDIAtable)[9] <- "intMax"
  DaDIAtable <- DaDIAtable[order(DaDIAtable[,1]),]
  row.names(DaDIAtable) <- 1:nrow(DaDIAtable)
  write.csv(DaDIAtable, file = "DaDIAtable.csv")
  
} else if(num.samples > 1){
  #DDA guided DIA SWATH Extraction (multi-sample)
  print("Extracting DDA features ...")
  setwd(DDA.directory)
  dda_file <- list.files(pattern = ".mzXML")
  dda_data <- readMSData(dda_file, mode = "onDisk")
  dda_data <- findChromPeaks(dda_data, param = cwpDDA, SnowParam()) #DDA feature extraction
  dda_data_filtered <- filterMsLevel(dda_data, msLevel = 1L)
  xsetDDA <- as(dda_data_filtered, 'xcmsSet')
  print("Finished DDA feature extraction")
  print(Sys.time() - start_time)
  
  print("Extracting SWATH features ...")
  setwd(DIA.directory)
  swath_file <- list.files(pattern = ".mzXML")
  swath_data <- readMSData(swath_file, mode = "onDisk")
  swath_data <- filterEmptySpectra(swath_data)
  swath_data <- findChromPeaks(swath_data, param = cwpDIA, SnowParam()) #DIA feature extraction
  swath_data_filtered <- filterMsLevel(swath_data, msLevel = 1L)
  xsetSWATH <- as(swath_data_filtered, 'xcmsSet')
  print("Finished SWATH feature extraction")
  print(Sys.time() - start_time)
  
  xsetSWATH@peaks <- cbind(xsetSWATH@peaks, DIA.unique)
  colnames(xsetSWATH@peaks)[ncol(xsetSWATH@peaks)] <- "DDA_DIA"
  DIAtable <- as.data.frame(xsetSWATH@peaks) #generate data frame with DIA features
  
  print("Generating dereplicated DDA feature list ...")
  dereplicatedDDAtable <- data.frame(matrix(ncol = 11, nrow = 0)) #generate data frame with dereplicated DDA features
  rawDDAtable <- as.data.frame(xsetDDA@peaks)
  colnames(dereplicatedDDAtable) <- colnames(rawDDAtable)
  for(m in (1:nrow(xsetDDA@peaks))) {
    mass.lower.limit <- rawDDAtable$mz[m] * (1 - mass.tol * 1e-6)
    mass.upper.limit <- rawDDAtable$mz[m] * (1 + mass.tol * 1e-6)
    rt.lower.limit <- rawDDAtable$rt[m] - rt.tol
    rt.upper.limit <- rawDDAtable$rt[m] + rt.tol
    temp <- dereplicatedDDAtable[dereplicatedDDAtable$mz >= mass.lower.limit & dereplicatedDDAtable$mz <= mass.upper.limit,]
    temp <- temp[temp$rt >= rt.lower.limit & temp$rt <= rt.upper.limit,]
    if(nrow(temp) == 0) {
      dereplicatedDDAtable[nrow(dereplicatedDDAtable) + 1,] = rawDDAtable[m,]
    }
  }
  dereplicatedDDAtable <- cbind(dereplicatedDDAtable, FALSE)
  colnames(dereplicatedDDAtable)[ncol(dereplicatedDDAtable)] <- "present_in_DIA"
  print("Finished generating dereplicated DDA feature list")
  print(Sys.time() - start_time)
  
  print("Rescuing SWATH features ...")
  rescue <- foreach(n = (1:length(swath_file)), .packages = c("xcms", "dplyr")) %dopar% {
    #label DDA feature present in DIA as TRUE
    print(n)
    dereplicatedDDAtable$present_in_DIA <- FALSE
    ssDIA <- DIAtable[DIAtable$sample == n, ] # generate sample stratified DIA data frame
    for(i in 1:nrow(dereplicatedDDAtable)){
      mass.lower.limit <- dereplicatedDDAtable$mz[i] * (1 - mass.tol * 1e-6)
      mass.upper.limit <- dereplicatedDDAtable$mz[i] * (1 + mass.tol * 1e-6)
      rt.lower.limit <- dereplicatedDDAtable$rt[i] - rt.tol
      rt.upper.limit <- dereplicatedDDAtable$rt[i] + rt.tol
      short.list <- ssDIA[ssDIA$mz >= mass.lower.limit & ssDIA$mz <= mass.upper.limit,]
      short.list <- short.list[short.list$rt >= rt.lower.limit & short.list$rt <= rt.upper.limit,]
      if(nrow (short.list) > 0){
        dereplicatedDDAtable$present_in_DIA[i] <- TRUE
      }
    }
    
    #rescue features in DIA guided by DDA
    xrawSWATH <- xcmsRaw(filepaths(xsetSWATH)[n],profstep=0)
    uniqueDDAtable <- dereplicatedDDAtable[dereplicatedDDAtable$present_in_DIA == FALSE, ]
    is.inDIA <- logical(length = nrow(uniqueDDAtable))
    inDIA.matrix <- data.frame(matrix(nrow = nrow(uniqueDDAtable), ncol = ncol(DIAtable)))
    colnames(inDIA.matrix) <- colnames(DIAtable)
    for (j in 1:nrow(uniqueDDAtable)){
      mass.lower.limit <- uniqueDDAtable$mz[j] - mass.const.tol
      mass.upper.limit <- uniqueDDAtable$mz[j] + mass.const.tol
      rt.lower.limit <- uniqueDDAtable$rt[j] - rt.tol
      rt.upper.limit <- uniqueDDAtable$rt[j] + rt.tol
      # filter the features out of the retention time range
      if(rt.lower.limit > tail(xrawSWATH@scantime, n=1) | rt.upper.limit > tail(xrawSWATH@scantime, n=1)) next  
      if(rt.lower.limit < xrawSWATH@scantime[1]+1){
        rt.lower.limit <- xrawSWATH@scantime[1]+1
      }
      if(rt.lower.limit < 1){
        rt.lower.limit <- 1
      }
      if(rt.upper.limit > tail(xrawSWATH@scantime, n=1)){
        rt.upper.limit <- tail(xrawSWATH@scantime, n=1) -1
      }
      # filter the features out of the m/z range
      if(mass.lower.limit < xrawSWATH@mzrange[1]) next()
      if(mass.upper.limit > xrawSWATH@mzrange[2]) next()
      mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
      RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
      eeic <- getEIC(xrawSWATH, mzrange=mzRange, rtrange=RTRange)
      
      eic.matrix <- as.data.frame(eeic@eic[["xcmsRaw"]][[1]][,"intensity"])
      eic.matrix <- eic.matrix[is.na(eic.matrix[,1])==FALSE,]
      peak.int <- max(eic.matrix)
      if(is.na(peak.int)) next
      is.inDIA[j] <- TRUE
      #Put rescued features in DaDIA table
      inDIA.matrix[j,1]  <- uniqueDDAtable$mz[j]
      inDIA.matrix[j,2]  <- uniqueDDAtable$mz[j]
      inDIA.matrix[j,3]  <- uniqueDDAtable$mz[j]
      inDIA.matrix[j,4]  <- uniqueDDAtable$rt[j]
      inDIA.matrix[j,5]  <- uniqueDDAtable$rt[j]
      inDIA.matrix[j,6]  <- uniqueDDAtable$rt[j]
      inDIA.matrix[j,9]  <- peak.int
      inDIA.matrix[j,11] <- n
      inDIA.matrix[j,12] <- DDA.aid
    }
    inDIA.matrix <- inDIA.matrix[is.na(inDIA.matrix$mz)==FALSE,]
    return(inDIA.matrix)
  }
  
  for (t in 1:length(rescue)) {
    xsetSWATH@peaks <- rbind(xsetSWATH@peaks, as.matrix(rescue[[t]]))
  }
  print("Finished SWATH feature rescue")
  print(Sys.time() - start_time)
  
  xsetSWATH@peaks <- xsetSWATH@peaks[order(xsetSWATH@peaks[,11]),]
  DaDIAtable <- as.data.frame(xsetSWATH@peaks)
  for(n in (1:length(swath_file))){
    setwd(DIA.directory)
    sampleOutput <- DaDIAtable[DaDIAtable$sample == n, ]
    sampleOutput <- sampleOutput[order(sampleOutput[,1]),]
    colnames(sampleOutput)[9] <- "intMax"
    row.names(sampleOutput) <- 1:nrow(sampleOutput)
    write.csv(sampleOutput, file = paste(n,"DaDIAtable.csv",sep = "_"))
  }
  print("Output individual sample data")
  
  #ALIGNMENT 
  print("Aligning sample features ...")
  xsetSWATH <- group(xsetSWATH, bw = bw, minfrac = minfrac, mzwid = mzwid, minsamp = 1, max = max)
  xsetSWATH <- retcor(xsetSWATH, method = "obiwarp", profStep = 1)
  xsetSWATH <- group(xsetSWATH, bw = bw, minfrac = minfrac, mzwid = mzwid, minsamp = 1, max = max)
  xsetSWATH <- fillPeaks(xsetSWATH)
  XCMt <- data.frame(xsetSWATH@groups)
  xcmI <- groupval(xsetSWATH, value = quantitative.method)
  featureTable <- cbind(XCMt$mzmed, XCMt$rtmed, XCMt$rtmin, XCMt$rtmax, xcmI)
  colnames(featureTable)[1:4] <- c("mz", "RT", "RTmin", "RTmax")
  featureTable <- featureTable[order(featureTable[,1]),]
  featureTable <- cbind(featureTable, 1:nrow(featureTable))
  colnames(featureTable)[ncol(featureTable)] <- "ID"
  featureTable <- as.data.frame(featureTable)
  #Output
  write.csv(featureTable, file = "alignedDaDIAtable.csv")
  print("Alignment finished")
  print(Sys.time() - start_time)
  
  if(plot.DaDIA){
    print("Plotting DaDIA features ...")
    setwd(DIA.directory)
    dir.create("DaDIA_EIC")
    setwd("DaDIA_EIC")
    plot.matrix <- featureTable[,5:(ncol(featureTable)-1)]
    xrawList <- list()
    for(n in 1:length(swath_file)){
      xrawList[n] <- xcmsRaw(filepaths(xsetSWATH)[n],profstep=0)
    }
    for(k in 1:nrow(plot.matrix)){
      rt.lower.limit <- featureTable$RT[k] - plot.DaDIA.rttol
      rt.upper.limit <- featureTable$RT[k] + plot.DaDIA.rttol
      mass.lower.limit <- featureTable$mz[k] - plot.DaDIA.mztol
      mass.upper.limit <- featureTable$mz[k] + plot.DaDIA.mztol
      maxIndex <- as.numeric(which.max(plot.matrix[k,]))
      png(file = paste0(featureTable$mz[k],"_", featureTable$RT[k],".png"), width = 480, height = 480)
      eic <- plotEIC(xrawList[[maxIndex]], mzrange = c(mass.lower.limit, mass.upper.limit), 
                     rtrange = c(rt.lower.limit,rt.upper.limit))
      dev.off()
    }
    print("Finished DaDIA feature plotting")
    print(Sys.time() - start_time)
  }
  setwd(DIA.directory)
}

#Functions--------------------------------------------------------------------------------------------
#Dot product function
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
#Single sample MS2 matching
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
#Multi sample MS2 matching
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

#MS2 deconvolution-----------------------------------------------------------------------------------
if(feature.annotation == TRUE){
  if(num.samples == 1){
    #DDA & DIA Extraction single sample
    setwd(DIA.directory)
    swath_setting_file <- list.files(pattern = ".txt") 
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
    
    # This code is for DIA-AIF
    # length.pre.num <- length(swath_data@featureData@data[["precursorScanNum"]])-1
    # precursor.peaknum <- 0:length.pre.num
    # precursor.peaknum[seq(1,length(precursor.peaknum),2)] <- NA
    # swath_data@featureData@data[["precursorScanNum"]] <- precursor.peaknum
    # swath_data@featureData@data[["precursorMZ"]] <- rep(swath_setting[,4],length.out=length(swath_data))
    # precursor.charge <- rep(c(0), times = length.pre.num + 1L)
    # precursor.charge[seq(1,length(precursor.charge),2)] <- NA
    # swath_data@featureData@data[["precursorCharge"]] <- precursor.charge
    
    
    fData(swath_data)[, c("isolationWindowTargetMZ",
                          "isolationWindowLowerOffset",
                          "isolationWindowUpperOffset",
                          "msLevel", "retentionTime")]
    #view isolation window of mz
    isolationWindowLowerMz(swath_data)
    isolationWindowUpperMz(swath_data)
    #list the number of spectra that are recorded in each pocket/isolation window
    table(isolationWindowTargetMz(swath_data))
    
    #MS2 spectral assignment from DDA
    dda_spectra <- matchMS2(dda_data, DaDIAtable, expandRt = rt.tol, expandMz = mass.const.tol, ppm = mass.tol)
    DaDIAtable <- cbind(DaDIAtable, FALSE)
    colnames(DaDIAtable)[ncol(DaDIAtable)] <- "MS2_match"
    MS2_Spectra_Table <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(MS2_Spectra_Table) <- c("ID", "PrecursorMZ", "MS2mz", "MS2int", "PeaksCount", "Source")
    for (i in 1:nrow(DaDIAtable)) {
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
          DaDIAtable$MS2_match[i] <- TRUE
          MS2_Spectra_Table[nrow(MS2_Spectra_Table) + 1,] = list(i, 
                                                                 finalSpectra@precursorMz,
                                                                 paste(finalSpectra@mz, sep = ",", collapse = ","),
                                                                 paste(finalSpectra@intensity, sep = ",", collapse = ","),
                                                                 finalSpectra@peaksCount,
                                                                 "DDA")
        }
      }
    }
    
    cwp <- CentWaveParam(snthresh = 3, noise = 10, ppm = 10,
                         peakwidth = c(5,60))
    #performs a peak detection, separately for all spectra belonging to the same isolation window and adds them to the chromPeaks() matrix of the result object
    swath_data <- findChromPeaksIsolationWindow(swath_data, param = cwp) 
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
        tmpMatch <- which((DaDIAtable$mz > tmpDIAtable$mzmin[i]) & 
                            (DaDIAtable$mz < tmpDIAtable$mzmax[i]) & 
                            (DaDIAtable$rt > tmpDIAtable$rtmin[i]) & 
                            (DaDIAtable$rt < tmpDIAtable$rtmax[i]) &
                            (DaDIAtable$MS2_match == FALSE))
        if(length(tmpMatch) > 0){
          for(j in 1:length(tmpMatch)){
            DaDIAtable$MS2_match[tmpMatch[j]] <- TRUE
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
    #DDA & DIA Extraction multi sample
    print("Extracting SWATH pocket information ...")
    setwd(DIA.directory)
    swath_setting_file <- list.files(pattern = ".txt") 
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
    # #Visualizatoin of isolation window of mz
    # fData(swath_data)[, c("isolationWindowTargetMZ",
    #                       "isolationWindowLowerOffset",
    #                       "isolationWindowUpperOffset",
    #                       "msLevel", "retentionTime")]
    # #view isolation window of mz
    # isolationWindowLowerMz(swath_data)
    # isolationWindowUpperMz(swath_data)
    # #list the number of spectra that are recorded in each isolation window
    # table(isolationWindowTargetMz(swath_data))
    
    ##MS2 spectral assignment from DDA
    print("Matching MS2 spectra using DDA MS2 scans ...")
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
        }
      }
    }
    print("Finished DDA MS2 spectra matching")
    print(Sys.time() - start_time)
    
    print("Deconvoluting SWATH MS2 scans ...")
    cwp <- CentWaveParam(snthresh = 3, noise = 10, ppm = 10,
                         peakwidth = c(5,60))
    #performs a peak detection, separately for all spectra belonging to the same isolation window and adds them to the chromPeaks() matrix of the result object
    swath_data <- findChromPeaksIsolationWindow(swath_data, param = cwp) 
    chromPeakData(swath_data) #lists identified peaks (both MS1 and MS2) 
    table(chromPeakData(swath_data)$isolationWindow) #count the number of chromatographic peaks identified within each isolation window
    swath_spectra <- reconstructChromPeakSpectra(swath_data, minCor = 0.2, BPPARAM = SnowParam())
    print("Finished deconvoluting SWATH MS2 scans")
    print(Sys.time() - start_time)
    
    print("Matching MS2 spectra using SWATH MS2 scans ...")
    tmpDIAtable <- as.data.frame(chromPeaks(swath_data, msLevel = 1L))
    mz.diff <- tmpDIAtable[, "mz"] * mass.tol / 1e6 
    tmpDIAtable[, "mzmin"] <- tmpDIAtable[, "mz"] - mz.diff
    tmpDIAtable[, "mzmax"] <- tmpDIAtable[, "mz"] + mz.diff
    tmpDIAtable[, "rtmin"] <- tmpDIAtable[, "rt"] - rt.tol
    tmpDIAtable[, "rtmax"] <- tmpDIAtable[, "rt"] + rt.tol
    combSpectra <- list()
    combSpectra[nrow(featureTable) + 1] <- NULL 
    for (i in 1:nrow(DIAtable)) {
      currSpectra <- swath_spectra@listData[[i]]
      if(currSpectra@peaksCount > 0){
        tmpMatch <- which((featureTable$mz > tmpDIAtable$mzmin[i]) & 
                            (featureTable$mz < tmpDIAtable$mzmax[i]) & 
                            (featureTable$RT > tmpDIAtable$rtmin[i]) & 
                            (featureTable$RT < tmpDIAtable$rtmax[i]) &
                            (featureTable$MS2_match == FALSE))
        if(length(tmpMatch) > 0){
          for(j in 1:length(tmpMatch)){
            if(is.null(combSpectra[[tmpMatch[j]]])){
              combSpectra[[tmpMatch[j]]] <- list(currSpectra)
            }else{
              combSpectra[[tmpMatch[j]]] <- c(combSpectra[[tmpMatch[j]]], currSpectra)
            }
          }
        }
      }
    }
    print("Finished SWATH MS2 spectra matching")
    print(Sys.time() - start_time)
    print("Combining SWATH MS2 spectra cross different samples ...")
    combined_Spectra <- foreach(t = 1:nrow(featureTable), .packages = c("xcms", "MSnbase")) %dopar% {
      if(is.null(combSpectra[[t]])){
        return(NULL)
      }else{
        if(length(combSpectra[[t]]) > 1){
          combined <- consensusSpectrum(
            combSpectra[[t]],
            mzd = 0,
            minProp = 0.1,
            intensityFun = stats::median,
            mzFun = stats::median,
            ppm = 20,
            weighted = FALSE)
          return(combined)
        }else{
          combined <- combSpectra[[t]][[1]]
          return(combined)
        }
      }
    }
    print("Finished SWATH MS2 spectra consensus")
    print(Sys.time() - start_time)
    print("Assigning SWATH MS2 spectra to features ...")
    for (i in 1:nrow(featureTable)) {
      if(is.null(combined_Spectra[[i]]) == FALSE){
        finalSpectra <- combined_Spectra[[i]]
        featureTable$MS2_match[i] <- TRUE
        MS2_Spectra_Table[nrow(MS2_Spectra_Table) + 1,] = list(i, 
                                                               finalSpectra@precursorMz,
                                                               paste(finalSpectra@mz, sep = ",", collapse = ","),
                                                               paste(finalSpectra@intensity, sep = ",", collapse = ","),
                                                               finalSpectra@peaksCount,
                                                               "SWATH")
      }
    }
    print("Finished assigning SWATH MS2 spectra")
    print(Sys.time() - start_time)
  }
}

#Metabolites annotation
if(feature.annotation == TRUE){
  # load msp database
  print("Loading annotation library ...")
  library(CAMERA)
  library('metaMS')
  setwd(DIA.directory)
  if(RDS){
    database<-readRDS(db.name)
  }else{
    database <- read.msp(db.name, only.org = FALSE,
                         org.set = c('C','H','N','O','P','S','F','Cl','Br','I'), noNumbers = NULL)
  }
  print("Finished library import")
  print(Sys.time() - start_time)
  MS2_Spectra_Table <- cbind(MS2_Spectra_Table, 0)
  colnames(MS2_Spectra_Table)[ncol(MS2_Spectra_Table)] <- "Annotation"
  MS2_Spectra_Table <- cbind(MS2_Spectra_Table, 0)
  colnames(MS2_Spectra_Table)[ncol(MS2_Spectra_Table)] <- "DPscore"
  
  if(num.samples == 1){
    #Metabolites annotation single sample
    DaDIAtable <- cbind(DaDIAtable, 0)
    colnames(DaDIAtable)[ncol(DaDIAtable)] <- "Annotation"
    DaDIAtable <- cbind(DaDIAtable, 0)
    colnames(DaDIAtable)[ncol(DaDIAtable)] <- "DPscore"
    for(x in 1:nrow(MS2_Spectra_Table)){
      premass.Q <- MS2_Spectra_Table[x, 2] #query precursor ion mass
      ms2.Q <- data.frame(m.z = strsplit(MS2_Spectra_Table[x, 3], ",")[[1]],
                          int = strsplit(MS2_Spectra_Table[x, 4], ",")[[1]]) #query MS2 input, ncol = 2, m.z & int
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
        DaDIAtable[MS2_Spectra_Table$ID[x], 14] <- dot.product
      }
      
      # Dp score threshold, Dp score >= 0.7 , match_No >= 6 (used in GNPS identification)
      output <- output[output[,2] >= dot.product.threshold,]
      output <- output[output[,3] >= match.number.threshold,]
      
      if(nrow(output)==0) {feature.identity <- 'unknown'}
      if(nrow(output)> 0) {
        output <- output[order(-output[,2]),] # sort by scores
        feature.identity <- output[1,1] # Rank 1, std name
      } 
      MS2_Spectra_Table[x,7] <- feature.identity
      DaDIAtable[MS2_Spectra_Table$ID[x], 13] <- feature.identity
    }
    if(adduct_isotope.annotation) {
      xsa<-xsAnnotate(xsetSWATH)    
      anF <- groupFWHM(xsa, perfwhm = 0.6)
      anI <- findIsotopes(anF, mzabs = 0.01)
      anIC <- groupCorr(anI, cor_eic_th = 0.75)
      anFA <- findAdducts(anIC, polarity="positive")
      peaklist <- getPeaklist(anFA)
      peaklist <- peaklist[order(peaklist$mz),]
      DaDIAtable <- cbind(DaDIAtable, peaklist$isotopes)
      colnames(DaDIAtable)[ncol(DaDIAtable)] <- "Isotopes"
      DaDIAtable <- cbind(DaDIAtable, peaklist$adduct)
      colnames(DaDIAtable)[ncol(DaDIAtable)] <- "Adduct"
      DaDIAtable <- cbind(DaDIAtable, as.numeric(peaklist$pcgroup))
      colnames(DaDIAtable)[ncol(DaDIAtable)] <- "pcgroup"
    }
    write.csv(DaDIAtable, file = "annotated_output.csv")
    
  } else if(num.samples > 1){
    #Metabolite annotation multi sample
    print("Performing dot product annotation ...")
    featureTable <- cbind(featureTable, 0)
    colnames(featureTable)[ncol(featureTable)] <- "Annotation"
    featureTable <- cbind(featureTable, 0)
    colnames(featureTable)[ncol(featureTable)] <- "DPscore"
    d <- foreach(x = 1:nrow(MS2_Spectra_Table)) %dopar% {
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
      
      outVector <- as.numeric(c(0,0))
      #record dot product
      if(nrow(output > 0)){
        output <- output[order(-output[,2]),]
        dot.product <- output[1,2]
        outVector[1] <- dot.product
      }
      
      # Dp score threshold, Dp score >= 0.7 , match_No >= 6 (used in GNPS identification)
      output <- output[output[,2] >= dot.product.threshold,]
      output <- output[output[,3] >= match.number.threshold,]
      
      if(nrow(output)==0) {feature.identity <- 'unknown'}
      if(nrow(output)> 0) {
        output <- output[order(-output[,2]),] # sort by scores
        feature.identity <- output[1,1] # Rank 1, std name
      }
      outVector[2] <- feature.identity
      return(outVector)
    }
    for(x in 1:length(d)){
      MS2_Spectra_Table[x,7] <- d[[x]][2]
      featureTable$Annotation[MS2_Spectra_Table$ID[x]] <- d[[x]][2]
      MS2_Spectra_Table[x,8] <- d[[x]][1]
      featureTable$DPscore[MS2_Spectra_Table$ID[x]] <- d[[x]][1]
    }
    print("Dot product annotation finished")
    print(Sys.time() - start_time)
    
    if(adduct_isotope.annotation){
      print("Performing CAMERA adduct & isotope annotation ...")
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
      print("Finished CAMERA annotation")
      print(Sys.time() - start_time)
    }
    write.csv(featureTable, file = "annotated_output.csv")
  }
}

if(export.mgf){
  print("Exporting individual mgf files ...")
  dir.create("SWATHmgf")
  setwd("SWATHmgf")
  for(y in 1:length(combined_Spectra)){
    if(is.null(combined_Spectra[[y]]) == FALSE){
      writeMgfData(combined_Spectra[[y]], con = paste0(featureTable$mz[y], "_", 
                                                       featureTable$RT[y], "_", "SWATH.mgf"))
    }
  }
  
  setwd(DIA.directory)
  dir.create("DDAmgf")
  setwd("DDAmgf")
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
        writeMgfData(finalSpectra, con = paste0(featureTable$mz[i], "_", 
                                                featureTable$RT[i], "_", "DDA.mgf"))
      }
    }
  }
  print("Finished exporting individual mgf files")
  print(Sys.time() - start_time)
}

if(combine.mgf){
  print("Combining all mgf files ...")
  setwd(DIA.directory)
  setwd("DDAmgf")
  DDAmgfs <- list.files()
  mgf <- read.delim(DDAmgfs[1],stringsAsFactors = FALSE, header = FALSE)
  for(m in 2:length(DDAmgfs)){
    mgf[nrow(mgf)+1,] <- NA
    mgftmp <- read.delim(DDAmgfs[m],stringsAsFactors = FALSE, header = TRUE)
    colnames(mgftmp) <- colnames(mgf) 
    mgf <- rbind(mgf, mgftmp)
  }
  setwd(DIA.directory)
  setwd("SWATHmgf")
  SWATHmgfs <- list.files()
  for(m in 1:length(SWATHmgfs)){
    mgf[nrow(mgf)+1,] <- NA
    mgftmp <- read.delim(SWATHmgfs[m],stringsAsFactors = FALSE, header = TRUE)
    colnames(mgftmp) <- colnames(mgf) 
    mgf <- rbind(mgf, mgftmp)
  }
  setwd(DIA.directory)
  write.table(mgf, file = "combined_mgf.mgf", row.names = FALSE, col.names = FALSE, 
              quote = FALSE, na = "")
  print("Finished combining mgf files")
  print(Sys.time() - start_time)
}

end_time <- Sys.time()
print(end_time - start_time)
#Output identification results  
print("DaDIA workflow completed")


#clean up the cluster
stopImplicitCluster()
