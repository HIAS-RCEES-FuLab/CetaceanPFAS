#selection of SID
a1 <- readLines("Sample.txt")
scan <- grep('</scan>',a1)

#full scan
a1_sid <- a1
lowmz <- grep('lowMz="49',a1)
for (i in 1:length(lowmz)) {
  c <- lowmz[i]
  a1_sid <- a1_sid[-((c-7):(c+10))]
  lowmz <- lowmz-18
}
scan1 <- grep('</scan>',a1_sid)
writeLines(a1_sid,"Sample_full.txt")

#sid=80
a1_sid <- a1
highmz <- grep('highMz="1515',a1)
for (i in 1:length(highmz)) {
  c <- highmz[i]
  a1_sid <- a1_sid[-((c-8):(c+9))]
  highmz <- highmz-18
}
scan2 <- grep('</scan>',a1_sid)
writeLines(a1_sid,"Sample_SID.txt")

#readMSdata
library(readMzXmlData)
MS1 <- readMzXmlFile("Sample_full.mzXML")
MS2 <- readMzXmlFile("Sample_SID.mzXML")

RT <- as.numeric()
for (i in 1:length(MS2)) {
  rt <- MS2[[i]][["metaData"]][["retentionTime"]]
  RT <- c(RT,rt)
}

#C2F5
#select scan number
scan <- which(RT>470 & RT<483)
Inten <- as.numeric()  #extracted ion chromatogram
for (i in min(scan):max(scan)) {
  a <- abs(MS2[[i]][["spectrum"]][["mass"]]-118.993)/118.993*1000000<5
  peakint <- sum(MS2[[i]][["spectrum"]][["intensity"]][a])
  Inten <- c(Inten,peakint)
}
plot(RT[scan], Inten, type="h")
frag_scan <- scan[which.max(Inten)]
frag_RT <- MS2[[frag_scan]][["metaData"]][["retentionTime"]]

median(scan)
mz <- MS1[[median(scan)]]$spectrum$mass
intensity <- MS1[[median(scan)]]$spectrum$intensity

#region-specific deconvolution
Parameter <- as.numeric()
for (j in 1:length(mz)) {
  selmz <- mz[intensity==max(intensity)]  #select mz with maximum intensity
  
  PeakInt <- as.numeric()  #extracted ion chromatography
  for (i in min(scan):max(scan)) {
    a <- abs(MS1[[i]][["spectrum"]][["mass"]]-selmz)/selmz*1000000<5
    peakint <- sum(MS1[[i]][["spectrum"]][["intensity"]][a])
    PeakInt <- c(PeakInt,peakint)
  }
  #plot(RT[scan], PeakInt, type="h")
  
  cor_parameter <- cor(PeakInt,Inten)  #correlation coefficient
  cos_sim <-  sum(PeakInt*Inten)/(sqrt(sum(PeakInt^2))*sqrt(sum(Inten^2)))  #cosine simility
  
  pre_scan <- scan[which.max(PeakInt)]
  pre_RT <- MS1[[pre_scan]][["metaData"]][["retentionTime"]]
  delta_RT <- abs(pre_RT-frag_RT)
  
  para <- c(selmz,cor_parameter,cos_sim,delta_RT)
  Parameter <- rbind(Parameter,para)  #combine parameters
  
  mz_1 <- which(abs((mz-selmz)/selmz) < 0.00005)  #next j
  intensity[mz_1] <- 0
  if (max(intensity)<max(Inten)*0.5) break
}
colnames(Parameter) <- c("mz","Pearson","cosine","delta_RT")

write.csv(Parameter,"C2F5.csv")


#C3F7
#select scan number
scan <- which(RT>602.5 & RT<606.7)
Inten <- as.numeric()  #extracted ion chromatogram
for (i in min(scan):max(scan)) {
  a <- abs(MS2[[i]][["spectrum"]][["mass"]]-168.9894)/168.9894*1000000<5
  peakint <- sum(MS2[[i]][["spectrum"]][["intensity"]][a])
  Inten <- c(Inten,peakint)
}
plot(RT[scan], Inten, type="h")
frag_scan <- scan[which.max(Inten)]
frag_RT <- MS2[[frag_scan]][["metaData"]][["retentionTime"]]

median(scan)
mz <- MS1[[median(scan)]]$spectrum$mass
intensity <- MS1[[median(scan)]]$spectrum$intensity

#region-specific deconvolution
Parameter <- as.numeric()
for (j in 1:length(mz)) {
  selmz <- mz[intensity==max(intensity)]  #select mz with maximum intensity
  
  PeakInt <- as.numeric()  #extracted ion chromatography
  for (i in min(scan):max(scan)) {
    a <- abs(MS1[[i]][["spectrum"]][["mass"]]-selmz)/selmz*1000000<5
    peakint <- sum(MS1[[i]][["spectrum"]][["intensity"]][a])
    PeakInt <- c(PeakInt,peakint)
  }
  #plot(RT[scan], PeakInt, type="h")
  
  cor_parameter <- cor(PeakInt,Inten)  #correlation coefficient
  cos_sim <-  sum(PeakInt*Inten)/(sqrt(sum(PeakInt^2))*sqrt(sum(Inten^2)))  #cosine simility
  
  pre_scan <- scan[which.max(PeakInt)]
  pre_RT <- MS1[[pre_scan]][["metaData"]][["retentionTime"]]
  delta_RT <- abs(pre_RT-frag_RT)
  
  para <- c(selmz,cor_parameter,cos_sim,delta_RT)
  Parameter <- rbind(Parameter,para)  #combine parameters
  
  mz_1 <- which(abs((mz-selmz)/selmz) < 0.00005)  #next j
  intensity[mz_1] <- 0
  if (max(intensity)<max(Inten)*0.5) break
}
colnames(Parameter) <- c("mz","Pearson","cosine","delta_RT")

write.csv(Parameter,"C3F7.csv")


