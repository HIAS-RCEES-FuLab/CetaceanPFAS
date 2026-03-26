#extract full scan data
a1 <- readLines("Sample.txt")
scan <- grep('</scan>',a1)
a1_sid <- a1
lowmz <- grep('lowMz="49',a1)
for (i in 1:length(lowmz)) {
  c <- lowmz[i]
  a1_sid <- a1_sid[-((c-7):(c+10))]
  lowmz <- lowmz-18
}
scan1 <- grep('</scan>',a1_sid)
writeLines(a1_sid,"Sample_full.txt")

#peak picking
library(xcms)
file <- "Sample_full.mzXML"
data <- readMSData(file, mode = "onDisk")
table(msLevel(data))

cwp <- CentWaveParam(snthresh = 10, noise = 500, ppm = 5, peakwidth = c(5, 20))
dda <- findChromPeaks(data, param = cwp)
pk <- chromPeaks(dda)
write.csv(pk,"Sample_full.csv")


#homologue analysis
pk_full <- read.csv("Sample_full.csv")
pk_full <- subset(pk_full,AKMD>0.8 & AKMD <1.1)
pk_full <- subset(pk_full,mz>200)

pk <- as.numeric()
i=1
while (i<length(pk_full$mz)) {
  delta_m <- abs((pk_full$mz-pk_full$mz[i])/pk_full$mz[i])
  a <- which(delta_m <0.000005)
  pk_2 <- pk_full[i,]
  pk <- rbind(pk,pk_2)
  i=i+length(a)
}  #delete duplicated mz

CF2 <- 49.9968
nCF2 <- CF2*1:20
PK_homo <- c(1:20)
for (i in 1:length(pk$mz)) {
  delta_m <- abs(pk$mz-pk$mz[i])
  a <- which(abs(delta_m-nCF2[1])<0.001 |
             abs(delta_m-nCF2[2])<0.001 |
             abs(delta_m-nCF2[3])<0.001 |
             abs(delta_m-nCF2[4])<0.001 |
             abs(delta_m-nCF2[5])<0.001 |
             abs(delta_m-nCF2[6])<0.001 |
             abs(delta_m-nCF2[7])<0.001 |
             abs(delta_m-nCF2[8])<0.001 |
             abs(delta_m-nCF2[9])<0.001 |
             abs(delta_m-nCF2[10])<0.001 |
             abs(delta_m-nCF2[11])<0.001 |
             abs(delta_m-nCF2[12])<0.001 |
             abs(delta_m-nCF2[13])<0.001 |
             abs(delta_m-nCF2[14])<0.001 |
             abs(delta_m-nCF2[15])<0.001 |
             abs(delta_m-nCF2[16])<0.001 |
             abs(delta_m-nCF2[17])<0.001 |
             abs(delta_m-nCF2[18])<0.001 |
             abs(delta_m-nCF2[19])<0.001 |
             abs(delta_m-nCF2[20])<0.001)
    if (length(a)>4) {
    pk_homo <- c(pk$mz[i],pk$mz[a])
    PK_homo <- rbind(PK_homo,pk_homo)
    pk <- pk[-a,]
  } else {
    pk <- pk
  }
  
}
write.csv(PK_homo,"homologues.csv")
