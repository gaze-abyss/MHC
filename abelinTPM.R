setwd('C:/Users/hjwang/Desktop/rawdataxtalpi/HLAthena_raw')
ucsc2stable <- read.csv("ucsc2stable ID.csv",header = T)
abelinTPM <- read.csv("C:/Users/hjwang/Desktop/AbelinTPM.csv",header = T)

ucsc2stable$TPM <- NA
for (i in 1:length(abelinTPM[,1])) {
  if (length(which(abelinTPM[i,1] == ucsc2stable$UCSC.Stable.ID)) != 0) {
    sit = which(abelinTPM[i,1] == ucsc2stable$UCSC.Stable.ID)
    ucsc2stable$TPM[sit] = abelinTPM$mean[i]
  }else{
    next
  }
}

write.csv(ucsc2stable,file = "TPMstable.csv")

abelinHLA <- c("A0101","A0201","A0203","A0204","A0207","A0301","A2402","A2902","A3101","A6802","B3501","B4402","B4403","B5101","B5401","B5701")
TPMdata <- read.csv("C:/Users/hjwang/Desktop/TPMstable.csv",header = T)
setwd('C:/Users/hjwang/Desktop/HLAthena_clean')
for (i in abelinHLA) {
  singledata = paste0(i,"clean.csv")
  data = read.csv(singledata,header = T,row.names = 1)
  for (o in 1:length(data[,1])) {
    if (length(which(data$geneSymbol[o] == abelinTPM$锘縯ranscript_id)) != 0) {
      sit = which(data$geneSymbol[o] == abelinTPM$锘縯ranscript_id)
      data$exp[o] = abelinTPM$mean[sit]
    }else{
      next
    }
  }
  filename <- paste0(i,"_TPM.csv")
  write.csv(data,file = filename)
}
