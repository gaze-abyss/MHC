library(readxl)
library(stringr)
#read data
setwd('C:/Users/hjwang/Desktop')
stable_data <- read.csv("homo_pfams_202003101726.csv",header = T)
raw_data <- read.csv("systeMHC.csv",header = T)
raw_data <- raw_data[,-1]
#find context
for (i in 1:length(raw_data[,1])) {
  if (is.na(raw_data[i,1])) {
    raw_data[i,4:6] <- NA
  }else{
    sit <- which(stable_data[,4] == raw_data[i,1])
    sit <- sit[which.max(nchar(stable_data[sit,10]))]
    if (length(which(stable_data[,4] == raw_data[i,1]))==0) {
      raw_data[i,4:6] <- NA
    }else{
      raw_data[i,4] <- stable_data[sit,10]
      protein_len <- nchar(stable_data[sit,10])
      peptide_start <- regexec(raw_data[i,2],stable_data[sit,10])[[1]][1]
      peptide_len <- nchar(raw_data[i,2])
      if (peptide_start > 30 & (protein_len - (peptide_len + peptide_start)) > 29) {
        raw_data[i,5] <- substring(stable_data[sit,10],peptide_start - 30,peptide_start - 1)
        raw_data[i,6] <- substring(stable_data[sit,10],peptide_start + peptide_len, peptide_start + peptide_len + 29)
      }else{
        if (peptide_start < 30 & (protein_len - (peptide_len + peptide_start)) < 29) {
          oddupx <- Reduce('paste0',rep('X',31 - peptide_start))
          oddupamino <- substring(stable_data[sit,10],1,peptide_start - 1)
          raw_data[i,5] <- paste0(oddupx,oddupamino)
          
          oddlowx <- Reduce('paste0',rep('X',30 - (protein_len - (peptide_len + peptide_start))))
          oddlowamino <- substring(stable_data[sit,10],peptide_start + peptide_len, protein_len)
          raw_data[i,6] <- paste0(oddlowamino,oddlowx)
        }else{
          if (peptide_start <= 30 & (protein_len - (peptide_len + peptide_start)) > 29) {
            oddupx <- Reduce('paste0',rep('X',31 - peptide_start))
            oddupamino <- substring(stable_data[sit,10],1,peptide_start - 1)
            raw_data[i,5] <- paste0(oddupx,oddupamino)
            raw_data[i,6] <- substring(stable_data[sit,10],peptide_start + peptide_len, peptide_start + peptide_len + 29)
          }else{
            raw_data[i,5] <- substring(stable_data[sit,10],peptide_start - 30,peptide_start - 1)
            oddlowx <- Reduce('paste0',rep('X',30 - (protein_len - (peptide_len + peptide_start - 1))))
            oddlowamino <- substring(stable_data[sit,10],peptide_start + peptide_len, protein_len)
            raw_data[i,6] <- paste0(oddlowamino,oddlowx)
          }
        }
      }
    }
  }
}
#clean data
sit <- which(is.na(raw_data$proteinSequence) == T)
raw_data <- raw_data[-sit,]
write.csv(raw_data,"systeMHC_context.csv")
