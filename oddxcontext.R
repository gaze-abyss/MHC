setwd('C:/Users/hjwang/Desktop/rawdataxtalpi/HLAthena_raw')
stable_data <- read.csv("C:/Users/hjwang/Desktop/stable_data.csv",header = T)

library(readxl)
library(stringr)

sheet_HLAthena <- excel_sheets("41587_2019_322_MOESM3_ESM.xlsx")
path <- "41587_2019_322_MOESM3_ESM.xlsx"

for (o in 2:length(sheet_HLAthena)) {
  
  data <- lapply(excel_sheets(path)[o],read_excel,path = path)
  raw_data <- matrix(NA,nrow = length(data[[1]]$entry_name),ncol = 7)
  colnames(raw_data) <- c("geneSymbol","peptidSequence","ensemblID","proteinSequence","upper","lower","exp")
  raw_data[,1] <- data[[1]]$entry_name
  raw_data[,1] <- substring(raw_data[,1],1,8)
  raw_data[,2] <- data[[1]]$sequence
  raw_data[,2] <- toupper(raw_data[,2])
  
  for (i in 1:length(raw_data[,1])) {
    if (is.na(raw_data[i,1])) {
      raw_data[i,3:6] <- NA
    }else{
      sit <- which(stable_data[,12] == raw_data[i,1])
      sit <- sit[which.max(nchar(stable_data[sit,11]))]
      if (length(which(stable_data[,12] == raw_data[i,1]))==0) {
        raw_data[i,3:6] <- NA
      }else{
        raw_data[i,3] <- stable_data$transcript_stable_id[sit]
        raw_data[i,4] <- stable_data[sit,11]
        protein_len <- nchar(stable_data[sit,11])
        peptide_start <- regexec(raw_data[i,2],stable_data[sit,11])[[1]][1]
        peptide_len <- nchar(raw_data[i,2])
        if (peptide_start > 30 & (protein_len - (peptide_len + peptide_start)) > 29) {
          raw_data[i,5] <- substring(stable_data[sit,11],peptide_start - 30,peptide_start - 1)
          raw_data[i,6] <- substring(stable_data[sit,11],peptide_start + peptide_len, peptide_start + peptide_len + 29)
        }else{
          if (peptide_start < 30 & (protein_len - (peptide_len + peptide_start)) < 29) {
            oddupx <- Reduce('paste0',rep('X',31 - peptide_start))
            oddupamino <- substring(stable_data[sit,11],1,peptide_start - 1)
            raw_data[i,5] <- paste0(oddupx,oddupamino)
            
            oddlowx <- Reduce('paste0',rep('X',30 - (protein_len - (peptide_len + peptide_start))))
            oddlowamino <- substring(stable_data[sit,11],peptide_start + peptide_len, protein_len)
            raw_data[i,6] <- paste0(oddlowamino,oddlowx)
          }else{
            if (peptide_start <= 30 & (protein_len - (peptide_len + peptide_start)) > 29) {
              oddupx <- Reduce('paste0',rep('X',31 - peptide_start))
              oddupamino <- substring(stable_data[sit,11],1,peptide_start - 1)
              raw_data[i,5] <- paste0(oddupx,oddupamino)
              raw_data[i,6] <- substring(stable_data[sit,11],peptide_start + peptide_len, peptide_start + peptide_len + 29)
            }else{
              raw_data[i,5] <- substring(stable_data[sit,11],peptide_start - 30,peptide_start - 1)
              oddlowx <- Reduce('paste0',rep('X',30 - (protein_len - (peptide_len + peptide_start - 1))))
              oddlowamino <- substring(stable_data[sit,11],peptide_start + peptide_len, protein_len)
              raw_data[i,6] <- paste0(oddlowamino,oddlowx)
            }
          }
        }
      }
    }
  }
  na_sit <- which(is.na(raw_data[,4])==TRUE)
  raw_data <- raw_data[-na_sit,]
  raw_data <- raw_data[!duplicated(raw_data),]
  notfind <- Reduce('paste0',rep('X',32))
  if (length(which(raw_data[,5] == notfind))==0) {
    raw_data <- raw_data
  }else{
  not_sit <- which(raw_data[,5] == notfind)
  raw_data <- raw_data[-not_sit,]
  }
  filename <- paste0(excel_sheets(path)[o],"clean",".csv")
  write.csv(raw_data,filename)
}


