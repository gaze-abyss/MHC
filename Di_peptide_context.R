library(readxl)
library(stringr)

setwd('C:/Users/hjwang/Desktop/DI Marco_raw')
stable_data <- read.csv("stable_data.csv",header = T)
sheet_Pearson <- excel_sheets("JI_1700938_Supplemental_Table_1.xlsx")
path <- "JI_1700938_Supplemental_Table_1.xlsx"
data <- lapply(excel_sheets(path)[1],read_excel,path = path)

raw_data <- matrix(NA,nrow = length(data[[1]]$Sequence),ncol = 7)
colnames(raw_data) <- c("UniprotID","peptidSequence","ensemblID","proteinSequence","upper","lower","HLA")
raw_data[,1] <- data[[1]]$`Uniprot ID`
raw_data[,1] <- toupper(raw_data[,1])
raw_data[,2] <- data[[1]]$Sequence
raw_data[,2] <- toupper(raw_data[,2])
raw_data[,7] <- data[[1]]$HLA

for (i in 1:length(raw_data[,1])) {
  if (is.na(raw_data[i,1])) {
    raw_data[i,3:6] <- NA
  }else{
    sit <- which(stable_data[,12] == raw_data[i,1])
    sit <- sit[which.max(nchar(stable_data[sit,11]))]
    if (length(which(stable_data[,12] == raw_data[i,1]))==0) {
      raw_data[i,3:6] <- NA
    }else{
      raw_data[i,3] <- stable_data[sit,2]
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
write.csv(raw_data,"DiMarco_context.csv")

raw_data <- read.csv("DiMarco_context.csv",header = T)
raw_data <- raw_data[,-1]
HLA <- unique(raw_data[,7])
for (i in HLA) {
  sit <- which(raw_data[,7] == i)
  classify_rawdata <- matrix(NA,length(sit),7)
  classify_rawdata <- raw_data[sit,]
  filename <- paste0(i,".csv")
  write.csv(classify_rawdata,filename)
}
setwd('C:/Users/hjwang/Desktop/DI Marco_raw')
for (i in HLA) {
  filename <- paste0(i,".csv")
  classify_rawdata <- read.csv(filename,header = T)
  classify_rawdata <- classify_rawdata[,-1]
  na_sit <- which(is.na(classify_rawdata[,4])==TRUE)
  if(length(na_sit) == 0){
    classify_rawdata <- classify_rawdata
  }else{
    classify_rawdata <- classify_rawdata[-na_sit,]
  }
  classify_rawdata <- classify_rawdata[!duplicated(classify_rawdata),]
  notfind <- Reduce('paste0',rep('X',32))
  not_sit <- which(classify_rawdata[,5] == notfind)
  if(length(not_sit)==0){
    classify_rawdata <- classify_rawdata
  }else{
    classify_rawdata <- classify_rawdata[-not_sit,]
  }
  filename <- paste0(i,"clean",".csv")
  write.csv(classify_rawdata,filename)
}

