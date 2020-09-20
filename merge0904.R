library(readxl)
setwd('C:/Users/hjwang/Desktop')
data <- read.csv("systeMHC_context.csv",header = T)
data <- data[,-1]
hla <- read.csv("syshla.csv",header = F)
data[,8] <- "systeMHC"
colnames(data) <- c("geneSymbol","peptidSequence","ensemblID","proteinSequence","upper","lower","exp","Source")

setwd('C:/Users/hjwang/Desktop/新建文件夹 (2)')
for (i in hla[,1]) {
  filename1 <- paste0(i,".xlsx")
  path <- filename1
  data_stable <- lapply(excel_sheets(path)[1],read_excel,path = path)
  colnames(data_stable) <- c("delect","geneSymbol","peptidSequence","ensemblID","proteinSequence","upper","lower","exp","Source")
  raw_data <- matrix(NA,nrow = length(data_stable[[1]]$peptidSequence),ncol = 8)
  colnames(raw_data)<-c("geneSymbol","peptidSequence","ensemblID","proteinSequence","upper","lower","exp","Source")
  raw_data[,1] <- data_stable[[1]]$geneSymbol
  raw_data[,2] <- data_stable[[1]]$peptidSequence
  raw_data[,3] <- data_stable[[1]]$ensemblID
  raw_data[,4] <- data_stable[[1]]$proteinSequence
  raw_data[,5] <- data_stable[[1]]$upper
  raw_data[,6] <- data_stable[[1]]$lower
  raw_data[,7] <- data_stable[[1]]$exp
  raw_data[,8] <- data_stable[[1]]$Source
  
  sit = which(data$exp == i)
  data2 <- data[sit,1:8]
  data2[,7] <- NA
  
  data3 <- rbind(raw_data,data2)
  filename2 <- paste0(i,".csv")
  write.csv(data3,file = filename2)
  
  data4 <- data3[!duplicated(data3[,c(2,4)]),]
  filename2 <- paste0(i,"_duplicated.csv")
  write.csv(data4,file = filename2)
}

library(readxl)
setwd('C:/Users/hjwang/Desktop')
hla <- read.csv("aq.csv",header = F)
hla[1,1] <- "A0202"
hla <- hla[-31,1]
hla <- hla[-42]
hla <- hla[-42]
setwd('C:/Users/hjwang/Desktop/新建文件夹 (2)')
for (i in hla) {
  filename1 <- paste0(i,".xlsx")
  path <- filename1
  data_stable <- lapply(excel_sheets(path)[1],read_excel,path = path)
  colnames(data_stable[[1]]) <- c("delect","geneSymbol","peptidSequence","ensemblID","proteinSequence","upper","lower","exp","Source")
  raw_data <- matrix(NA,nrow = length(data_stable[[1]]$peptidSequence),ncol = 8)
  colnames(raw_data)<-c("geneSymbol","peptidSequence","ensemblID","proteinSequence","upper","lower","exp","Source")
  raw_data[,1] <- data_stable[[1]]$geneSymbol
  raw_data[,2] <- data_stable[[1]]$peptidSequence
  raw_data[,3] <- data_stable[[1]]$ensemblID
  raw_data[,4] <- data_stable[[1]]$proteinSequence
  raw_data[,5] <- data_stable[[1]]$upper
  raw_data[,6] <- data_stable[[1]]$lower
  raw_data[,7] <- data_stable[[1]]$exp
  raw_data[,8] <- data_stable[[1]]$Source
  
  
  filename2 <- paste0(i,".csv")
  write.csv(raw_data,file = filename2)
  
  data4 <- raw_data[!duplicated(raw_data[,c(2,4)]),]
  filename2 <- paste0(i,"_duplicated.csv")
  write.csv(data4,file = filename2)
}
