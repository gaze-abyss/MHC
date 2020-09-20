setwd('C:/Users/hjwang/Desktop')
mulu <- read.csv("mulu.csv",header = FALSE)
mulu[1,] <- "A0101"

hlethena <- mulu[1:95,]
data2 <- matrix(NA,nrow = 139,ncol = 18)
data2[,1] <- mulu[,1]

for (i in hlethena) {
  filename <- paste0("C:/Users/hjwang/Desktop/HLAthena_clean/",i,"clean.csv")
  data <- read.csv(filename,header = T,row.names = 1)
  data$length <- nchar(data$peptidSequence)
  sit = which(hlethena == i)
  for (o in 7:23) {
    data2[sit,o-5] <- length(which(data$length == o))
  }
}


pearson <- mulu[96:122,]
for (i in pearson) {
  filename <- paste0("C:/Users/hjwang/Desktop/pearson_clean/",i,"clean.csv")
  data <- read.csv(filename,header = T,row.names = 1)
  data$length <- nchar(data$peptidSequence)
  sit = which(pearson == i) + 95
  for (o in 7:23) {
    data2[sit,o-5] <- length(which(data$length == o))
  }
}

dimarco <- mulu[123:139,]
for (i in dimarco) {
  filename <- paste0("C:/Users/hjwang/Desktop/DiMarco_clean/",i,"clean.csv")
  data <- read.csv(filename,header = T,row.names = 1)
  
  data$length <- nchar(data$peptidSequence)
  sit = which(dimarco == i) + 122
  for (o in 7:23) {
    data2[sit,o-5] <- length(which(data$length == o))
  }
}

write.csv(data2,"cleanstat.csv")

setwd('C:/Users/hjwang/Desktop')
mulu <- read.csv("mulu.csv",header = FALSE)
mulu[1,] <- "A0101"

abelin <- mulu[1:16,]
data3 <- matrix(NA,nrow = 16,ncol = 18)
data3[,1] <- abelin

for (i in abelin) {
  filename <- paste0("C:/Users/hjwang/Desktop/HLAthena_clean/abelin_tpm/",i,"_TPM.csv")
  data <- read.csv(filename,header = T,row.names = 1)
  data$length <- nchar(data$peptidSequence)
  sit = which(data3[,1] == i)
  for (o in 7:23) {
    data3[sit,o-5] <- length(which(data$length == o))
  }
}
write.csv(data3,"abelincleanstat.csv")
