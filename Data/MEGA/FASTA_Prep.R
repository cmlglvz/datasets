library(tidyverse)
library(ShortRead)
library(DECIPHER)
library(Biostrings)

Tax.sum <- function(OTU.Table, Tax.Table, Tax.lvl ){
  z <- NULL
  y <- NULL
  for (i in 1:length(unique(Tax.Table[colnames(OTU.Table),Tax.lvl]))) {
    if (length(OTU.Table[,which(Tax.Table[colnames(OTU.Table),Tax.lvl]==unique(Tax.Table[colnames(OTU.Table),Tax.lvl])[i])])!=length(rownames(OTU.Table))) {
      z <- which(Tax.Table[colnames(OTU.Table),Tax.lvl]==unique(Tax.Table[colnames(OTU.Table),Tax.lvl])[i])
      y <- cbind(y, apply(OTU.Table[,which(Tax.Table[colnames(OTU.Table),Tax.lvl]==unique(Tax.Table[colnames(OTU.Table),Tax.lvl])[i])], 1, function(x) sum(x)))
    } else { 
      y <- cbind(y, OTU.Table[,which(Tax.Table[colnames(OTU.Table),Tax.lvl]==unique(Tax.Table[colnames(OTU.Table),Tax.lvl])[i])])
    }
  }
  colnames(y) <- unique(Tax.Table[colnames(OTU.Table),Tax.lvl])
  invisible((y))
}

rltv.Otu.Table <- function(x){
  x.Data.rltv <- NULL
  for (i in 1:dim(x)[1]) {
    x.Data.rltv <- rbind(x.Data.rltv, x[i,]/apply(x, 1, function(x) sum(x))[i])
  }
  rownames(x.Data.rltv) <- rownames(x)
  invisible(x.Data.rltv)
}

APwATs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/APwATs.csv", sep = ";", dec = ".", skip = 0)
wASVs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/wASVs.csv", sep = ";", dec = ".", skip = 0)
rownames(wASVs) <- c("C1_2017.08", "C1_2018.02", "C1_2018.08", "C2_2017.08", "C2_2018.02", "C2_2018.08", "C3_2017.08", "C3_2018.02", "C3_2018.08", "C4_2017.08", "C4_2018.02", "C4_2018.08", "F1_2017.08", "F1_2018.02", "F1_2018.08", "F2_2017.08", "F2_2018.02", "F2_2018.08", "F3_2017.08", "F3_2018.02", "F3_2018.08", "F4_2017.08", "F4_2018.02", "F4_2018.08", "H1_2017.08", "H1_2018.02", "H1_2018.08", "H2_2017.08", "H2_2018.02", "H2_2018.08", "H3_2017.08", "H3_2018.02", "H3_2018.08", "H4_2017.08", "H4_2018.02", "H4_2018.08", "P1_2018.02", "P1_2018.08", "P2_2018.02", "P3_2018.02", "P4_2018.02")
wASVs <- wASVs[, -1]
wTaxa <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/xTXs.csv", sep = ";", dec = ".", skip = 0)
rownames(wTaxa) <- wTaxa[, 2]
wTaxa <- wTaxa[, -1]

uSha <- APwATs %>% filter(Cha == 1 & Fla == 1 & Hu == 1 & Pc == 1)
shaseq <- uSha[, 2]
ShaASV <- wASVs %>% select(all_of(shaseq)) #equivalente a seqtab.nochim
write.csv2(ShaASV, file = "C:/Users/Camilo/Dropbox/R/MEGA/seqtab_shared.csv")
ShaTaxa <- wTaxa %>% filter(Seq %in% all_of(shaseq))
write.csv2(ShaTaxa, file = "C:/Users/Camilo/Dropbox/R/MEGA/shared_taxa.csv")
Shared <- as.data.frame(t(ShaASV))
Shared <- add_column(Shared, Seq = all_of(rownames(Shared)), .before = "C1_2017.08")
Shared <- inner_join(ShaTaxa, Shared, by = "Seq")
ASV <- rownames(ShaTaxa)
Shared <- add_column(Shared, ASV = ASV, .before = "Seq")
rownames(Shared) <- ASV
write.csv2(Shared, file = "C:/Users/Camilo/Dropbox/R/MEGA/shared_complete_table.csv")
