library(tidyverse)
library(vegan)
library(viridis)
library(ggvegan)
library(ggpubr)

APwATs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/heatmapASV/master/Data/APwATs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
oASVs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/heatmapASV/master/Data/oASVs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
xTXs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/heatmapSPP/master/Data/xTXs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
rownames(oASVs) <- oASVs[, 1]
oASVs <- oASVs[, -1]
rownames(xTXs) <- xTXs[, 2]
xTXs <- xTXs[, -1]

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

Shared <- APwATs %>% filter(Cha == 1 & Fla == 1 & Hu == 1 & Pc == 1)
ShaSeqs <- APwATs[, 2]
