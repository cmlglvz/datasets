library(tidyverse)
library(dada2)
library(ShortRead)
library(DECIPHER)
library(Biostrings)

wASVs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/wASVs.csv", sep = ";", dec = ".", skip = 0, fill = TRUE)
xTXs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/xTXs.csv", sep = ";", dec = ".", skip = 0)
rownames(wASVs) <- c("C1_2017.08", "C1_2018.02", "C1_2018.08", "C2_2017.08", "C2_2018.02", "C2_2018.08", "C3_2017.08", "C3_2018.02", "C3_2018.08", "C4_2017.08", "C4_2018.02", "C4_2018.08", "F1_2017.08", "F1_2018.02", "F1_2018.08", "F2_2017.08", "F2_2018.02", "F2_2018.08", "F3_2017.08", "F3_2018.02", "F3_2018.08", "F4_2017.08", "F4_2018.02", "F4_2018.08", "H1_2017.08", "H1_2018.02", "H1_2018.08", "H2_2017.08", "H2_2018.02", "H2_2018.08", "H3_2017.08", "H3_2018.02", "H3_2018.08", "H4_2017.08", "H4_2018.02", "H4_2018.08", "P1_2018.02", "P1_2018.08", "P2_2018.02", "P3_2018.02", "P4_2018.02")
wASVs <- wASVs[, -1]
write.csv2(wASVs, file = "D:/Documents/GitHub/datasets/datasets/Data/eAnalisis/wASVs.csv")
rownames(xTXs) <- xTXs[, 1]
xTXs <- xTXs[, -1]
seqs <- colnames(wASVs)
wTXs <- filter(xTXs, Seq %in% all_of(seqs))
uniques <- wTXs[, -1]
uniques$OTU <- paste(">", uniques$OTU, sep = "")
rownames(uniques) <- c()
write.csv2(uniques, file = "C:/Users/Camilo/Dropbox/R/Analisis/fasta_prep_ppe.csv")
#edit with notepad++