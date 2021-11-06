library(tidyverse)
library(hrbrthemes)
library(viridis)
library(plotly)
library(d3heatmap)
library(vegan)
library(heatmaply)
library(htmlwidgets)

dTSDiv <- read.csv2("/Users/Artemis/Dropbox/R/eAnalisis/dTSDiv.csv", header = T, sep = ";", dec = ",", skip = 0)
rownames(dTSDiv) <- dTSDiv[, 1]
dTSDiv <- dTSDiv[, -c(1, 43:48)]
dTSDiv <- t(dTSDiv)
htmap <- t(dTSDiv)
htmap <- decostand(htmap, method = "normalize")
htmap <- as.matrix(t(dTSDiv))
heatmap(htmap, Colv = NA, Rowv = NA, scale="column")
htmp <- heatmaply(normalize(dTSDiv), 
                  Colv = NA, 
                  Rowv = NA, 
                  seriate = "none", 
                  xlab = "PPE Phylum", 
                  ylab = "Samples", 
                  main = "Normalized data", 
                  file = paste0(getwd(), "heatmap.html"))
htmp
saveWidget(htmp, file=paste0(getwd(), "/htmp.html"))
