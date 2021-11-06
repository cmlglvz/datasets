library(tidyverse)
library(vegan)
library(BiodiversityR)
library(viridis)
library(ggvegan)
library(ggpubr)

APwATs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/APwATs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
oASVs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/oASVs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
xTXs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/xTXs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
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
ShaSeqs <- Shared[, 2]
ShaASVs <- select(oASVs, all_of(ShaSeqs))
ShaSPP <- Tax.sum(ShaASVs, xTXs, 10)
ShaSPP <- ShaSPP[, -17]
write.csv2(ShaSPP, file = "/Users/Artemis/Documents/GitHub/datasets/Data/eAnalisis/ShaSPP.csv")

rltv.Otu.Table <- function(x){
  x.Data.rltv <- NULL
  for (i in 1:dim(x)[1]) {
    x.Data.rltv <- rbind(x.Data.rltv, x[i,]/apply(x, 1, function(x) sum(x))[i])
  }
  rownames(x.Data.rltv) <- rownames(x)
  invisible(x.Data.rltv)
}

rShaSPP <- rltv.Otu.Table(ShaSPP)
write.csv2(rShaSPP, file = "/Users/Artemis/Documents/GitHub/datasets/Data/eAnalisis/relative_ShaSPP.csv")

Whittaker <- radfit(rShaSPP)
plot(Whittaker)

rank <- rankabundance(hellSPP)
ranked <- rownames(rank)[1:20]
shllSPP <- as.data.frame(hellSPP, all_of(ranked))
nmds2 <- metaMDS(shllSPP, autotransform = FALSE)

shllSPP <- mutate(shllSPP, 
                  Site = c(rep("Cha", 12), rep("Fla", 12), rep("Hu", 12), rep("Pc", 5)), 
                  Color = c(rep("#E68431", 12), rep("#D8A64F", 12), rep("#183D84", 12), rep("#80E0FB", 5))
                  )

fort <- fortify(nmds2)
p1 <- ggplot() +
  geom_point(data = subset(fort, Score == 'sites'),
             mapping = aes(x = NMDS1, y = NMDS2),
             colour = c(rep("#E68431", 12), rep("#D8A64F", 12), rep("#183D84", 12), rep("#80E0FB", 5)),
             alpha = 0.8, 
             shape = c(rep(15, 12), rep(17, 12), rep(16, 12), rep(18, 5)), 
             size = 2.5) + 
  geom_segment(data = subset(fort, Score == 'species'),
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.015, "npc"),
                             type = "closed"),
               colour = "darkgray",
               size = 0,
               alpha = 0) +
  geom_text(data = subset(fort, Score == 'sites'),
            mapping = aes(label = c(rep("Cha", 12), rep("Fla", 12), rep("Hu", 12), rep("Pc", 5)),
                          x = NMDS1 * 1.05, 
                          y = NMDS2 * 1.15), 
            check_overlap = TRUE, 
            alpha = 0) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = 0.8, colour = "gray") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 0.8, colour = "gray") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

p2 <- ggplot() +
  geom_point(data = subset(fort, Score == 'sites'),
             mapping = aes(x = NMDS1, y = NMDS2),
             colour = "black",
             alpha = 0) +
  geom_segment(data = subset(fort, Score == 'species'),
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.015, "npc"),
                             type = "closed"),
               colour = "#F71B46",
               size = 0.8) +
  geom_text(data = subset(fort, Score == 'species'),
            mapping = aes(label = Label, x = NMDS1 * 1.1, y = NMDS2 * 1.1)) + 
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = 0.8, colour = "gray") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 0.8, colour = "gray") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

#create a multi-panel plot with one column
ggarrange(p1, p2, ncol = 1)





