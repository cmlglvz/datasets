library(tidyverse)
library(reshape2)
library(hrbrthemes)
library(viridis)

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

ShaASVs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/ShaASVs.csv", header = TRUE, sep = ";", dec = ".", row.names = 1, skip = 0, fill = TRUE)
ShaTXs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/ShaTXs.csv", header = TRUE, sep = ";", dec = ".", row.names = 1, skip = 0, fill = TRUE)
ShaRAs <- as.data.frame(t(ShaASVs))
ShaRAs <- mutate(ShaRAs, 
                 Cha = rowSums(ShaRAs[1:12]), 
                 Fla = rowSums(ShaRAs[13:24]), 
                 Hu = rowSums(ShaRAs[25:36]), 
                 Pc = rowSums(ShaRAs[37:41]), 
                 Total = rowSums(ShaRAs[1:41])
                 )
sRAs <- ShaRAs[, c(42:45)]
sRAs <- rltv.Otu.Table(sRAs)
sRAs <- apply(sRAs, 2, function(x) x * 100) %>% as.data.frame()
sRAs <- as.data.frame(t(sRAs))

ChlTXs <- filter(ShaTXs, Division == "Chlorophyta")
ChlRAs <- select(ShaRAs, all_of(ChlTXs$Seq))
colnames(ChlRAs) <- ChlTXs$OTU
trying <- as.data.frame(t(ChlRAs))
trying <- mutate(trying, 
                 Cha = rowSums(trying[1:12]), 
                 Fla = rowSums(trying[13:24]), 
                 Hu = rowSums(trying[25:36]), 
                 Pc = rowSums(trying[37:41]), 
                 Total = rowSums(trying[1:41])
                 )
aChl <- mutate(ChlRAs, Muestra = rownames(ChlRAs), 
               Site = as.factor(c(rep("Cha", 12), rep("Fla", 12), rep("Hu", 12), rep("Pc", 5))), 
               .before = "ASV_1")

OchTXs <- filter(ShaTXs, Division == "Ochrophyta")
OchASVs <- select(ShaASVs, all_of(OchTXs$Seq))
colnames(OchASVs) <- OchTXs$OTU
OchASVs <- mutate(OchASVs, Site = as.factor(c(rep("Cha", 12), rep("Fla", 12), rep("Hu", 12), rep("Pc", 5))), .before = "ASV_18")

HapTXs <- filter(ShaTXs, Division == "Haptophyta")
HapASVs <- select(ShaASVs, all_of(HapTXs$Seq))
colnames(HapASVs) <- HapTXs$OTU
HapASVs <- mutate(HapASVs, Site = as.factor(c(rep("Cha", 12), rep("Fla", 12), rep("Hu", 12), rep("Pc", 5))), .before = "ASV_38")

CKTXs <- filter(ShaTXs, Division == "Cryptophyta" | Division == "Katablepharidophyta" | Division == "Cryptophyta:nucl")
CKASVs <- select(ShaASVs, all_of(CKTXs$Seq))
colnames(CKASVs) <- CKTXs$OTU
CKASVs <- mutate(CKASVs, Site = as.factor(c(rep("Cha", 12), rep("Fla", 12), rep("Hu", 12), rep("Pc", 5))), .before = "ASV_31")

bChl <- gather(aChl, key = "OTU", value = "Rel.Abun", -c(1, 2))
ggChl <- ggplot(bChl, aes(x = variable, y = value, fill = Site)) + 
  geom_boxplot() + 
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) + 
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) + 
  theme_ipsum() + 
  theme(
    legend.position = "none", 
    plot.title = element_text(size = 11)
    ) + 
  ggtitle("Chlorophyta relative abundance") +
  xlab("")
ggChl

gg <- ggplot(bChl, aes(x = variable, y = value, fill = Site)) + 
  geom_boxplot(colour = "black", position = position_dodge(0.5)) + 
  geom_vline(xintercept = c(1.5,2.5,3.5), colour = "grey85", size = 1.2) +
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), legend.position = "right", 
        axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 12, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.key=element_blank()) + 
  labs(x= "", y = "Relative Abundance (%)", fill = "Site") + 
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, option = "D")
gg

vln <- ggplot(bChl, aes(x = variable, y = value, fill = Site)) +
  geom_violin() + 
  scale_fill_viridis(discrete = TRUE, alpha = 1, option = "A") +
  theme_ipsum() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 11)
    ) + 
  ggtitle("Violin chart") + 
  xlab("")
vln

















