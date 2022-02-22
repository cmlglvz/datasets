library(tidyverse)

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

aus.pres <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/APwATs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
ppe.abun <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/wASVs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
ID <- c("C1A17", "C1F18", "C1A18", "C2A17", "C2F18", "C2A18", "C3A17", "C3F18", "C3A18", "C4A17", "C4F18", "C4A18", "F1A17", "F1F18", "F1A18", "F2A17", "F2F18", "F2A18", "F3A17", "F3F18", "F3A18", "F4A17", "F4F18", "F4A18", "H1A17", "H1F18", "H1A18", "H2A17", "H2F18", "H2A18", "H3A17", "H3F18", "H3A18", "H4A17", "H4F18", "H4A18", "P1F18", "P1A18", "P2F18", "P3F18", "P4F18")
rownames(ppe.abun) <- ID
ppe.abun <- ppe.abun[, -1] #Abundancia de ASVs rarefaccionadas filtradas por ASV de PPEs de interes
#eASVs <- wASVs[-c(1:12, 25:41), ] #Solo contamos las muestras del Sitio de interes (no es necesario)
colnames(aus.pres)[3:43] <- ID
taxa <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/xTXs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
rownames(taxa) <- taxa[, 2]
taxa <- taxa[, -1] #La identificación se realiza contra todas las asignaciones, las ASV de los sitios mandan

div.tot <- Tax.sum(ppe.abun, taxa, 5) %>% as.data.frame()
colnames(div.tot)[c(3,6)] <- c("A", "B")
div.tot <- div.tot %>% 
  mutate(Cryptophyta = rowSums(div.tot[c(3, 6)])) %>% 
  select(Chlorophyta, Ochrophyta, Cryptophyta, Haptophyta, Katablepharidophyta, Rhodophyta)

#Todas las ASVs de PPE presentes en todas las muestras, pero no exclusivas para un sitio específico
shared <- aus.pres %>% filter(Cha == 1 & Fla == 1 & Hu == 1 & Pc == 1) #ASV presentes en todos los sitios (al menos en una de las muestras correspondientes) "al mismo tiempo"
#print(all(all_of(shared$Seq)%in%all_of(aus.pres$Seq))) #Nos aseguramos que las secuencias únicas efectivamente se encuentren en todas las secuencias presentes en Huasco
sha.abun <- select(ppe.abun, all_of(shared$Seq)) #Abundancia para las 158 ASV compartidas por los sitios
div.sha <- Tax.sum(sha.abun, taxa, 5) %>% as.data.frame()
colnames(div.sha)[c(3,6)] <- c("A", "B")
div.sha <- div.sha %>% 
  mutate(Cryptophyta = rowSums(div.sha[c(3,6)])) %>% 
  select(Chlorophyta, Ochrophyta, Cryptophyta, Haptophyta, Katablepharidophyta)
write.csv2(div.tot, file = "../datasets/Data/Revisited/ppe_total_division_abundance.csv")
write.csv2(div.sha, file = "../datasets/Data/Revisited/ppe_shared_division_abundance.csv")

sha <- shared$Seq
diff.abun <- select(ppe.abun, -sha)
div.diff <- Tax.sum(diff.abun, taxa, 5) %>% 
  as.data.frame()
colnames(div.diff)[c(3,6)] <- c("A", "B")
div.diff <- div.diff %>% 
  mutate(Cryptophyta = rowSums(div.diff[c(3,6)])) %>% 
  select(Chlorophyta, Ochrophyta, Cryptophyta, Haptophyta, Katablepharidophyta, Rhodophyta)
write.csv2(div.diff, file = "../datasets/Data/Revisited/ppe_not_shared_division_abundance.csv")

cha <- aus.pres %>% filter(Cha == 1)
unique.cha <- cha %>% filter(Fla == 0 & Hu == 0 & Pc == 0)
cha.abun <- select(ppe.abun, all_of(unique.cha$Seq))
write.csv2(cha.abun, file = "../datasets/Data/Revisited/unique_chanaral_asv.csv")
div.cha <- Tax.sum(cha.abun, taxa, 5) %>% as.data.frame()
colnames(div.cha)[4] <- "Cryptophyta"
div.cha <- div.cha[, c(3,1,4,2,5,6)]
write.csv2(div.cha, file = "../datasets/Data/Revisited/ppe_unique_chanaral_division_abundance.csv")

fla <- aus.pres %>% filter(Fla == 1)
unique.fla <- fla %>% filter(Cha == 0 & Hu == 0 & Pc == 0)
fla.abun <- select(ppe.abun, all_of(unique.fla$Seq))
write.csv2(fla.abun, file = "../datasets/Data/Revisited/unique_flamenco_asv.csv")
div.fla <- Tax.sum(fla.abun, taxa, 5) %>% as.data.frame()
div.fla <- div.fla[, c(3,1,6,2,5,4)]
write.csv2(div.fla, file = "../datasets/Data/Revisited/ppe_unique_flamenco_division_abundance.csv")

hu <- aus.pres %>% filter(Hu == 1)
unique.hu <- hu %>% filter(Cha == 0 & Fla == 0 & Pc == 0)
hu.abun <- select(ppe.abun, all_of(unique.hu$Seq))
write.csv2(hu.abun, file = "../datasets/Data/Revisited/unique_huasco_asv.csv")
div.hu <- Tax.sum(hu.abun, taxa, 5) %>% as.data.frame()
div.hu <- div.hu[, c(1,2,4,3)]
write.csv2(div.hu, file = "../datasets/Data/Revisited/ppe_unique_huasco_division_abundance.csv")

pc <- aus.pres %>% filter(Pc == 1)
unique.pc <- pc %>% filter(Cha == 0 & Fla == 0 & Hu == 0)
pc.abun <- select(ppe.abun, all_of(unique.pc$Seq))
write.csv2(pc.abun, file = "../datasets/Data/Revisited/unique_punta_choros_asv.csv")
div.pc <- Tax.sum(pc.abun, taxa, 5) %>% as.data.frame()
div.pc <- div.pc[, c(6,1,2,4,5,3)]
write.csv2(div.pc, file = "../datasets/Data/Revisited/ppe_unique_punta_choros_division_abundance.csv")

#Proximo dataframe fue editado externamente = merge de los 4 df anteriores
div.unique <- read.csv2(file = "Data/Revisited/ppe_unique_division_abundance.csv", header = TRUE, sep = ";", dec = ".", row.names = 1, skip = 0)
div.diff2 <- read.csv2(file = "Data/Revisited/ppe_not_shared_division_abundance_minus_unique.csv", 
                       header = TRUE, 
                       sep = ";",
                       dec = ".", 
                       row.names = 1, 
                       skip = 0)
div.diff2 <- div.diff2[, -c(1:12)]
div.comp <- read.csv2(file = "Data/Revisited/ppe_composite_total_division_abundance.csv", 
                      header = TRUE, 
                      sep = ";", 
                      dec = ".", 
                      row.names = 1, 
                      skip = 0)
div.comp <- div.comp %>% 
  mutate(Unique = rowSums(div.comp[c(7:12)])) %>% 
  mutate(Other = rowSums(div.comp[c(13:18)]))
div.comp <- div.comp[, -c(7:18)]

relative.comp <- rltv.Otu.Table(div.comp)
apply(relative.comp, 1, function(x) sum(x))[1:41]

tiff("Relative_Abundance_Phylum_Composite_PPE_ASV.tiff", width = 10, height = 8, units = 'in', res = 600)
par(mar = c(5.1,4.1,4.1,2.1), oma = c(0,0,0,0))
barplot(t(relative.comp), 
        border = NA, 
        ylab = "Relative Abundance", 
        ylim = c(0,1), 
        axes = TRUE, 
        col = c("#440154", "#20A387", "#FFBE17", "#1C3B74", "#F94144", "#95D840", "#777B81", "#ADB5BD"), 
        las = 2, 
        cex.names = 0.8, 
        cex.axis = 0.9)
dev.off()

tiff("Relative_Abundance_Phylum_Composite_PPE_ASV_Legend.tiff", width = 5, height = 7, units = 'in', res = 600)
plot.new()
par(mar = c(0,0,0,0), oma = c(0,0,0,0))
legend("center", 
       legend = colnames(relative.comp), 
       cex = 1, 
       ncol = 1, 
       fill = c("#440154", "#20A387", "#FFBE17", "#1C3B74", "#F94144", "#95D840", "#777B81", "#ADB5BD"), 
       x.intersp = 0.1, 
       xjust = 0.1, 
       yjust = 0.3, 
       y.intersp = 1, 
       bty = "n", 
       adj = 0, 
       text.width = 0.1, 
       pt.cex = 0.1)
dev.off()

relative.total <- rltv.Otu.Table(div.tot)
tiff("Relative_Abundance_Phylum_Total_PPE_ASV.tiff", width = 10, height = 8, units = 'in', res = 600)
par(mar = c(5.1,4.1,4.1,2.1), oma = c(0,0,0,0))
barplot(t(relative.total), 
        border = NA, 
        ylab = "Relative Abundance", 
        ylim = c(0,1), 
        axes = TRUE, 
        col = c("#440154", "#20A387", "#FFBE17", "#1C3B74", "#F94144", "#95D840"), 
        las = 2, 
        cex.names = 0.8, 
        cex.axis = 0.9)
dev.off()



uShaDiv <- Tax.sum(uShaASVs, xTXs, 5) %>% as.data.frame()
uShaDiv <- uShaDiv %>% rename(s.nuc = Cryptophyta)
uShaDiv <- uShaDiv %>% 
  mutate(Cryptophyta = rowSums(uShaDiv[c(3, 6)])) %>% 
  select(Chlorophyta, Ochrophyta, Cryptophyta, Haptophyta, Katablepharidophyta)
write.csv2(uShaDiv, file = "C:/Users/Camilo/Dropbox/R/eAnalisis/uShaDiv.csv")
ShaASVs <- wASVs %>% select(all_of(Sha$Seq)) 
ShaASVs <- ShaASVs %>% select(-all_of(uSha$Seq)) #Del total restamos las ASV únicas para el sitio, recordemos que no estan presentes únicamente en Chañaral, por eso presentan abundancia en otros sitios (muestras)
ShaDiv <- Tax.sum(ShaASVs, xTXs, 5) %>% as.data.frame() %>% rename(s.nuc = Cryptophyta)
ShaDiv <- ShaDiv %>% 
  mutate(Cryptophyta = rowSums(ShaDiv[c(3, 6)])) %>% 
  select(Chlorophyta, Ochrophyta, Cryptophyta, Haptophyta, Katablepharidophyta, Rhodophyta)
write.csv2(ShaDiv, file = "C:/Users/Camilo/Dropbox/R/eAnalisis/ShaDiv.csv")
MergedDiv <- ShaDiv %>% 
  mutate(nChl = ShaDiv$Chlorophyta, 
         sChl = uShaDiv$Chlorophyta, 
         nOch = ShaDiv$Ochrophyta, 
         sOch = uShaDiv$Ochrophyta, 
         nCryp = ShaDiv$Cryptophyta, 
         sCryp = uShaDiv$Cryptophyta, 
         nHap = ShaDiv$Haptophyta,
         sHap = uShaDiv$Haptophyta,
         nKata = ShaDiv$Katablepharidophyta, 
         sKata = uShaDiv$Katablepharidophyta, 
         nRho = ShaDiv$Rhodophyta)
MergedDiv <- MergedDiv[, -c(1:6)]
MergedDiv <- rltv.Otu.Table(MergedDiv)
apply(MergedDiv, 2, function(x) sum(x))[1:11]
write.csv2(MergedDiv, file = "C:/Users/Camilo/Dropbox/R/eAnalisis/Relative_SharedMergedDiv.csv", na = "NA")
rMergedDiv <- read.csv2("C:/Users/Camilo/Dropbox/R/eAnalisis/Relative_SharedMergedDiv.csv", header = TRUE, sep = ";", dec = ".", skip = 0, na.strings = "NA")
rownames(rMergedDiv) <- rMergedDiv[, 1]
rMergedDiv <- rMergedDiv[, -1]
rMergedDiv <- as.data.frame(t(rMergedDiv))
rMergedDiv <- rMergedDiv %>% 
  mutate(Division = as.factor(c("Chlorophyta", 
                                "Chlorophyta", 
                                "Ochrophyta", 
                                "Ochrophyta", 
                                "Cryptophyta", 
                                "Cryptophyta", 
                                "Haptophyta",
                                "Haptophyta", 
                                "Katablepharidophyta", 
                                "Katablepharidophyta", 
                                "Rhodophyta")
  )
  )  %>% 
  mutate(Color = c("#3FA03F", "#005F00", "#7A3F8C", "#3B004D", "#3F3F99", "#000052", "#D93F3F", "#AC0000", "#50BFBC", "#007572", "#FFD466")
         )
write.csv2(rMergedDiv, file = "C:/Users/Camilo/Dropbox/R/eAnalisis/rSharedMergedDiv.csv")
bplt <- t(MergedDiv)
bplt <- bplt[, -c(1:36)] #If colSums(bplt) = 1 for each col you can continue

tiff("Relative_Abundance_PPE_Shared.tiff", width = 10, height = 8, units = 'in', res = 600)
par(mar = c(5.1,4.1,4.1,2.1), oma = c(1,1,0,0))
barplot(bplt, border = NA, ylab = "Relative Abundance", ylim = c(0,1), axes = TRUE, col = rMergedDiv$Color, las = 2, cex.names = 0.8, cex.axis = 0.9)
dev.off()

tiff("Relative_Abundance_PuntaChoros_Legend.tiff", width = 5, height = 7, units = 'in', res = 600)
plot.new()
par(mar = c(0,0,0,0), oma = c(0,0,0,0))
legend("center", legend = colnames(MergedDiv), cex = 1.2, ncol = 1, fill = rMergedDiv$Color, x.intersp = 0.2, xjust = 0.2, yjust = 0.3, y.intersp = 1.1, bty = "n", adj = 0, text.width = 0.1, pt.cex = 0.1)
dev.off()

ShaTXs <- filter(xTXs, Seq %in% all_of(uSha$Seq))
write.csv2(ShaTXs, "D:/Documents/GitHub/datasets/datasets/Data/eAnalisis/ShaTXs.csv")

