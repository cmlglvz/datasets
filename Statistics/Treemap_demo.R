library(treemap)
library(tidyverse)
library(RColorBrewer)
library(d3treeR)
library(htmlwidgets)

fTXs <- read.csv("/Users/Artemis/Documents/Proyecto_Doctorado/Pipeline/fTXs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
fTXs <- fTXs[, -1]
iASVs <- read.csv("/Users/Artemis/Documents/Proyecto_Doctorado/Pipeline/iASVs.csv", header = T, sep = ";", dec = '.', skip = 0) %>% 
  rename(OTU = X)
iASVs <- mutate(iASVs, Total = rowSums(iASVs[2:44]))
iTVs <- inner_join(iASVs, fTXs, by = "OTU") %>% relocate(Seq, .after = OTU)
iSN <- iTVs[,1]
iTVs <- iTVs[,-1]
rownames(iTVs) <- iSN
write.csv2(iTVs, file = "/Users/Artemis/Documents/Proyecto_Doctorado/Treemap/iTVs.csv")
kTVs <- iTVs %>% select(Division, Class, Order, Family, Genus, Species, Total)
#%>% group_by(Division) %>% summarise(Total = sum(Total)) %>% as.data.frame()
png("Treemap_Division_PPE_ASV.png", width = 8, height = 7, units = 'in', res = 600)
treemap(kTVs, 
        index = "Division", 
        vSize = "Total", 
        type = "index", 
        palette = c("#064c55", "#d00000", "#9d0208", "#f8961e", "#20704a", "#2cb58e", "#40b557", "#ded2ba", "#f1c453")
        )
dev.off()

kols <- colorRampPalette(c("#064c55", "#d00000", "#9d0208", "#f8961e", "#20704a", "#2cb58e", "#40b557", "#ded2ba", "#002b3d", "#f1c453"))(45)
p <- treemap(kTVs,
             index = c("Division","Class"),
             vSize = "Total",
             type = "index",
             palette = kols,
             bg.labels = c("white"),
             align.labels = list(c("center", "center"), 
                                 c("right", "bottom")
                                 )
             )
inter <- d3tree2(p, rootname = "General")
saveWidget(inter, file = paste0(getwd(), "/HtmlWidget/interactiveTreemap.html"))

pTVs <- kTVs %>% filter(Division == "Chlorophyta")
png("Treemap_Class_Chlorophtya_ASV.png", width = 8, height = 7, units = 'in', res = 600)
treemap(pTVs,
        index = "Class", 
        vSize = "Total", 
        type = "index", 
        palette = c("#54478c", "#2c699a", "#048ba8", "#0db39e", "#16db93", "#83e377", "#b9e769", "#efea5a", "#f1c453", "#f29e4c")
        )
dev.off()

hASVs <- read.csv("/Users/Artemis/Documents/Proyecto_Doctorado/Pipeline/iASVs.csv", header = T, sep = ";", dec = '.', skip = 0)
hSN <- hASVs[, 1]
hASVs <- hASVs[, -1]
rownames(hASVs) <- hSN
hASVs <- as.matrix(hASVs)
heatmap(hASVs)

#Modelos
group <- c(rep("group-1",4),rep("group-2",2),rep("group-3",3))
subgroup <- paste("subgroup" , c(1,2,3,4,1,2,1,2,3), sep="-")
value <- c(13,5,22,12,11,7,3,1,23)
data <- data.frame(group,subgroup,value)

# Custom labels:
treemap(data, index=c("group","subgroup"),     vSize="value", type="index",
        
        fontsize.labels=c(15,12),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("white","orange"),    # Color of labels
        fontface.labels=c(2,1),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels=c("transparent"),              # Background color of labels
        align.labels=list(
          c("center", "center"), 
          c("right", "bottom")
        ),                                   # Where to place labels in the rectangle?
        overlap.labels=0.5,                      # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        inflate.labels=F,                        # If true, labels are bigger when rectangle is bigger.
        
)


data2 <- as.matrix(mtcars)














