library(tidyverse)
library(vegan)
library(BiodiversityR)
library(ggvegan)
library(ggpubr)
library(ggfortify)

metal <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/shASVs.csv", header = TRUE, sep = ";", dec = ".", skip = 0, fill = TRUE)
metal <- mutate(metal, 
                ID = c("C1A17", "C1F18", "C1A18", "C2A17", "C2F18", "C2A18", "C3A17", "C3F18", "C3A18", "C4A17", "C4F18", "C4A18", "F1A17", "F1F18", "F1A18", "F2A17", "F2F18", "F2A18", "F3A17", "F3F18", "F3A18", "F4A17", "F4F18", "F4A18", "H1A17", "H1F18", "H1A18", "H2A17", "H2F18", "H2A18", "H3A17", "H3F18", "H3A18", "H4A17", "H4F18", "H4A18", "P1F18", "P1A18", "P2F18", "P3F18", "P4F18"), 
                .after = "X")
rownames(metal) <- metal[, 2]
metal <- metal[, -c(1,2)]

metal.env <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/EnvMe.csv", header = TRUE, sep = ";", dec = ".", skip = 0, fill = TRUE)
metal.env <- mutate(metal.env, 
                    ID = c("C1A17", "C1F18", "C1A18", "C2A17", "C2F18", "C2A18", "C3A17", "C3F18", "C3A18", "C4A17", "C4F18", "C4A18", "F1A17", "F1F18", "F1A18", "F2A17", "F2F18", "F2A18", "F3A17", "F3F18", "F3A18", "F4A17", "F4F18", "F4A18", "H1A17", "H1F18", "H1A18", "H2A17", "H2F18", "H2A18", "H3A17", "H3F18", "H3A18", "H4A17", "H4F18", "H4A18", "P1F18", "P1A18", "P2F18", "P3F18", "P4F18"), 
                    .after = "X")
rownames(metal.env) <- metal.env[, 2]
metal.env <- metal.env[, -c(1,2)]

metal.env$Site <- as.factor(metal.env$Site)

long.me <- pivot_longer(data = metal.env, 
                        cols = c("Vanadium", "Iron", "Niquel", "Copper", "Zinc", "Molybdenum"), 
                        names_to = "Metals", 
                        values_to = "Concentration")

long.uno <- pivot_longer(data = metal.env, 
                        cols = !Site, 
                        names_to = "Metals", 
                        values_to = "Concentration")

metal.env$Site <- as.factor(metal.env$Site)
metal.env$Vanadium <- as.integer(metal.env$Vanadium)

cu <- ggplot(metal.env, aes(x = Site, y = Copper)) + geom_boxplot() + geom_jitter(aes(color = Site))
fe <- ggplot(metal.env, aes(x = Site, y = Iron)) + geom_boxplot() + geom_jitter(aes(color = Site))

ggsave("copper_boxplot.tiff", plot = cu)
ggsave("iron_boxplot.tiff", plot = fe)

metal.dist <- vegdist(metal)
attach(metal.env)
metal.sim <- anosim(metal.dist, Site)
summary(metal.sim)
plot(metal.sim)
cu.sim <- anosim(metal.dist, Copper)
summary(cu.sim)
plot(cu.sim)
fe.sim <- anosim(metal.dist, Iron)

pca <- rda(metal.env[,c(3,5)])
summary(pca)
autoplot(pca, arrows = TRUE)

rda <- rda(metal~Copper+Iron, data = metal.env)
rda
ordiplot(rda, type = "text")
