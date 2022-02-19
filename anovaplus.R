library(tidyverse)
library(vegan)
library(BiodiversityR)
library(ggvegan)
library(ggpubr)
library(ggfortify)
library(pvclust)
library(clustsig)

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

#NMDS
metal.hel <- decostand(metal, method = "hellinger")
nmds1 <- metaMDS(metal.hel, autotransform = FALSE)
ordiplot(nmds1, type = "t")
autoplot(nmds1)
fort <- fortify(nmds1)
p1 <- ggplot() +
  geom_point(data = subset(fort, Score == 'sites'),
             mapping = aes(x = NMDS1, y = NMDS2),
             colour = "black",
             alpha = 0.5) +
  geom_segment(data = subset(fort, Score == 'species'),
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.015, "npc"),
                             type = "closed"),
               colour = "darkgray",
               size = 0,
               alpha = 0) +
  geom_text(data = subset(fort, Score == 'species'),
            mapping = aes(label = Label, x = NMDS1 * 1.1, y = NMDS2 * 1.1),
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
               colour = "darkgray",
               size = 0.8) +
  geom_text(data = subset(fort, Score == 'species'),
            mapping = aes(label = Label, x = NMDS1 * 1.1, y = NMDS2 * 1.1)) +
  
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = 0.8, colour = "gray") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 0.8, colour = "gray") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggarrange(p1, p2, ncol = 1)

summary(metal.env)
adonis(metal ~ Site, data = metal.env)
p3 <- ggplot() +
  geom_point(data = subset(fort, Score == 'sites'),
             mapping = aes(x = NMDS1, y = NMDS2, colour = metal.env$Site),
             alpha = 0.5) +
  geom_segment(data = subset(fort, Score == 'species'),
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.015, "npc"),
                             type = "closed"),
               colour = "darkgray",
               size = 0,
               alpha = 0) +
  geom_text(data = subset(fort, Score == 'species'),
            mapping = aes(label = Label, x = NMDS1 * 1.1, y = NMDS2 * 1.1),
            alpha = 0) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = 0.8, colour = "gray") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 0.8, colour = "gray") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.8, 0.2)) +
  scale_color_discrete("Site")
ggarrange(p3, p2, ncol = 1)

menv <- metal.env %>% mutate_at(c("Vanadium", "Iron", "Niquel", "Copper", "Zinc", "Molybdenum"), ~(scale(.)))
menv
cufe <- menv[, c(3,5)]
tcufe <- as.data.frame(t(cufe))
clstr <- pvclust(tcufe, method.hclust = "complete", method.dist = "correlation", nboot = 1000, parallel = TRUE)
clust <- hclust(cufe, method = "complete")
me.env <- metal.env[, c(3,5)]
sp <- simprof(me.env, method.distance = "braycurtis")
pl.color <- simprof.plot(sp)






