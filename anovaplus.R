library(tidyverse)
library(readxl)
library(glue)
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
lookup <- mutate(metal.env, 
                 Samples = c("C1A17", "C1F18", "C1A18", "C2A17", "C2F18", "C2A18", "C3A17", "C3F18", "C3A18", "C4A17", "C4F18", "C4A18", "F1A17", "F1F18", "F1A18", "F2A17", "F2F18", "F2A18", "F3A17", "F3F18", "F3A18", "F4A17", "F4F18", "F4A18", "H1A17", "H1F18", "H1A18", "H2A17", "H2F18", "H2A18", "H3A17", "H3F18", "H3A18", "H4A17", "H4F18", "H4A18", "P1F18", "P1A18", "P2F18", "P3F18", "P4F18"), 
                 .after = "X")
rownames(metal.env) <- metal.env[, 2]
metal.env <- metal.env[, -c(1,2)]
lookup <- lookup[, -1]
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

metal.dist <- vegdist(metal, method = "bray")
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

metal.hel <- decostand(x = metal, method = "hellinger")
metal.dist <- vegdist(x = metal.hel, method = "bray")
write.csv2(x = as.matrix(metal.dist), file = "metal_dist.csv")
pcoa <- cmdscale(metal.dist, eig = TRUE, add = TRUE)
positions <- pcoa$points
positions
percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)
pround <- format(round(percent_explained[1:2], digits = 1), nsmall = 1, trim = TRUE)
colnames(positions) <- c("PCoA1", "PCoA2")
labs <- c(glue("PCo 1({pround[1]}%)"), 
          glue("PCo 2({pround[2]}%)")
          )
positions %>% as_tibble(rownames = "samples")
ggplot(positions, aes(x = PCoA1, y = PCoA2)) + 
  geom_point() + 
  labs(x = labs[1], y = labs[2])

tibble(pe = cumsum(percent_explained), 
       axis = 1:length(percent_explained)) %>% 
  ggplot(aes(x = axis, y = pe)) + 
  geom_line()

dist.df <- read.csv2(file = "metal_dist.csv", header = TRUE, sep = ";", dec = ",", skip = 0, row.names = 1)
metadata <- inner_join(metal.env, dist.df, by = "ID")
metadata
dist.df <- as.matrix(metal.dist) %>% as.data.frame()

vdist <- vegdist(metal.hel, method = "bray", diag = TRUE)
npmanova <- adonis(vdist~Site, data = metal.env, permutations = 1000)
str(npmanova)
npmanova$aov.tab$`Pr(>F)`[1]
p.adonis <- npmanova[["aov.tab"]][["Pr(>F)"]][1] #lo mismo que el método anterior

#Otro modelo
oth.model <- adonis(vdist~Iron*Copper, 
                    data = metal.env, 
                    permutations = 1000)


#Subset
pairwise_p <- numeric()

chafla.dist <- vegdist(metal.hel[c(1:24),], method = "bray", diag = TRUE)
chafla.env <- metal.env[c(1:24),]
chafla.test <- adonis(chafla.dist~Site, 
                 data = chafla.env, 
                 permutations = 1000)
pairwise_p["Cha_Fla"] <- chafla.test[["aov.tab"]][["Pr(>F)"]][1]

hupc.dist <- vegdist(metal.hel[c(25:41),], method = "bray", diag = TRUE)
hupc.env <- metal.env[c(25:41),]
hupc.test <- adonis(hupc.dist~Site, 
                      data = hupc.env, 
                      permutations = 1000)
pairwise_p["Hu_Pc"] <- hupc.test[["aov.tab"]][["Pr(>F)"]][1]

chahu.dist <- vegdist(metal.hel[c(1:12, 25:36),], method = "bray", diag = TRUE)
chahu.env <- metal.env[c(1:12, 25:36),]
chahu.test <- adonis(chahu.dist~Site, 
                      data = chahu.env, 
                      permutations = 1000)
pairwise_p["Cha_Hu"] <- chahu.test[["aov.tab"]][["Pr(>F)"]][1]

flapc.dist <- vegdist(metal.hel[c(13:24, 37:41),], method = "bray", diag = TRUE)
flapc.env <- metal.env[c(13:24, 37:41),]
flapc.test <- adonis(flapc.dist~Site, 
                      data = flapc.env, 
                      permutations = 1000)
pairwise_p["Fla_Pc"] <- flapc.test[["aov.tab"]][["Pr(>F)"]][1]

p.adjust(pairwise_p, method = "BH") < 0.05

set.seed(19910420)
nmds <- metaMDS(vdist)
str(nmds)
nmds$points
scores(nmds) %>% #funcion vegan similar al objeto anterior 
  as_tibble(rownames = "Samples") %>% 
  inner_join(., lookup, by = "Samples") %>% 
  mutate(Cu.Level = if_else(Copper > 10, "Higher", "Lower")) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = Cu.Level, size = Iron)) + 
  geom_point(aes(shape = Site)) + 
  coord_fixed() + 
  scale_color_manual(name = NULL, 
                     breaks = c("Higher", "Lower"), 
                     values = c("#D90429", "#34A0A4"), 
                     labels = c("Higher Cu concentration", 
                                "Lower Cu concentration")
                     ) + 
  scale_shape_manual(name = NULL, 
                     breaks = c("Cha", "Fla", "Hu", "Pc"),
                     values = c(0, 1, 2, 21), 
                     labels = c("Chañaral", 
                                "Flamenco", 
                                "Huasco", 
                                "Punta de Choros")
                     ) + 
  theme_minimal()
