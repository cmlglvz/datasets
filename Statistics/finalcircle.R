library(tidyverse)
library(viridis)
library(RColorBrewer)

# DF tidy up and management for better visualization
TSEj <- read.csv("/Users/Artemis/Documents/Proyecto_Doctorado/TSEj.csv", header = T, sep = ";", dec = '.', skip = 0)
snames <- TSEj[,1]
TSEji <- TSEj[, -1]
TSEj <- apply(TSEji, 2, function(x) as.numeric(as.integer(x))) %>% as.data.frame()
rownames(TSEj) <- snames
TSSelect <- select(TSEj, Chlorophyta, Ochrophyta, Metazoa, Cryptophyta, Cryptophyta.nucl, Haptophyta, Katablepharidophyta, Rhodophyta, Stramenopiles_X, Alveolata_X, Cercozoa) 
TSSelect <- mutate(TSSelect, Total = rowSums(TSSelect))
TSEj <- TSSelect %>% mutate(Muestra = c(paste0("Sample ", 1:43)), 
                        Sitio = as.factor(c(rep("Cha", 12), rep("Fla", 12), rep("Hu", 12), rep("Ho", 2), rep("Pc", 5))), 
                        .before = "Chlorophyta") %>% 
  mutate(Chlorophyta = Chlorophyta/Total*100, 
         Ochrophyta = Ochrophyta/Total*100, 
         Metazoa = Metazoa/Total*100, 
         Cryptophyta = Cryptophyta/Total*100,
         Cryptophyta.nucl = Cryptophyta.nucl/Total*100,
         Haptophyta = Haptophyta/Total*100,
         Katablepharidophyta = Katablepharidophyta/Total*100,
         Rhodophyta = Rhodophyta/Total*100, 
         Stramenopiles_X = Stramenopiles_X/Total*100, 
         Alveolata_X = Alveolata_X/Total*100, 
         Cercozoa = Cercozoa/Total*100) %>%
  gather(key = "Taxa", value = "Abundance", -c(1,2))

ebar <- 2 # Set a number of 'empty bar'
n0bsType <- nlevels(as.factor(TSEj$Taxa))
plus <- data.frame(matrix(NA, ebar*nlevels(TSEj$Sitio)*n0bsType, ncol(TSEj)))
colnames(plus) <- colnames(TSEj)
plus$Sitio <- rep(levels(TSEj$Sitio), each = ebar*n0bsType)
TSEj <- rbind(TSEj, plus)
TSEj <- TSEj %>% arrange(Sitio, Muestra)
TSEj$ID <- rep(seq(1, nrow(TSEj)/n0bsType), each = n0bsType)

# Get the name and the y position of each label
label_data <- TSEj %>% group_by(ID, Muestra) %>% summarize(tot = sum(Abundance))
number_of_bar <- nrow(label_data) # calculate the ANGLE of the labels
angle <-  90 - 360 * (label_data$ID-0.5) /number_of_bar # substract 0.5 to center the names with their respective bars, extreme right=1 and extreme left=0
# calculate the alignment of labels: right or left
label_data$hjust <- ifelse(angle < -90, 1, 0) # If I am on the left part of the plot, my labels have currently an angle < -90
label_data$angle <- ifelse(angle < -90, angle+180, angle) # flip angle BY to make them readable

# prepare a data frame for base lines
base_data <- TSEj %>% 
  group_by(Sitio) %>% 
  summarize(start = min(ID), end = max(ID) - ebar) %>% 
  rowwise() %>% 
  mutate(title = mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1, ]

# Make the plot (extra colours "#0B525B", "#E9D8A6", "#8943C7", "#75EC53", "#3FB7F8", "#BB3E03", "#B0AADA", "#E9F961", "#778899", "#00609F", "#F0FFF0", "#FFDAB9", "#7FFFD4", "#B0F961", "#F5005A", "Hilomonadea" = "#808000", "Apusomonadidae" = "#DCDCDC", "Lobosa" = "#AEAE6E", "Metamonada" = "#9B2226", "Foraminifera" = "#FFFFF0") "#3E1F47", "#FF1493", "#ADFF2F", "#FF4500", "#191970", "#B0AADA", "#E0FFFF", "#D6A218", "#94D2BD", "#638C50", "#E96E9B", "#DBD18C", "#7D8CB8", "#F0A280", "#048BA8", "#E78888", "#FFFACD" / scale_fill_manual(values = c("Chlorophyta" = "#0AFF0A", "Ochrophyta" = "#0000F5", "Metazoa" = "#20E99F", "Cryptophyta" = "#FF4500", "Haptophyta" = "#FF0000", "Katablepharidophyta" = "#54478C", "Cryptophyta.nucl" = "#FF00FF", "Rhodophyta" = "#002C3D", "Stramenopiles_X" = "#FA8072", "Alveolata_X" = "#E9D8A6", "Cercozoa" = "#3FB7F8")))
circulo <- ggplot(TSEj) +
  # Add the stacked bar
  geom_bar(aes(x = as.factor(ID), y = Abundance, fill = Taxa), stat = "identity", alpha = 1) + 
  scale_fill_viridis(discrete = TRUE) +
  ylim(-15000, max(label_data$tot, na.rm = TRUE)) + 
  theme_minimal() + 
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(), 
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  ) + 
  coord_polar() + 
  # Add labels on top of each bar
  geom_text(data = label_data, aes(x = ID, y = tot+10, label = Muestra, hjust = hjust), color = "black", fontface = "bold", alpha = 0.6, size = 5, angle = label_data$angle, inherit.aes = FALSE) + 
  # Add base line information
  geom_segment(data = base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha = 0.8, size = 0.6, inherit.aes = FALSE) + 
  geom_text(data = base_data, aes(x = title, y = -18, label = Sitio), hjust = c(1.25, 1.35, 0, -0.2, -0.55), colour = "black", alpha = 0.8, size = 4, fontface = "bold", inherit.aes = FALSE)
circulo

# Save at png
ggsave(circulo, file = "output.png", width = 15, height = 15, dpi = 600)
