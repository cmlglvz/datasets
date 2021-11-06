# DF tidy up and management for better visualization
TSEj <- read.csv("/Users/Artemis/Documents/Proyecto_Doctorado/TSEj.csv", header = T, sep = ";", dec = '.', skip = 0)
NeoTSEj <- mutate(TSEj,
                  Muestra = c(paste0("Sample ", 1:43)),
                  Sitio = as.factor(c(rep("Cha", 12), rep("Fla", 12), rep("Hu", 12), rep("Ho", 2), rep("Pc", 5))), 
                  .before = "Chlorophyta")
tidyTSEj <- NeoTSEj %>% gather(key = "Taxa", value = "Abundance", -c(1,2))
plus <- data.frame(matrix(NA, ebar*nlevels(NeoTSEj$Sitio), ncol(NeoTSEj)))
colnames(plus) <- colnames(NeoTSEj)
plus$Sitio <- rep(levels(NeoTSEj$Sitio), each=ebar)
NeoTSEj <- rbind(NeoTSEj, plus)
NeoTSEj <- NeoTSEj %>% arrange(Sitio)
NuTSEj <- mutate(NeoTSEj,
                  id = as.numeric(seq(1, nrow(NeoTSEj))),
                 .before = "Muestra")
# Transform data in a tidy format (long format)
tidyTSEj <- TSEj %>% gather(key = "Taxa", value = colnames(TSEj), -c(1,2))

# Get the name and the y position of each label
ebar <- 3 # Set a number of 'empty bar'
ldata <- NuTSEj
number_of_bar <- nrow(ldata) # calculate the ANGLE of the labels
angle <-  90 - 360 * (ldata$id-0.5) /number_of_bar # substract 0.5 to center the names with their respective bars, extreme right=1 and extreme left=0
# calculate the alignment of labels: right or left
ldata$hjust<-ifelse(angle < -90, 1, 0) # If I am on the left part of the plot, my labels have currently an angle < -90
ldata$angle<-ifelse(angle < -90, angle+180, angle) # flip angle BY to make them readable

# prepare a data frame for base lines
bdata <- NuTSEj %>% 
  group_by(Sitio) %>% 
  summarize(start=min(id), end=max(id) - ebar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
cuadricula <- bdata
cuadricula$end <- cuadricula$end[ c( nrow(cuadricula), 1:nrow(cuadricula)-1)] + 1
cuadricula$start <- cuadricula$start - 1
cuadricula <- cuadricula[-1,]

circulo <- ggplot(NuTSEj, aes(x=as.factor(id), y=Chlorophyta, fill=Sitio)) + # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x = as.factor(id), y = Chlorophyta, fill = Sitio), stat = "identity", alpha = 0.8) + 
  ylim(-15000,50000) + # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  theme_minimal() + # Custom the theme: no axis title and no cartesian grid
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2, 4), "cm")) + # Adjust the margin to make in sort labels are not truncated!
  coord_polar() + # This makes the coordinate polar instead of cartesian.
  #Add the labels, using the label_data dataframe previously created
  geom_text(data = ldata, 
            aes(x = id, 
                y = Chlorophyta+10, 
                label = Muestra, 
                hjust = hjust), 
            color = "black", 
            fontface = "bold", 
            alpha = 0.6, 
            size = 2.5, 
            angle = ldata$angle, 
            inherit.aes = FALSE) + 
  geom_text(data = bdata, 
            aes(x = title, 
                y = -18, 
                label = Sitio), 
            hjust = c(1.25, 1.35, 0, -0.2, -0.55), 
            colour = "black", 
            alpha = 0.5, 
            size = 4, 
            fontface = "bold", 
            inherit.aes = FALSE)
circulo