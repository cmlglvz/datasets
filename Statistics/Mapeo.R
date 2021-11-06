library(tidyverse)
library(leaflet)
library(rgdal)
library(rgeos)
library(maptools)
library(maps)
library(broom)
library(htmlwidgets)

map <- leaflet() %>% addTiles()
saveWidget(map, file=paste0(getwd(), "/backgroundMapBasic.html"))
download.file("http://thematicmapping.org/downloads/TM_WORLD_BORDERS_SIMPL-0.3.zip" , destfile = paste0(getwd(), "/world_shape_file.zip"))
system("unzip /Users/Artemis/Documents/Proyecto_Doctorado/world_shape_file.zip")
my_spdf <- readOGR( 
  dsn = "/Users/Artemis/Documents/Proyecto_Doctorado/", 
  layer = "TM_WORLD_BORDERS_SIMPL-0.3",
  verbose = FALSE
)
#Base version
par(mar=c(0,0,0,0))
plot(my_spdf, col = "#002C3D", bg = "#EDF2F4", lwd = 0.25, border = 0)
#ggplot2 version (recommended)
spdf_fortified <- tidy(my_spdf, region = "NAME")
ggplot() +
  geom_polygon(data = spdf_fortified, aes(x = long, y = lat, group = group), fill = "#002C3D", color = "#EDF2F4") +
  theme_void()

mapdos <- leaflet() %>% 
  addTiles() %>% 
  setView(lng = -70.667222, lat = -26.224278, zoom = 6) %>% 
  addProviderTiles("Stamen.Toner") %>% 
  addCircleMarkers(lng = -70.667222, lat = -26.224278, color = "#F77F00", radius = 30, opacity = 1, label = "Cha") %>% 
  addCircleMarkers(lng = -70.7, lat = -26.548611, color = "#FCB740", radius = 30, opacity = 1, label = "Fla") %>% 
  addCircleMarkers(lng = -71.210278, lat = -28.336944, color = "#023E8A", radius = 30, opacity = 1, label = "Hu") %>% 
  addCircleMarkers(lng = -71.4975, lat = -29.223333, color = "#57E3FF", radius = 30, opacity = 1, label = "Hu") %>% 
  addMeasure(position = "bottomright", primaryLengthUnit = "kilometers") %>% 
  addScaleBar(position = "bottomleft")
mapdos
saveWidget(mapdos, file=paste0(getwd(), "/bgMapAll.html"))

Chile <- map_data("world")%>% filter(region == "Chile")
Ciudades <- world.cities %>% filter(country.etc == "Chile")
ggplot() +
  geom_polygon(data = Chile, aes(x = long, y = lat, group = group), fill = "grey", alpha = 0.3) + 
  geom_point(data = Ciudades, aes(x = long, y = lat)) + 
  theme_void() + 
  coord_map() 

