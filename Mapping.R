library(tidyverse)
library(sf)
library(tmap)    # for static and interactive maps
#library(leaflet) # for interactive maps
#install.packages("osmdata")
#install.packages("osmplotr")
#library(osmdata)
#library(osmplotr)
library(spData)
library(sp)

theme_set(theme_classic())
#https://www.youtube.com/watch?v=WsbkVkBLkdQ (Mapping change in samples and osmplotr/osmdata info)
#https://bookdown.org/nicohahn/making_maps_with_r5/docs/tmap.html

#tmap.csv can be just site co-ordinates or you could add columns with site parameters if mapping change in samples
map_bug <- read_csv("tmap.csv") %>%
  filter(Site %in% c("EAS1", "FRY1", "SPC1", "LLW1", "CRO1", "ROL1"))

bugs2 <-st_as_sf(map_bug, coords = c("x","y"))

#Create bounding box based off of coordinates of sites
my_bbox <- c(xmin = min(map_bug$x) - 0.4 ,
     ymin = min(map_bug$y) - 0.1,
     xmax = max(map_bug$x) + 0.6 ,
     ymax = max(map_bug$y) +0.2)

my_bbox

tmap_mode("plot")

#### Map for inlay

#create box for region the map represents
#source: https://geocompr.robinlovelace.net/adv-map.html
region = st_bbox( c(xmin =   -82.61103, ymin = 36.90756,
                    xmax = -81.63782, ymax =37.74699))  %>%
  st_as_sfc()

region

# Represent appalachian region on map
# source: https://hub.arcgis.com/datasets/ARCgov::arc-boundary-2008/explore?location=37.375411%2C-82.619823%2C5.99
appalachia <- st_read("app_coal_regions/app_coal_regions.shp")

#Map of US with region of sampling represented
state_bbox <- c(xmin =   -89.782417, ymin = 27.991440,
                xmax = -71.196289, ymax = 45.519394)

 map_us <- tm_shape(us_states, bbox = state_bbox)+
  tm_borders( alpha = 0.5)+
    tm_shape(appalachia) +
    tm_fill(alpha = 0.5, col = "grey36") +
   tm_shape(region) + 
   tm_borders(lwd = 3, alpha = 1, col = "black") 
  
map_us

### site map

# to add state names add coordinates with state name and add as layer
# add to x = goes west, add to y = goes N
statenames <- tribble(~"x", ~"y", ~"name",
                      -81.7, 36.989084, "Virginia",
                      -81.25,37.45, "West Virginia",
                      -82.6, 37.40000, "Kentucky",
                      ) %>%
  st_as_sf(coords = c("x","y"), crs = 4326)

# Represent counties
#51 = VA, 54 = WV; source: https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
counties <- st_read("cb_2018_us_county_5m/cb_2018_us_county_5m.shp") %>%
  filter(COUNTYNS %in% c("01550061", "01496656", "01497431", "01497376","01689162", "	
01209192", "01550029", "01550009" , "01550047", "01497573")) 

#add state lines
States <- st_read("cb_2018_us_state_5m/cb_2018_us_state_5m.shp") 

# add rivers; https://www.weather.gov/gis/Rivers
rivers <- st_read("rivers/rv16my07.shp") %>%
  filter(TYPE %in% c("O", "P", "Q","R"))

#Map it
sitemap <- tm_shape(rivers, bbox = my_bbox) +
  tm_lines(col = "deepskyblue") +
  tm_shape(States) +
  tm_borders(lwd = 3)+
  tm_text("NAME") +
tm_shape(counties)+
    tm_borders()+
    tm_text("NAME", size = 1) +
    tm_shape(appalachia) +
    tm_fill(alpha = 0.4, col = "grey36") +
  tm_shape(bugs2, projection = "longlat") + 
    tm_symbols(size = 1, col = "Stream Type") +
  tm_shape(statenames) +
  tm_text("name", size= 2) +
  tm_layout(legend.position = c("RIGHT", "bottom"),
            legend.text.size = 1.3,
            legend.title.size = 1.6,
              inner.margins = 0.1,
              title.position = c("center", "TOP"),
            title.size = 25) +
  tm_scale_bar(position = c("left", "BOTTOM"))+
  tm_compass(position = c("right", "bottom"))
   
sitemap

# function to add zoomed out map to site map 
#print(map_us, vp = grid::viewport(0.25, 0.8, width = 0.3, height = 0.6))

tmap_save(sitemap, filename = "sitemap.png",
          dpi=100, insets_tm = map_us, 
          insets_vp = grid::viewport(0.15, 0.76, width = 0.25, height = 0.5))


# objects for experimenting with map layers
#roads <- extract_osm_objects(key = "highway", bbox = my_bbox, sf = T)
#waterways <- extract_osm_objects(key = "waterway", bbox = my_bbox, sf = T)
