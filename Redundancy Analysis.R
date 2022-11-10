#Redundancy Analysis
# https://r.qcbs.ca/workshop10/book-en/exploration.html
library(vegan)
library(dplyr)
library(tidyverse)

com <- read.csv("bugid.csv")  %>%
  filter(Season == "Spring")%>%
  #filter(Site != "ROL1" ) %>% # For Spring
  select(-c(Season, Site))
com <- com[, colSums(com != 0, na.rm = TRUE) > 0]

chem <- read.csv("chem.f21-s22.notrib.csv") %>%
  filter(season == "Spring") %>%
  #filter(Site != "SPC6" ) %>% #For Fall
  select(  do.mgl:so4.hco3) 
chem <- chem[, colSums(chem != 0, na.rm = TRUE) > 0] 
springrow <-  read.csv("chem.f21-s22.notrib.csv") %>%
  filter(season == "Spring") %>%
  #filter(Site != "SPC6" ) %>% #For Fall
  select( Stream) 

hab <- read.csv("habitatmasteradj.csv") %>%
  #filter(Site != "EAS1" & Site != "EAS9" & Site !="LLW6" ) %>%
  select(-c(Site:smallcobble)) 
  

#Examine Community Data

sum(com == 0)
# Calculate proportion of zeros in the dataset
sum(com == 0)/(nrow(com) * ncol(com))

# Apply Hellinger transformation to correct for the double
# zero problem 77% zeros, double zero problem
spe.hel <- decostand(com, method = "hellinger")

# Examine environmental data
# We can visually look for correlations between variables:
heatmap(abs(cor(chem)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
# Scale and center variables
chem.z <- decostand(chem, method = "standardize",  na.rm =TRUE)
chem.z <- chem.z %>%
  select_if(~ ! any(is.na(.)))

# Variables are now centered around a mean of 0
round(apply(chem.z, 2, mean), 1)
apply(chem.z, 2, sd)

chem.z$Stream <- springrow$Stream
# Again for habitat
habheat <- heatmap(abs(cor(hab)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

habheat
# Scale and center variables
hab.z <- decostand(hab, method = "standardize", na.rm =TRUE)
hab.z <- hab.z[, colSums(hab.z != 0, na.rm = TRUE) > 0]

# Variables are now centered around a mean of 0
round(apply(hab.z, 2, mean), 1)
apply(hab.z, 2, sd)

#Begin RDA
# Model the effect of all environmental variables on fish
# community composition

spe.rda <- rda(spe.hel ~ ., data = chem.z)
summary(spe.rda.signif)
test <- summary(spe.rda.signif)
centroids <- as.data.frame(test$centroids)
centroids <- rownames_to_column(centroids, "VALUE")

ggplot(centroids, aes(RDA1, RDA2, color = VALUE))+
  geom_point()

# Forward selection of variables:
fwd.sel <- ordiR2step(rda(spe.hel ~ 1, data = chem.z), # lower model limit (simple!)
                      scope = formula(spe.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) # change to TRUE to see the selection process!

fwd.sel$call
spe.rda.signif <- rda(formula = spe.hel ~ Stream + ca.mg + ba.ugl + mg.mgl + sc.uScm + 
      co.ugl + u.ugl + cd.ugl + ni.ugl + al.ugl + temp.C + cl.mgl, 
    data = chem.z) # when grouped by stream

#spe.rda.signif <- rda(formula = spe.hel ~ ba.ugl + v.ugl + ti.ugl + na.mgl + so4.hco3 + 
                        k.mgl + ca.mg + cl.mgl + ac.uScm + sr.ugl + mn.ugl + npoc.mgl, 
                      data = chem.z)

#spe.rda.signif <- rda(formula = spe.hel ~ pfines + D50 + avgcancov + avgwetwidth + avgembedd, data = hab.z)

spe.rda.signif <- rda(formula = spe.hel ~ pfines + D16 + avgcancov + D84 + avgwetwidth, data = hab.z) # habitat adj fines 2-4 and largecobble <4000

#spe.rda.signif <- rda(formula = spe.hel ~ ppebbles + pfines + D16, data = hab.z)# habitat adj fines 2-4 and largecobble <bedrock

anova.cca(spe.rda.signif, step = 1000)
RsquareAdj(spe.rda.signif)

# Type 1 scaling
ordiplot(spe.rda.signif, scaling = 1, type ='text')

# Type 2 scaling
ordiplot(spe.rda.signif, scaling = 2, type = "text")

# Custom triplot code!

## extract % explained by the first 2 axes
perc <- round(100*(summary(spe.rda.signif)$cont$importance[2, 1:2]), 2)

## extract scores - these are coordinates in the RDA space
sc_si <- scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(spe.rda.signif, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(spe.rda.signif, display="bp", choices=c(1, 2), scaling=1)

## Custom triplot, step by step

# Set up a blank plot with scaling, axes, and labels
plot(spe.rda.signif,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-1,1), 
     ylim = c(-1,1),
     # label the plot (title, and axes)
     main = "Triplot RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)

# add points for site scores
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = "steelblue", # fill colour
       cex = 1.2) # size
# add points for species scores
points(sc_sp, 
       pch = 22, # set shape (here, square with a fill colour)
       col = "black",
       bg = "#f2bd33", 
       cex = 1.2)
points("sites", color = Stream)
# add text labels for species abbreviations
text(sc_sp + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = rownames(sc_sp), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)
# add arrows for effects of the expanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp[,1], sc_bp[,2], # end them at the score value
       col = "red", 
       lwd = 3)
# add text labels for arrows
text(x = sc_bp[,1] -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 1, 
     font = 2)

#### Sokol's Variation partitioning tutorial
library("geosphere")
library(sp)

spring.xy <- read.csv("spring_coord_notrib.csv") %>%
  select(-c("Site"))

site.loc.sp = sp::SpatialPointsDataFrame(coords = data.frame(spring.xy$x, 
                                                             spring.xy$y),
                                         data = spring.xy,
                                         proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

# projected to WGS 1984 UTM N 17
site.loc.sp = spTransform(site.loc.sp, 
                          CRS("+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
site.loc.sp$x_utm = site.loc.sp@coords[,1] 
site.loc.sp$y_utm = site.loc.sp@coords[,2]

# could have also converted directly to sf

# check formatting and projection
str(site.loc.sp)
#View(head(site.loc.sp@data))
plot(site.loc.sp)

site.loc.sf = st_as_sf(site.loc.sp, coords = c("long", "lat"), # convert sp to sf
                       crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

site.xy = as.data.frame(site.loc.sf) # convert sf to df, needed to remove point geometry for writing to .txt
site.xy = site.xy %>%
  select( x_utm,
         y_utm) 

spring.xy <- site.xy %>%
  select(c(x_utm, y_utm))

mod.pcnm <- pcnm( dist(spring.xy) )
vectors.pcnm <- data.frame(mod.pcnm$vectors)

ordisurf(spring.xy, scores(mod.pcnm, choi=1), bubble = 4, main = "PCNM 1")

ordisurf(spring.xy, scores(mod.pcnm, choi=13), bubble = 4, main = "PCNM 13")

com.hel <- decostand( com, "hel")

# Space
d.space <- data.frame( spring.xy, vectors.pcnm)
#combine all spatial variables, x, y, and PCNMs

d.space.scaled <- data.frame( scale(d.space) )
#center spatial variables on 0, and standardize
# null model with intercept
mod.0 <- rda( com.hel ~ 1, data = d.space.scaled)

# model with all spatial variables included
mod.1 <- rda( com.hel ~ ., data = d.space.scaled)

 #stepwise selection of the best model
mod.best <- ordiR2step(mod.0, scope = mod.1 )
summary(mod.best)

S.keepers <- names( mod.best$terminfo$ordered )

# Chem
mod.0 <- rda(com.hel ~ 1, data = chem.z)
mod.1 <- rda(com.hel ~ ., data = chem.z)

#stepwise selection of the best model
mod.best <- ordiR2step(mod.0, scope = mod.1 )
E.keepers <- names(mod.best$terminfo$ordered)
E.keepers

#Habitat
mod.0 <- rda(com.hel ~ 1, data = hab.z)
mod.1 <- rda(com.hel ~ ., data = hab.z)

#stepwise selection of the best model
mod.best <- ordiR2step(mod.0, scope = mod.1 )
H.keepers <- names(mod.best$terminfo$ordered)
H.keepers

# put them all together
d.C <- chem.z[,E.keepers]
d.S <- d.space.scaled[,S.keepers]
d.H <- hab.z[,H.keepers]
mod.varpart <- varpart(com.hel, d.C, d.H)

png("Pcoa.C&H.allstream.spring.png", width = 600, height = 500)
plot(mod.varpart)
title(main = "All Streams (Spring 2022)", sub = "X1 = Water Quality, X2 = Habitat")
dev.off()


library(vegan)
library(ggplot2)
library()
install.packages("remotes")
remotes::install_github("jfq3/ggordiplots")
data("dune")
data("dune.env")
dune.hel <- decostand(dune, method = "hellinger")
ord <- rda(dune.hel)
ggordiplots::gg_ordiplot(ord, groups = dune.env$Management, pt.size = 3)

ord <- rda(spe.hel)
ggordiplots::gg_ordiplot(ord, groups = chem.z$Stream)
ggordiplots::gg_ordiplot()
library("digest")

