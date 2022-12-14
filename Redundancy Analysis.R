#Redundancy Analysis
# https://r.qcbs.ca/workshop10/book-en/exploration.html
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyverse)

com <- read.csv("bugid.csv")  %>%
  filter(Season == "Spring") %>%
  #filter(Stream == "CRO") %>%
  select(-c(Season, Site, Stream))
com <- com[, colSums(com != 0, na.rm = TRUE) > 0]

chem <- read.csv("chem.f21-s22.notrib.reduced.csv")%>%
  filter(season == "Spring") %>%
  #filter(Stream == "ROL (MI)") %>%
  #filter(Stream == "CRO") %>%
  select(-c("Stream", Site:season))
chem <- chem[, colSums(chem != 0, na.rm = TRUE) > 0] 
springrow <-  read.csv("chem.f21-s22.notrib.csv") %>%
  filter(season == "Spring") %>%
  select(Stream) 

hab <- read.csv("habitatmaster.csv") %>%
  filter(Stream == "LLW (MI)") %>%
  select(-c(Stream:smallcobble)) 
  

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
chem <- chem %>% select(do.mgl:so4.hco3)
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
hab <- select(hab, -c(Stream))
hab.z <- decostand(hab, method = "standardize", na.rm =TRUE)
hab.z <- hab.z[, colSums(hab.z != 0, na.rm = TRUE) > 0]

# Variables are now centered around a mean of 0
round(apply(hab.z, 2, mean), 1)
apply(hab.z, 2, sd)

hab.z$Stream <- springrow$Stream

env <- left_join(chem, hab, by = "Site") %>%
  select(-c(Stream.x:season, Stream.y:smallcobble))
env.z <- decostand(env, method = "standardize", na.rm =TRUE)
env.z <- env.z[, colSums(env.z != 0, na.rm = TRUE) > 0]

round(apply(env.z, 2, mean), 1)
apply(env.z, 2, sd)
# row 32 NA needs removed
env.z <- env.z[-c(32),]
spe.hel <- spe.hel[-c(32),]


# VARIATION PARTITIONING
#### Sokol's Variation partitioning tutorial
library("geosphere")
library(sp)
library(sf)

spring.xy <- read.csv("spring_coord_notrib.csv") %>%
#filter(Site != "EAS9") %>%
#filter(Site != "FRY8")%>%
#filter(Site != "FRY9")%>%
#filter(Site != "ROL7") %>%
  #filter(Site != "SPC6") %>%
  filter(Stream == "SPC (MI)") %>%
  select(-c("Stream","Site"))

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

#specified knts, implications?
ordisurf(spring.xy, scores(mod.pcnm, choi=1), bubble = 4,knots = 5, main = "PCNM 1")

ordisurf(spring.xy, scores(mod.pcnm, choi=13), bubble = 4, main = "PCNM 13")

com <- read.csv("bugid.csv")  %>%
  filter(Season == "Spring") %>%
  filter( Stream == "SPC (MI)") %>%
  select(-c(Stream:Site))
com.hel <- decostand( com, "hel")

# Space
d.space <- data.frame( spring.xy, vectors.pcnm)
#combine all spatial variables, x, y, and PCNMs

d.space.scaled <- data.frame( scale(d.space) )
#center spatial variables on 0, and standardize
# null model with intercept
mod.0 <- rda( com.hel ~ 1, data = d.space.scaled)
plot(mod.1)
anova.cca(mod.1)

# model with all spatial variables included
mod.1 <- rda( com.hel ~ ., data = d.space.scaled)

 #stepwise selection of the best model
mod.best <- ordiR2step(rda( com.hel ~1, data = d.space.scaled), scope = formula(mod.1),  
                       direction = "forward",
                       R2scope = FALSE, # can't surpass the "full" model's R2
                       pstep = 1000,
                       trace = FALSE)
summary(mod.best)

plot(mod.best)

S.keepers <- names( mod.best$terminfo$ordered )
S.keepers

# Chem
chem <- read.csv("chem.f21-s22.notrib.reduced.csv")%>%
  filter(season == "Spring") %>%
  filter( Stream == "SPC (MI)") %>%
  select(do.mgl:so4.hco3)
chem <- chem[, colSums(chem != 0, na.rm = TRUE) > 0] 
chem.z <- decostand(chem, method = "standardize",  na.rm =TRUE)
chem.z <- chem.z %>%
  select_if(~ ! any(is.na(.)))

# Variables are now centered around a mean of 0
round(apply(chem.z, 2, mean), 1)
apply(chem.z, 2, sd)

#chem.z$Stream <- springrow$Stream

mod.0 <- rda(com.hel ~ 1, data = chem.z)
mod.1 <- rda(com.hel ~ ., data = chem.z)

#stepwise selection of the best model
mod.best <- ordiR2step(mod.0, scope = mod.1, R2scope = FALSE)
C.keepers <- names(mod.best$terminfo$ordered)
C.keepers

#Habitat
hab <- read.csv("habitatmaster.csv") %>%
  filter( Stream == "SPC (MI)") %>%
select(-c(Stream:smallcobble)) 
hab.z <- decostand(hab, method = "standardize", na.rm =TRUE)
hab.z <- hab.z[, colSums(hab.z != 0, na.rm = TRUE) > 0]

# Variables are now centered around a mean of 0
round(apply(hab.z, 2, mean), 1)
apply(hab.z, 2, sd)

#hab.z$Stream <- springrow$Stream

mod.0 <- rda(com.hel ~ 1, data = hab.z)
mod.1 <- rda(com.hel ~ ., data = hab.z)

#stepwise selection of the best model
mod.best <- ordiR2step(mod.0, scope = mod.1 )

mod.best <- ordiR2step(rda( com.hel ~1, data = hab.z), scope = formula(mod.1),  
                       direction = "forward",
                       R2scope = FALSE, # can't surpass the "full" model's R2
                       pstep = 1000,
                       trace = FALSE)
H.keepers <- names(mod.best$terminfo$ordered)
H.keepers

# put them all together
d.C <- chem.z[,C.keepers]
d.S <- d.space.scaled[,S.keepers]
d.H <- hab.z[,H.keepers]
spc <- varpart(com.hel, d.C, d.H, d.S)
plot(spc)
rol

C.keepers

#setwd("/Users/melaniemcmillan/Desktop/McMillan_R/Spatial_Comm_Comp_F21-S22")
png("varpart.C&H&S.ex.spring.png", width = 150, height = 200)
plot(spc, Xnames = c("Water \nQuality","Habitat", "Space"))
title(main = "SPC (MI)", sub = "Drivers: average vegetative protection (L), latitude, PCNM2, \ncalcium, nitrate+nitrite")
dev.off()

#RDA with species
#Source: https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html
library(BiodiversityR) # also loads vegan
library(ggplot2)
library(readxl)
library(ggsci)
library(ggrepel)
library("ggforce")


com.hel <- disttransform(com, method='hellinger')
Ordination.model2 <- rda(com.hel ~ Stream, data=chem.z, scaling="species")
summary(Ordination.model2)

plot2 <- ordiplot(Ordination.model2, choices=c(1,2))
plot2

sites.long2 <- sites.long(plot2, env.data=chem.z)
head(sites.long2)

species.long2 <- species.long(plot2)
species.long2

axis.long2 <- axis.long(Ordination.model2, choices=c(1, 2))
axis.long2

spec.envfit <- envfit(plot2, env=chem.z)
spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(plot2, spec.data=spec.data.envfit)
species.long2

species.long3 <- species.long2[species.long2$r >= 0.6, ]
species.long3

plotgg2 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2, colour=Stream, shape=Stream), 
             size=5) +
  geom_segment(data=species.long3, 
               aes(x=0, y=0, xend=axis1*4, yend=axis2*4), 
               colour="red", size=0.7, arrow=arrow()) +
  geom_text_repel(data=species.long3, 
                  aes(x=axis1*4, y=axis2*4, label=labels),
                  colour="red") +
  theme_classic() +
  ggsci::scale_colour_npg() +
  coord_fixed(ratio=1)

plotgg2 

png("PCOA.fall.species.chem.png", width = 550, height = 550)
plot(plotgg2)
dev.off()

#RDA, CCA, CAP, dbRDA
#https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html
com <- read.csv("bugid.csv") %>%
  filter(Season == "Spring") %>%
  #filter(Stream == "CRO (R)") %>%
  dplyr::select(-c(Stream:Site))
com <- com[, colSums(com != 0, na.rm = TRUE) > 0]

chem <- read.csv("chem.f21-s22.notrib.reduced.csv")%>%
  filter(season == "Spring") %>%
  #filter(Stream == "ROL (MI)") %>%
  #filter(Stream == "CRO (R)") %>%
  dplyr::select(-c(Stream:season))
chem.z <- decostand(chem, method = "standardize",  na.rm =TRUE)
chem.z <- chem.z %>%
  select_if(~ ! any(is.na(.)))
round(apply(chem.z, 2, mean), 1)
apply(chem.z, 2, sd)
springrow <-  read.csv("chem.f21-s22.notrib.csv") %>%
  filter(season == "Spring") %>%
  dplyr::select(Stream) 
chem.z$Stream <- springrow$Stream

met <- read.csv("metrics.f21-s22.csv")%>%
  filter(Season == "Spring") %>%
  #filter(Stream == "ROL (MI)") %>%
  dplyr::select(-c(Stream:totind))
met.z <- decostand(met, method = "standardize",  na.rm =TRUE)
met.z <- met.z %>%
  select_if(~ ! any(is.na(.)))
round(apply(met.z, 2, mean), 1)
apply(met.z, 2, sd)
springrow <-  read.csv("metrics.f21-s22.csv") %>%
  filter(Season == "Spring") %>%
  dplyr::select(Stream) 
met.z$Stream <- springrow$Stream

hab <- read.csv("habitatmaster.csv") %>%
  select(-c(Stream:smallcobble)) 
hab.z <- decostand(hab, method = "standardize",  na.rm =TRUE)
hab.z <- hab.z %>%
  select_if(~ ! any(is.na(.)))
round(apply(hab.z, 2, mean), 1)
apply(hab.z, 2, sd)
hab.z$Stream <- springrow$Stream

all <- cbind(met, hab, chem) %>%
  select(pSC, pChiO , pP.less.Amph , 
           pEPT.less.HBL , pCling , pCG , pD , pPR , rich.D , 
           pSC.less.E , pCF , rich.INT , pSprawl , rich.P
         , pOligo, ca.mg , ba.ugl , mg.mgl , sc.uScm , 
           co.ugl , u.ugl , al.ugl , ni.ugl , cl.mgl , mn.ugl , temp.C , 
           no2no3.n.mgl, D10 , avgembedd , psmallcobble , 
           avgcancov , avgvegprotecR, avg.slope)
all.z <- decostand(all, method = "standardize",  na.rm =TRUE)
all.z <- all.z %>%
  select_if(~ ! any(is.na(.)))
round(apply(all.z, 2, mean), 1)
apply(all.z, 2, sd)

# add stream?

# hellinger transform the species dataset: gives low weights to rare species 
spe.hel <- decostand(com, "hellinger")

# Calculate distance matrix
bc<-vegdist(spe.hel, method="bray", binary=FALSE) 

# look at an unconstrained ordination first, it is always a good idea to look at both unconstrained and constrained ordinations
# set the seed: to reproduce the same result in the fture
set.seed(100)
bci.mds<-metaMDS(spe.hel, distance = "bray", k = 2)

# extract x and y coordinates from MDS plot into new dataframe, so you can plot with ggplot 
MDS_xy <- data.frame(bci.mds$points)
bci.mds$stress # 0.1853822

# colour by island
ggplot(MDS_xy, aes(MDS1, MDS2, col=met.z$Stream)) + geom_point() + theme_bw() + ggtitle('stress:0.185')

mod.0 <- rda(spe.hel ~ 1, data = met.z)
mod.1 <- rda(spe.hel ~ ., data = met.z)

#stepwise selection of the best model
simpleRDA <- ordiR2step(mod.0, scope = mod.1, R2scope = FALSE)
simpleRDA$call

simpleRDA <- rda(formula = spe.hel ~  pP.less.Amph + 
                    pD + pPR + 
                   pSC.less.E + pCF + pSprawl + ca.mg + mg.mgl + sc.uScm + 
                   no2no3.n.mgl + D10 + avgembedd + 
                   avgcancov + avg.slope , data = all.z)

simpleRDA <- rda(formula = spe.hel ~ D10 + avgembedd + psmallcobble + 
      avgcancov + avgvegprotecR, data = hab.z)# Stream group removed, Spring

simpleRDA <- rda(formula = spe.hel ~ Stream + ca.mg + ba.ugl + mg.mgl + sc.uScm + 
                   co.ugl + u.ugl + al.ugl + ni.ugl + cl.mgl + mn.ugl + temp.C + 
                   no2no3.n.mgl, data = chem.z)  # Stream group removed, Spring
  
simpleRDA <- rda(formula = spe.hel ~  pChiO + pP.less.Amph + pEPT.less.HBL + 
                   pCling + pCG + pD + pPR + rich.D + pSC.less.E + pCF + pEPT.less.H + 
                   rich.INT + pSprawl + pOligo, data = met.z)
  
summary(simpleRDA)
#rda(formula = spe.hel ~ ba.ugl + v.ugl + ti.ugl + na.mgl + so4.hco3 + 
#k.mgl + ca.mg + cl.mgl + sr.ugl + hco3.mgl + do.mgl + al.ugl + 
  #mn.ugl, data = chem.z) #Spring no stream
#rda(formula = spe.hel ~ Stream + ca.mg + ba.ugl + mg.mgl + sc.uScm + 
#co.ugl + u.ugl + al.ugl + ni.ugl + cl.mgl + mn.ugl + temp.C + 
  #no2no3.n.mgl, data = chem.z) #Spring with Stream

#Baetis=com$Baetis
#test <- cbind(Baetis, chem.z)
#mod <- lm(Baetis ~., data= test)
#summary(mod)

M.best <- names(simpleRDA$terminfo$ordered)
M.best
summary(simpleRDA)
screeplot(simpleRDA) #bstick not available for constrained ordinations

coef(simpleRDA)
R2 <- RsquareAdj(simpleRDA)$r.squared
R2 
R2adj <- RsquareAdj(simpleRDA)$adj.r.squared
R2adj 

#???lc???: orthogonal linear combinations of the explanatory variable (display=c("lc", "sp", "cn"))
# ???wa???: more robust to noise in the environmental variables but are a step between constrained towards unconstrained.(display="sp")
# Triplot: three different entities in the plot: sites, response variables and explanatory variables (arrowheads are on the explanatory variables)
# Scaling 1, wa = weighted sums of species
plot(simpleRDA, scaling=1, main="Triplot RDA matrix ~ env - scaling 1 - wa scores")
# arrows for species are missing, so lets add them without heads so they look different than the explanatory variables
spe.sc <- scores(simpleRDA, choices=1:2, scaling=1, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')

# Scaling 2
#wa
plot(simpleRDA, main="Triplot RDA matrix ~ env - scaling 2 - wa scores")
spe2.sc <- scores(simpleRDA, choices=1:2, display="sp") # scores() choices= indicates which axes are to be selected, make sure to specify the scaling if its different than 2 
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='red')

#IC
plot(simpleRDA, display=c("sp", "lc", "cn"), main="Triplot RDA matrix ~ env -scaling2-lc scores")
arrows(0,0,spe2.sc[,1],spe2.sc[,2], length=0, lty=1,col='red')

length()

library(ggord)
library(vcd)
plot <- ggord(simpleRDA, springrow$Stream, exp = 0.05, repel = TRUE, ext = 1.35, 
              size = 7, veclsz = 1, txt = 10, cols = c("chartreuse",
                                                       "chartreuse3",
                                                       "violetred1",
                                                       "violetred2",
                                                       "violetred3",
                                                       "violetred4")) + #facet = TRUE, xlims = -1, 1, ylims = -1,1
  theme(axis.title.x =element_text(size = 50), axis.text.x = element_text(size = 25),
        axis.title.y =element_text(size = 50), axis.text.y  = element_text(size = 25),
        strip.text.x = element_text(size = 30), 
        legend.text=element_text(size=20), legend.title=element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
        
  # looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  
plot

png("RDA.spring.allmetric.png", width = 1500, height = 1000)
plot
dev.off()
# CCA constrained by the same environmental data
simpleCCA <- cca(com ~ Stream + ca.mg + ba.ugl + mg.mgl + sc.uScm + 
                   co.ugl + u.ugl + al.ugl + ni.ugl + cl.mgl + mn.ugl + temp.C + 
                   no2no3.n.mgl, data = chem.z ) # Notice we are not using the hellinger transformation that downweights the importance of rare species

summary(simpleCCA)
# variation is now expressed as the mean squared contigency coefficient (biased and is not easily adjusted)
# species scores are represented as points
# site scores are averages of species scores 

screeplot(simpleCCA) # the first two axes are not as clear  

plot(simpleCCA, scaling=1, display=c('sp', 'lc', 'cn'), 
     main='Triplot CCA matrix ~ env -scaling 1')

# plot the CCA using ggplot (ggord package)
ggord(simpleCCA, chem.z$Stream) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # looking at the raw code, this is plotting the 'wa scores'

dbRDA <- capscale(com ~ Stream + ca.mg + ba.ugl + mg.mgl + sc.uScm + 
                    co.ugl + u.ugl + al.ugl + ni.ugl + cl.mgl + mn.ugl + temp.C + 
                    no2no3.n.mgl, data = chem.z , distance = "bray")

dbRDA
ggord(dbRDA, chem.z$Stream) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # looking at the raw code, this is plotting the 'wa scores'




