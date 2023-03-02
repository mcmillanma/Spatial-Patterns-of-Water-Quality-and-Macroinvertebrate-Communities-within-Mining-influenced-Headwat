#Redundancy Analysis
# https://r.qcbs.ca/workshop10/book-en/exploration.html
library(vegan)
library(ggplot2)
library(tidyverse)


com <- read.csv("bugid.csv")  %>%
  filter(Season == "Spring") %>%
  #filter(Stream == "CRO") %>%
  select(-c(Season, Site, Stream))
com <- com[, colSums(com != 0, na.rm = TRUE) > 0]

chem <- read.csv("chem.pca.csv")%>%
  filter(season == "Spring") %>%
#filter(Stream == "SPC (MI)") %>%
select(-c("Stream", Site:season))
chem <- chem[, colSums(chem != 0, na.rm = TRUE) > 0] 
springrow <-  read.csv("chem.pca.csv") %>%
  filter(season == "Spring") %>%
select(Stream) 

hab <- read.csv("hab.pca.csv") %>%
  #filter(Site != "EAS1") #%>%
  select(-c(Stream:Site))

#Examine Community Data

sum(com == 0)
# Calculate proportion of zeros in the dataset
sum(com == 0)/(nrow(com) * ncol(com))

# Apply Hellinger transformation to correct for the double
# zero problem 77% zeros, double zero problem
spe.hel <- decostand(com, method = "hellinger")

# Examine environmental data

# Scale and center variables
chem <- chem %>% 
  select(do.mgl:elements.pca) %>%
  select(-c(elements.pca))
  
chem.z <- decostand(chem, method = "standardize",  na.rm =TRUE)
chem.z <- chem.z %>%
  select_if(~ ! any(is.na(.)))

# check for colliniearity
library(car)

pairs(chem.z)
cor <- cor(chem.z)
#library (car)
model <- lm(sc.uScm ~ . , data = chem.z)
vif(model)

#alias <- alias( lm(sc.uScm ~ ., data = chem.z) )
#alias

# Variables are now centered around a mean of 0
round(apply(chem.z, 2, mean), 1)
apply(chem.z, 2, sd)

chem.z$Stream <- springrow$Stream

# Habitat
# Scale and center variables

hab.z <- decostand(hab, method = "standardize", na.rm =TRUE)
hab.z <- hab.z[, colSums(hab.z != 0, na.rm = TRUE) > 0]

# Variables are now centered around a mean of 0
round(apply(hab.z, 2, mean), 1)
apply(hab.z, 2, sd)

# check for colliniearity
pairs(hab.z)

model <- lm(ppebbles ~ ., data = hab.z)
vif(model)

alias <- alias( lm(avg.slope ~ ., data = hab.z) )
alias

cor <- cor(hab.z)
cor 

hab.z$Stream <- springrow$Stream


# VARIATION PARTITIONING
#### Sokol's Variation partitioning tutorial
library("geosphere")
library(sp)
library(sf)

spring.xy <- read.csv("spring_coord_notrib.csv") %>%
  filter(Stream == "EAS (R)") %>%
  filter(Site != "EAS9") %>%
  #filter(Site != "FRY8")%>%
  #filter(Site != "FRY9")%>%
  #filter(Site != "ROL7") %>%
  #filter(Site != "ROL7") %>%
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

plot <- ordisurf(spring.xy, scores(mod.pcnm, choi=2), bubble = 4, knots = 7)

png("PCNMall.png", width = 600, height = 600)
plot(plot, main = "SPC (MI) Spring 2022: PCNM 2")
dev.off()

com <- read.csv("bugid.csv")  %>%
  filter(Season == "Fall") %>%
  #filter(Site != "EAS1") %>%
  filter(Stream == "EAS (R)")%>%
  select(-c(Stream:Site))
com <- com[, colSums(com != 0, na.rm = TRUE) > 0] 
com.hel <- decostand( com, "hel")

# Space
d.space <- data.frame( spring.xy, vectors.pcnm)
#combine all spatial variables, x, y, and PCNMs

d.space.scaled <- data.frame( scale(d.space) ) %>%
  select(-c(x:y))
#center spatial variables on 0, and standardize
# null model with intercept
mod.0 <- rda( com.hel ~ 1, data = d.space.scaled)
plot(mod.0)

# model with all spatial variables included
mod.1 <- rda(com.hel ~ ., data = d.space.scaled)
summary(mod.1)

#stepwise selection of the best model
mod.best <- ordiR2step(rda( com.hel ~1, data = d.space.scaled), scope = formula(mod.1),  
                       direction = "backward",
                       R2scope = FALSE, # can't surpass the "full" model's R2
                       pstep = 1000,
                       trace = FALSE)
anova.cca(mod.best)
RsquareAdj(mod.best)

plot(mod.best)

S.keepers <- names( mod.best$terminfo$ordered )
S.keepers

# Chem
chem <- read.csv("chem.pca.csv")%>%
  filter(season == "Fall") %>%
  filter(Stream == "EAS (R)") %>%
  #filter(Site != "ROL7")%>%
  select(do.mgl:elements.pca)# %>%
#select(-c(elements.pca))

library(Hmisc)
model <- lm(sc.uScm ~ . , data = chem)
vif(model)
alias(model)
cro <- rcorr(as.matrix(chem), type = "spearman")
cro.r = data.frame(cro$P) %>%
  select("sc.uScm")

chem <- chem[, colSums(chem != 0, na.rm = TRUE) > 0] 
chem.z <- decostand(chem, method = "standardize",  na.rm =TRUE)
chem.z <- chem.z %>%
  select_if(~ ! any(is.na(.)))

# Variables are now centered around a mean of 0
round(apply(chem.z, 2, mean), 1)

apply(chem.z, 2, sd)

#chem.z$Stream <- springrow$Stream

mod.0 <-rda(com.hel ~ 1, data = chem.z)

mod.1 <- rda(com.hel ~ ., data = chem.z)
summary(mod.1)

#stepwise selection of the best model
mod.best <- ordiR2step(mod.0, scope = mod.1, direction = "backward",R2scope = FALSE)
anova.cca(mod.best)
RsquareAdj(mod.best)

C.keepers <- names(mod.best$terminfo$ordered)
C.keepers

#Habitat
hab <- read.csv("hab.pca.csv") %>%
filter(Stream == "EAS (R)") %>%
  filter(Site != "EAS9") %>%
  #filter(Site != "FRY8")%>%
  #ilter(Site != "FRY9")%>%
  #filter(Site != "ROL7") %>%
  #filter(Site != "SPC6") %>%
  select(-c(Stream:Site)) 
hab.z <- decostand(hab, method = "standardize", na.rm =TRUE)
hab.z <- hab.z[, colSums(hab.z != 0, na.rm = TRUE) > 0]

# Variables are now centered around a mean of 0
round(apply(hab.z, 2, mean), 1)
apply(hab.z, 2, sd)

#hab.z$Stream <- springrow$Stream

mod.0 <- rda(com.hel ~ 1, data = hab.z)
mod.1 <- rda(com.hel ~ ., data = hab.z)

#stepwise selection of the best model
mod.best <- ordiR2step(mod.0, scope = mod.1, direction = "backward",
                       R2scope = FALSE, # can't surpass the "full" model's R2
                       pstep = 1000,
                       trace = FALSE) 

anova.cca(mod.best)
RsquareAdj(mod.best)

H.keepers <- names(mod.best$terminfo$ordered)
H.keepers

mod.best

# put them all together
d.C <- chem.z[,C.keepers]
d.S <- d.space.scaled[,S.keepers]
d.H <- hab.z[,H.keepers]
all <- varpart(com.hel, d.C,d.H, d.S)
all
plot(all)

# [a+b] Chemistry without controlling for topography

S.keepers
C.keepers
H.keepers


#setwd("/Users/melaniemcmillan/Desktop/McMillan_R/Spatial_Comm_Comp_F21-S22")
png("varpart.C&H&S.springcro.png", width = 500, height = 500)
plot(all, Xnames = c("Water \nQuality","Habitat", "Space"))
title(main = "ROL (MI) Spring", sub = "Drivers: PCNM1&2, elements, and slope ")
dev.off()

###RDA with water chem and habitat data
library(vegan) #decostand

chem <- read.csv("chem.pca.csv") %>%
  filter(season == "Spring")
hab <- read.csv("hab.pca.csv")
env.s <- left_join(chem, hab, by = c("Stream", "Site"))  
env.f <- left_join(chem, hab, by = c("Stream", "Site")) 

env <- filter(env.s, Stream =="SPC (MI)") 
env <- select(env, -c(Stream:season))
env <- env[, colSums(env != 0, na.rm = TRUE) > 0]
env <- decostand(env, method = "standardize", na.rm =TRUE)

com <- read.csv("bugid.csv")  %>%
  filter(Season == "Spring") %>%
  filter(Stream == "SPC (MI)") %>%
  select(-c(Season, Site, Stream))
com <- com[, colSums(com != 0, na.rm = TRUE) > 0]

# hellinger transform the species dataset: gives low weights to rare species 
spe.hel <- decostand(com, "hellinger")

# Calculate distance matrix
bc<-vegdist(spe.hel, method="bray", binary=FALSE) 

# look at an unconstrained ordination first, it is always a good idea to look at both unconstrained and constrained ordinations
# set the seed: to reproduce the same result in the fture
set.seed(100)
bci.mds<-metaMDS(spe.hel, distance = "bray")
bci.mds

# extract x and y coordinates from MDS plot into new dataframe, so you can plot with ggplot 
MDS_xy <- data.frame(bci.mds$points)
bci.mds$stress # cro 0.0588, eas 0.0495, fry 0.0232, llw 0, rol 0.041, spc 0.0501

# colour by island
ggplot(MDS_xy, aes(MDS1, MDS2)) + geom_point() + theme_bw() + ggtitle('stress:0.0495')


# Spring ALl
simpleRDA <- rda(formula = spe.hel ~ sc.uScm + hardness.mgl+ npoc.mgl+ temp.C + 
                   do.mgl + nutrients.pca + pfines + avg.slope +D.pca + 
                   avgembedd + avgwetwidth+ Rip.pca, data = env.s)
# Fall All
simpleRDA <- rda(formula = spe.hel ~ sc.uScm + hardness.mgl+ npoc.mgl+ elements.pca + 
                   do.mgl  + pfines + avg.slope +D.pca + 
                  avgwetwidth+ Rip.pca, data = env.f)
#Spring CRO
simpleRDA <- rda(formula = spe.hel ~ nutrients.pca + hardness.mgl + plargecobble 
                 , data = env)
#Spring EAS
simpleRDA <- rda(formula = spe.hel ~ nutrients.pca + sc.uScm +avg.slope + Rip.pca, data = env)

#Spring FRY
simpleRDA <- rda(formula = spe.hel ~ elements.pca + avg.slope, data = env)

# Spring LLW
simpleRDA <- rda(formula = spe.hel ~  hardness.mgl , data = env)

# Spring ROL
simpleRDA <- rda(formula = spe.hel ~ npoc.mgl + LCF, data = env)

#Spring SPC
simpleRDA <- rda(formula = spe.hel ~ nutrients.pca + avgwetwidth, data = env)

summary(simpleRDA)


coef(simpleRDA)
R2 <- RsquareAdj(simpleRDA)$r.squared
R2 
R2adj <- RsquareAdj(simpleRDA)$adj.r.squared
R2adj 

#‘lc’: orthogonal linear combinations of the explanatory variable (display=c("lc", "sp", "cn"))
# ‘wa’: more robust to noise in the environmental variables but are a step between constrained towards unconstrained.(display="sp")
# Triplot: three different entities in the plot: sites, response variables and explanatory variables (arrowheads are on the explanatory variables)
# Scaling 1, wa = weighted sums of species
plot(simpleRDA, scaling=1, main="SPC (MI), Spring - scaling 1 - wa scores")


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
plot <- ggord(simpleRDA, env.s$Stream, exp = 0.07, repel = TRUE, ext = 1.25, 
             veclsz = 0.7, ptslab = TRUE, addsize = 4, cols = c("chartreuse",
                                                                     "chartreuse3",
                                                                     "violetred1",
                                                                     "violetred2",
                                                                     "violetred3",
                                                                     "violetred4")) + #facet = TRUE, xlims = -1, 1, ylims = -1,1
  labs(title = "All Streams Fall") +
  theme_classic() +
  annotate("text", x = -0.5, y = 0.8, label = "R2adj = 0.37") +
  annotate("text", x = -0.5, y = 0.7, label = "stress = 0.19") 

bci.mds$stress
R2adj
plot

# looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  
cro
eas
fry
llw
rol
spc

plot <- ggarrange(eas, fry, rol, spc)
plot

png("RDA.all.f.png", width = 600, height = 430)
plot(plot)
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
Ordination.model2 <- rda(formula = com.hel ~ pSH + pSC + pBurrow + pChiO + rich.CF + 
                           Hshannon + pCF + pOligo + pSprawl + pEPT + rich.D + pSwimm + 
                           pCG + pP + rich.SH + pTOL + pT + pT.less.H + rich.Swimm + 
                           pCling + p5dom, data = met.z)
summary(Ordination.model2)

plot2 <- ordiplot(Ordination.model2, choices=c(1,2))
plot2

sites.long2 <- sites.long(plot2, env.data=met.z)
head(sites.long2)

species.long2 <- species.long(plot2)
species.long2

axis.long2 <- axis.long(Ordination.model2, choices=c(1, 2))
axis.long2

spec.envfit <- envfit(plot2, env=met.z)
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
             aes(x=axis1, y=axis2), 
             size=5) +
  geom_segment(data=species.long3, 
               aes(x=0, y=0, xend=axis1*4, yend=axis2*4), 
               colour="red", size=0.7, arrow=arrow()) +
  geom_text_repel(data=species.long3, 
                  aes(x=axis1*4, y=axis2*4, label=labels),
                  colour="red") +
  theme_classic() +
  ggsci::scale_colour_npg() 

plotgg2 

png("PCOA.fall.species.chem.png", width = 550, height = 550)
plot(plotgg2)
dev.off()

#RDA, CCA, CAP, dbRDA
#https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html
com <- read.csv("bugid.csv") %>%
  #filter(Season == "Spring") %>%
  #filter(Stream == "CRO (R)") %>%
  dplyr::select(-c(Stream:Site))
com <- com[, colSums(com != 0, na.rm = TRUE) > 0]

chem <- read.csv("chem.f21-s22.notrib.reduced.csv")%>%
  #filter(season == "Spring") %>%
  #filter(Stream == "ROL (MI)") %>%
  #filter(Stream == "CRO (R)") %>%
  dplyr::select(-c(Stream:season))
chem.z <- decostand(chem, method = "standardize",  na.rm =TRUE)
chem.z <- chem.z %>%
  select_if(~ ! any(is.na(.)))
shapiro.test(chem.z$sc.uScm)
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
shapiro.test(met.z$pE)
round(apply(met.z, 2, mean), 1)
apply(met.z, 2, sd)
springrow <-  read.csv("metrics.f21-s22.csv") %>%
  filter(Season == "Spring") %>%
  dplyr::select(Stream) 
met.z$Stream <- springrow$Stream

hab <- read.csv("habitatmaster.csv") %>%
  select(-c(Stream:smallcobble)) 
hab.z <- decostand(hab, method = "log",  na.rm =TRUE)
hab.z <- hab.z %>%
  select_if(~ ! any(is.na(.)))
shapiro.test(hab.z$pfines)
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
shapiro.test(spe.hel$Allocapnia)

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

simpleRDA <- rda(formula = spe.hel ~ Stream + pSC + pChiO + pP.less.Amph + 
                   pEPT.less.HBL + pCling + pEPT + pCG + pD + pPR + rich.D + 
                   pSC.less.E + pCF + pEPT.less.H + rich.INT + pSprawl + rich.P + 
                   rich.PR + pOligo, data = met.z)

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

#‘lc’: orthogonal linear combinations of the explanatory variable (display=c("lc", "sp", "cn"))
# ‘wa’: more robust to noise in the environmental variables but are a step between constrained towards unconstrained.(display="sp")
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




