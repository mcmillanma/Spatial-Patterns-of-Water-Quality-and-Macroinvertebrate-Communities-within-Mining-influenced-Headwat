# AGU Poster Graphs
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)

# Fig 3 & 4

Similar <- read_csv("Similarity.csv") %>%
  filter(Season == "Spring")

library(scales)
percents <- Similar %>%
  group_by(Stream) %>%
  summarise(maxb = max(SIM.b), minb = min(SIM.b), maxsc = max(sc.delt)) %>%
  mutate(changeb = 1-minb)
percents$changeb <- percent(percents$changeb, accuracy=1)
  

ex <- Similar %>% ggplot(aes(x=dis.simm , y=SIM.b, color = Stream))+
  geom_point()  +
  facet_wrap("Stream") +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman", label.x = 0, label.y = 1, size = 10) +
  #stat_regline_equation(label.x = 0, label.y = 0.93, size =10 ) +
  scale_color_manual(values = c("chartreuse",
                                 "chartreuse3",
                                 "violetred1",
                                 "violetred2",
                                 "violetred3",
                                 "violetred4")) +
  theme_classic()+
  xlab("Distance Between Samples (m)") + # (uS/cm)
  ylab("Bray-Curtis Similarity") +
  #labs(title = "Bray-Curtis Community Similarity vs \n Change in Specific Conductance") +
  theme(plot.title = element_text(hjust = 0.5, size = 30, face = "bold"), 
        axis.title.x =element_text(size = 50), axis.text.x = element_text(size = 25),
        axis.title.y =element_text(size = 50), axis.text.y  = element_text(size = 25),
        strip.text.x = element_text(size = 30), 
        legend.text=element_text(size=20), legend.title=element_text(size=25)) +
  geom_text(data    = percents, mapping = aes(x = -Inf, y = -Inf, label = changeb),
    hjust   = -0.1,vjust   = -1, size = 10)
ex
#set size and save plot as png
png("ex.png", width = 1500, height = 1000)
plot(ex)
dev.off()

#Figure 5

metrics <- read_csv("metrics.f21-s22.csv")
chem <- read_csv("chem.f21-s22.notrib.csv") %>%
  filter(season == "Spring")

sc <- left_join(metrics, chem, by = c("Site","Season" = "season", "Stream", "dist.d")) %>%
  filter(Season == "Spring")

#pE, rich.EPT, rich.SC

rich.SC<- sc %>% ggplot(aes(x=sc.uScm , y=rich.SC, color = Stream))+
  geom_point()  +
  facet_wrap("Stream") +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman", label.y = 15, size = 10) +
  scale_color_manual(values = c("chartreuse",
                                "chartreuse3",
                                "violetred1",
                                "violetred2",
                                "violetred3",
                                "violetred4")) +
  theme_classic()+
  xlab("SC (uS/cm)") +
  ylab("Scraper Richness") +
  #labs(title = "Percent Ephemeroptera vs \n Specific Conductance") +
  theme(axis.title.x =element_text(size = 50), axis.text.x = element_text(size = 25),
        axis.title.y =element_text(size = 50), axis.text.y  = element_text(size = 25),
        strip.text.x = element_text(size = 30), 
        legend.text=element_text(size=20), legend.title=element_text(size=25)) 
  
pE                
rich.EPT
rich.SC
pCF
pSprawl
pP.less.Amph

fig4 <- ggarrange(pE, rich.EPT, rich.SC, nrow = 1, ncol= 3, common.legend = TRUE, 
                  legend = "right")

fig4 <-annotate_figure(fig4, top = text_grob("Select Compositional Metrics vs Specific Conductance by Stream", 
                              face = "bold", size = 60))
fig4

#set size and save plot as png
png("fig4.png", width = 3000, height = 1000)
plot(fig4)
dev.off()

# Fig 6

rich.EPT.g <- sc %>% ggplot(aes(x=sc.uScm , y=rich.EPT))+
  geom_point(aes(color = Stream))  +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman", label.x = 0, label.y = 9, size = 10) +
  stat_regline_equation(label.x = 0, label.y = 5, size = 10) +
  theme_classic()+
  xlab("Specific Conductance (uS/cm)") +
  ylab("Ephemeroptera, Plecoptera,\nTrichoptera Richness") +
  #labs(title = "Percent Ephemeroptera vs \n Specific Conductance") +
  #ylim(0,60) +
  theme(axis.title.x =element_text(size = 50), axis.text.x = element_text(size = 25),
        axis.title.y =element_text(size = 50), axis.text.y  = element_text(size = 25),
        legend.text=element_text(size=20), legend.title=element_text(size=25))
  

pE.g                
rich.EPT.g
rich.SC.g

fig5 <- ggarrange(pE.g, rich.EPT.g, rich.SC.g, nrow = 1, common.legend = TRUE, 
                  legend = "right")

fig5 <-annotate_figure(fig5, top = text_grob("Select Compositional Metrics vs Specific Conductance Globally", 
                                             face = "bold", size = 60))
fig5

#set size and save plot as png
png("fig5.png", width = 3000, height = 1000)
plot(fig5)
dev.off()

chem %>%
  summarise(max_sc = max(sc.uScm), min_sc = min(sc.uScm)) %>%
  mutate(range = max_sc - min_sc)

chem %>%
  group_by(Stream) %>%
  summarise(max_dist = max(dist.d))

#https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html
com <- read.csv("bugid.csv") %>%
  filter(Season == "Fall") %>%
  filter(Stream == "CRO (R)") %>%
  dplyr::select(-c(Stream:Site))
com <- com[, colSums(com != 0, na.rm = TRUE) > 0]

chem <- read.csv("chem.f21-s22.notrib.reduced.csv")%>%
  filter(season == "Spring") %>%
  filter(Stream == "CRO (R)") %>%
  dplyr::select(-c(Stream:season))

met <- read.csv("metrics.f21-s22.csv") %>%
  filter(Season == "Fall") %>%
  filter(Stream == "CRO (R)") %>%
  dplyr::select(-c(Stream:totind))
met.z <- decostand(met, method = "standardize",  na.rm =TRUE)
met.z <- met.z %>%
  select_if(~ ! any(is.na(.)))
round(apply(met.z, 2, mean), 1)
apply(met.z, 2, sd)

hab <- read.csv("habitatmaster.csv") %>%
  filter(Stream == "CRO (R)") %>%
  select(-c(Stream:smallcobble)) 

all <- cbind(chem, hab)
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
bci.mds

# extract x and y coordinates from MDS plot into new dataframe, so you can plot with ggplot 
MDS_xy <- data.frame(bci.mds$points)
bci.mds$stress # cro 0.0588, eas 0.0495, fry 0.0232, llw 0, rol 0.041, spc 0.0501

# colour by island
ggplot(MDS_xy, aes(MDS1, MDS2)) + geom_point() + theme_bw() + ggtitle('stress:0.0495')

mod.0 <- rda(spe.hel ~ 1, data = met.z)
mod.1 <- rda(spe.hel ~ ., data = met.z)

#stepwise selection of the best model
simpleRDA <- ordiR2step(mod.0, scope = mod.1, R2scope = FALSE)
simpleRDA$call

simpleRDA <- rda(formula = spe.hel ~ pP + pChiO, data = met.z)

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
spc <- ggord(simpleRDA, exp = 0.05, repel = TRUE, ext = 1.35, 
              size = 7, veclsz = 1, txt = 10, obslab = TRUE, cols = c("chartreuse",
                                                       "chartreuse3",
                                                       "violetred1",
                                                       "violetred2",
                                                       "violetred3",
                                                       "violetred4")) + #facet = TRUE, xlims = -1, 1, ylims = -1,1
  labs(title = "SPC (MI)", size = 50) +
  theme(axis.title.x =element_text(size = 50), axis.text.x = element_text(size = 25),
        axis.title.y =element_text(size = 50), axis.text.y  = element_text(size = 25),
        strip.text.x = element_text(size = 30), 
        legend.text=element_text(size=20), legend.title=element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_text(size = 50), ) 
  


# looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  
cro
eas
fry
llw
rol
spc

plot <- ggarrange(eas, fry, rol, spc)
plot

png("RDA.fall.spc.png", width = 600, height = 430)
plot(spc)
dev.off()
  

  
