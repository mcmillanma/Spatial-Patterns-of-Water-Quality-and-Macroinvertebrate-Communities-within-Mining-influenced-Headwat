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
chem <- read_csv("chem.f21-s22.notrib.reduced.csv") %>%
  filter(season == "Spring")

sc <- left_join(metrics, chem, by = c("Site","Season" = "season", "Stream", "dist.d")) %>%
  filter(Season == "Spring")

#pE, rich.EPT, rich.SC

VASCI<- sc %>% ggplot(aes(x=dist.d , y=VASCI, color = Stream))+
  geom_point()  +
  facet_wrap("Stream") +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman", label.y = 50) +
  scale_color_manual(values = c("chartreuse",
                                "chartreuse3",
                                "violetred1",
                                "violetred2",
                                "violetred3",
                                "violetred4")) +
  theme_classic()+
  xlab("Distance Downstream (m)") +
  ylab("VASCI") 
VASCI 
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
png("VASCIreg.png", width = 750, height = 600)
plot(VASCI)
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

mod.0 <- rda(com.hel ~ 1, data = met.z)
mod.1 <- rda(com.hel ~ ., data = met.z)

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

######################################################
######   VLWA   ############################################################
######################################################

chem <- read_csv("chem.pca.csv") %>%
  filter(season == "Spring")
hab <- read.csv("hab.pca.csv")
env <- left_join(chem, hab, by = c("Site", "Stream"))
env <- left_join(env, metrics, by = c("Site", "Stream",  "season" = "Season")) %>%
  select(-c(27:66))
env$group <- factor(env$stream, levels = c("Reference 1", "Reference 2", "Mined 1", "Mined 2", "Mined 3", "Mined 4"))
boxplot(env$VASCI ~ env$group)

metrics <- read_csv("metrics.reduced.csv") %>%
  filter(Season == "Spring")
all <- left_join(chem, metrics, by = c("stream", "Stream", "Impact", "Site"))
all$group <- factor(all$stream, levels = c("Reference 1", "Reference 2", "Mined 1", "Mined 2", "Mined 3", "Mined 4"))
boxplot(all$sc.uScm ~ all$group)

  
longterm <- env %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))

box <- all %>% #Dataframe from which data will be drawn
  ggplot(aes(x = group, y = sc.uScm, color = Impact)) +#Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0, size = 1) + #I like to set the color of boxplots to black with the alpha at 0 (fully transparent). I also like geom_jitter() but do not use it here for simplicity.
  geom_point(size = 3) +
  geom_point(data = longterm, aes(stream, sc.uScm), color = "black") +
  theme_classic() +
  #scale_colour_viridis_d()+ #viridis color blind friendly
  xlab("Stream (Spring)") +
  ylab("SC (uS/cm)") +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"))#+
#geom_text(data = letters.df, aes(x = Stream, y = Placement.Value, label = Letter), size = 4, color = "black", hjust = -1.25, vjust = -0.8, fontface = "bold")
box

png("SC.box.png", width = 850, height = 600)
plot(box)
dev.off()

Predper <- all %>% #Dataframe from which data will be drawn
  ggplot(aes(x =group, y = pPR, color = Impact)) +#Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0, size = 1) + #I like to set the color of boxplots to black with the alpha at 0 (fully transparent). I also like geom_jitter() but do not use it here for simplicity.
  geom_point(size = 3) +
  geom_point(data = longterm, aes(group,pPR), color = "black") +
  theme_classic() +
  #scale_colour_viridis_d()+ #viridis color blind friendly
  xlab("Stream (Spring)") +
  ylab("Percent Predators") +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"))#+
#geom_text(data = letters.df, aes(x = Stream, y = Placement.Value, label = Letter), size = 4, color = "black", hjust = -1.25, vjust = -0.8, fontface = "bold")
p5dom
EPTrich
Shan
Scraprich
Clingrich
Predper

plot <- ggarrange(p5dom, EPTrich, Shan, Scraprich, Clingrich, Predper)
plot

png("SC.box.png", width = 2550, height = 1200)
plot(plot)
dev.off()

plot <- all %>% ggplot(aes(x=sc.uScm , y=VASCI))+
  geom_point(aes(color = Impact))  +
  stat_smooth(method = "lm", linetype = 2, colour = "black") +
  geom_point(data = longterm, aes(x=sc.uScm , y=VASCI), color = "black") +
  facet_wrap("group") +
  theme_classic()+
  xlab("Specific Conductance (uS/cm)") +
  ylab("VASCI") +
  scale_y_continuous() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))
plot

png("SC.line.png", width = 850, height = 600)
plot(plot)
dev.off()

plot <- env %>% ggplot(aes(x=avgwetwidth , y=VASCI))+
  geom_point(aes(color = Impact))  +
  stat_smooth(method = "lm", linetype = 2, colour = "black") +
  geom_point(data = longterm, aes(x=avgwetwidth , y=VASCI), color = "black") +
  facet_wrap("group") +
  theme_classic()+
  xlab("Wetted Width") +
  ylab("VASCI") +
  scale_y_continuous() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))
plot

png("slope.png", width = 850, height = 600)
plot(plot)
dev.off()


  


