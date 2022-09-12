# PART 1) Simple regression of main metrics known to differ between streams by stream
# Part 2) Determine which habitat and water chemistry metrics are correlated to which insect metrics using NMDS. One NMDS graph with vectors of water quality and habitat on-top of raw insect community data.
# PART 3) Step-wise multiple regression of habitat and water quality to community similarity

#PART 1: metrics of interest: bugs) E.rich, E.rich.less.B, %E less B, EPT.rich, SC.rich, %E
# shannon, total richness, %5 Dominant Taxa, predator richness, scraper richness; 
# Water) SC, Se, SO4:HCO3; Ca:Mg ; Habitat) watershed size, stream gradient, LCF

install.packages("lmerTest")
library(lmerTest)
library(ggplot2)
library(ggpubr)

metrics <- read.csv("metrics.f21.csv")
chem <- read.csv("chem.f21-S22.csv")
chem <- filter( chem, season == "Fall") 
reg2 <- left_join(metrics, chem, by = c("Site" = "site.id"))

scxEPTrich <- ggplot(reg2,aes(x=sc.uScm,y=rich.E, group = Stream))+
  geom_point(aes(color=Stream, pch = Stream))+
  geom_line(aes(color=Stream, pch = Stream)) +
  stat_cor(method = "spearman",
           label.y = 8) +
  stat_regline_equation(label.y = 0.6) +
  ylim(0,8) +
  xlim(0, 2500) +
  theme_classic(base_size = 15, )+
  xlab("SC (uS/cm)") +
  ylab("Mayfly Richness") +
  ggtitle("Mayfly Richness vs Specific Conductance")+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  facet_wrap("Stream") +
  theme(legend.position = "none") +
  theme_classic()

scxEPTrich


#PART 2 
#install.packages("glmulti")
#library(glmulti)
library(dplyr)
library(car)

# prepare data to be used
sim <- read.csv("Similarity.f21.csv")
sim <- sim %>% filter(START %in% c("CRO-1", "EAS-1", "LLW-1", "ROL-1", "SPC-1", "FRY-2")) %>%
  filter(END != "FRY-1") %>%
  filter(END != "SPC-6") %>%
  filter(END != "FRY-8")
chem <- read.csv("chem.f21-S22.csv")
chem <- filter( chem, season == "Fall") 
multireg <- left_join(sim, chem, by = c("END" = "site.id")) 
multireg <- select(multireg, -c(1:5, 8:12, 28, 46:54))

which(is.na(multireg))

library(lme4)

#define intercept-only model
intercept_only <- lm( SIM ~ 1, data=multireg)
summary(intercept_only)

#define model with all predictors
colnames(multireg)
all <- lm(SIM ~ temp.C + hardness.mgl + no2no3.n.mgl + so4.mgl + k.mgl + ti.ugl +
            mn.ugl + zn.ugl + dis.simm + sc.uScm + tn.mgl + ortho.po4.p.mgl + hco3.mgl 
          + ca.mgl + v.ugl + co.ugl + as.ugl +do.mgl + ac.uScm + tp.mgl + npoc.mgl + 
            na.mgl + li.ugl + ph + alkalinity.mgl + nh3.n.mgl +
            cl.mgl + mg.mgl + al.ugl + se.ugl + fe.ugl + cu.ugl, data=multireg)
summary(all)
plot(all)

#perform forward stepwise regression
forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)
plot(forward)
#view results of forward stepwise regression
forward$anova
