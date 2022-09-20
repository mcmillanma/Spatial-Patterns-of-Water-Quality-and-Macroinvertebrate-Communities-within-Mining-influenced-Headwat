# PART 1) Simple regression of main metrics known to differ between streams by stream
# Part 2) Determine which habitat and water chemistry metrics are correlated to which insect metrics using NMDS. One NMDS graph with vectors of water quality and habitat on-top of raw insect community data.
# PART 3) Step-wise multiple regression of habitat and water quality to community similarity

#PART 1: metrics of interest: bugs) E.rich, E.rich.less.B, %E less B, EPT.rich, SC.rich, %E
# shannon, total richness, %5 Dominant Taxa, predator richness, scraper richness; 
# Water) SC, Se, SO4:HCO3; Ca:Mg ; Habitat) watershed size, stream gradient, LCF

#install.packages("lmerTest")
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(dplyr)

metrics <- read.csv("metrics.f21.csv")
chem <- read.csv("chem.f21-S22.csv")
chem <- chem %>%
  mutate("ratiocamg" = ca.mgl / mg.mgl) %>%
  mutate ("ratioso4hco3" = so4.mgl /hco3.mgl)
chem <- dplyr::filter(chem, season == "Fall") 
hab <- read.csv("habitatmaster.csv")
reg2 <- left_join(metrics, chem, by = "Site") 
reg2 <- left_join(reg2, hab, by = "Site")

# Top insect metrics: rich.SC, rich.EPT, rich.E.less.B, p5dom, Hshannon
# Top water quality metrics: sc.uScm, se.ugl, ca/mg, so4/hco3
# Top habitat metrics: LCF, avg.slope, avgembedd, pfines

EPTrichxse <- ggplot(reg2,aes(x=se.ugl, y=rich.EPT, group = Stream))+
  geom_point(aes(color=Stream, pch = Stream))+
  # geom_jitter() +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman",
           label.y = 25) +
  stat_regline_equation(label.y = 24) + #label.y = 
  #ylim(0,8) +
  #xlim(0, 2500) +
  theme_classic(base_size = 15, )+
  xlab("Se (ug/L)") +
  ylab("EPT richness") +
  ggtitle("EPT richness vs Selenium")+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  facet_wrap("Stream") +
  theme(legend.position = "none") +
  theme_classic()
EPTrichxse


#set size and save plot as png
png("EPTrichxse.png", width = 900, height = 450)
plot(EPTrichxse)
dev.off()

#create correlation matrix


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
