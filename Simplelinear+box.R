# PART 1) Simple regression of main metrics known to differ between streams by stream
# Part 2) Determine which habitat and water chemistry metrics are correlated to which insect metrics using NMDS. One NMDS graph with vectors of water quality and habitat on-top of raw insect community data.
# PART 3) Step-wise multiple regression of habitat and water quality to community similarity

#PART 1: metrics of interest: bugs) E.rich, E.rich.less.B, %E less B, EPT.rich, SC.rich, %E
# shannon, total richness, %5 Dominant Taxa, predator richness, scraper richness; 
# Water) SC, Se, SO4:HCO3; Ca:Mg ; Habitat) watershed size, stream gradient, LCF

# Boxplots. Top 3 that change between my streams pEPT.less.HBL, pE.less.B, and pT
# most expected to change with SC rich.SC, rich.EPT, rich.E.less.B, p5dom, Hshannon

library(ggpubr)

metrics <- read.csv("metrics.f21-s22.csv")
#distance.d <- read.csv("distance.csv")
chem <- read.csv("chem.f21-S22.notrib.csv") 
#chem <- left_join(distance.d, chem, by = c("site" = "Site")) 
hab <- read.csv("habitatmaster.csv")
reg <- left_join(metrics, chem, by = c("Site", "Season" = "season")) 
reg.s <- filter(reg, Season == "Spring")
reg.f <- filter(reg, Season == "Fall")


library(rstatix)
library("tidyverse")
library("egg") #The egg package contains one of my favorite themes, theme_article.
library("multcompView")
library("sigminer")

letters.df <- data.frame(multcompLetters(TukeyHSD(aov(rich.SC ~ Stream.y, 
                                                      data = reg))$Stream.y[,4])$Letters)

colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Stream.y <- rownames(letters.df) #Create column based on rownames

placement <- reg %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream.y) %>%
  summarise(quantile(rich.SC)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) #Merge dataframes

box <- reg %>%
  ggplot(aes(x= Stream.y, y= rich.SC, group = Stream.y, color = Stream.y)) +
  geom_boxplot( alpha = 0) +
  facet_wrap("Season") + 
  geom_text(data = letters.df, aes(x = Stream.y, y = Placement.Value, label = Letter),
            size = 4, color = "black", hjust = -1.25, vjust = -0.8, fontface = "bold")+
  theme_classic() +
  ylab("Scraper Richness") +
  xlab("Stream") +
  ggtitle("Scraper Richness by Stream") +
  guides(color = guide_legend(title = "Stream")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  stat_compare_means(method = "kruskal", label.x.npc = "left", label.y.npc = "bottom") 
box

box <- reg %>%
  ggplot(aes(x= Stream.y, y= rich.SC, group = Stream.y, color = Stream.y)) +
  geom_boxplot( alpha = 0) +
  facet_wrap("Season") + 
  #geom_text(data = letters.df, aes(x = Stream.y, y = Placement.Value, label = Letter),
            #size = 4, color = "black", hjust = -1.25, vjust = -0.8, fontface = "bold")+
  theme_classic() +
  ylab("Scraper Richness") +
  xlab("Stream") +
  ggtitle("Scraper Richness by Stream") +
  guides(color = guide_legend(title = "Stream")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  stat_compare_means(method = "kruskal", label.x.npc = "left", label.y.npc = "bottom") +
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.")

box

pairs <- compare_means(rich.SC ~ Stream.y, data = reg, 
              group.by = "Season", method = "wilcox.test")

png("box.tukeyhsd+aov.richElessBxstreamxseason.png", width = 900, height = 450)
plot(box)
dev.off()

# Top insect vs sc metrics of interest: rich.SC, rich.EPT, rich.E.less.B, p5dom, Hshannon
# Top insect vs water quality metrics: sc.uScm, se.ugl, ca/mg, so4/hco3; Hshannon and rich.EPT
# Top insect vs habitat metrics: LCF, avg.slope, avgembedd, pfines; rich.PR and rich.SC

library (ggplot2)

plot <- ggplot(chem ,aes(x=dist.d, y=sc.uScm, group = season, color = season))+
  geom_point()+
  # geom_jitter() +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman") +
  stat_regline_equation(label.x.npc = 0.5, label.y.npc = 0.95) +
  #ylim(0,8) +
  #xlim(0, 2500) +
  theme_classic() +
  xlab("Distance Downstream") +
  ylab("Specific Conductance (uS/cm)") +
  ggtitle("Specific Conductance (uS/cm) vs Distance Downstream") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  labs(color="Season") +
  facet_wrap("Stream")  
plot


plot <- ggplot(reg ,aes(x=sc.uScm, y=rich.E.less.B, group = Season, color = Season))+
  geom_point()+
  # geom_jitter() +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman") +
  stat_regline_equation(label.x.npc = 0.5, label.y.npc = 0.95) +
  #ylim(0,8) +
  #xlim(0, 2500) +
  theme_classic() +
  xlab("Specific Conductance (uS/cm)") +
  ylab("Richness of Ephemeroptera less Baetidae") +
  ggtitle("Richness of Ephemeroptera less Baetidae vs Specific Conductance (uS/cm)") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  facet_wrap("Stream.y")  

plot


#set size and save plot as png
png("scxdistd.season.png", width = 900, height = 450)
plot(plot)
dev.off()


#PART 2 MULTIPLE REGRESSION??
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
