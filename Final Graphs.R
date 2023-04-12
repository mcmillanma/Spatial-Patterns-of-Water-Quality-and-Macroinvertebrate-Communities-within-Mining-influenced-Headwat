# Final Figures

### NMDS ###
library(BiodiversityR) # also loads vegan
library(ggplot2)
library(goeveg)
library(ggord)
library(tidyverse)
library(ggrepel)
library(FSA) #dunnTest


bugs <- read_csv("bugid.csv") %>%
  filter(Season == "Fall")
bugs <- bugs[, colSums(bugs != 0, na.rm = TRUE) > 0]
bugs$group <- factor(bugs$Stream, levels = c("EAS (R)", "CRO (R)", 
                                             "FRY (L)", "SPC (L)", 
                                             "ROL (H)", "LLW (H)"))

data_1 <- bugs[,4:97] # Spring 93, Fall 97

nmds <- metaMDS(data_1, distance = "bray", k=2, autotransform=TRUE)

limited <- ordiselect(data_1, nmds, ablim = 0.15) # ablim = proportion of species with highest abundances
limited

scoresspecies <- as.data.frame(scores(nmds)$species) 
scoresspecies <- subset(scoresspecies, rownames(scoresspecies) %in% limited)
scoresspecies$taxa <- row.names(scoresspecies)

datascores <- as.data.frame(scores(nmds)$sites)
data_2 <- bugs[,98:99]
scores <- cbind(datascores, Stream = data_2$group)
scores <- cbind(scores, Impact = data_2$Impact)

stress <- nmds$stress
stress

#ggord(nmds, data_2$group, txt = limited) +
  #theme_classic()

plot <- ggplot(scores, aes(x= NMDS1, y = NMDS2)) +
         geom_point(aes(color = Impact, shape = Stream, size = 4)) +
  stat_ellipse(aes(data = Impact, color = Impact)) +
  geom_segment(data=scoresspecies, 
             aes(x=0, y=0, xend=NMDS1*1, yend=NMDS2*1), 
             colour="black", size=0.7, arrow=arrow()) +
  geom_text_repel(data=scoresspecies, 
                  aes(x=NMDS1*1.5, y=NMDS2*1.5, label=taxa),
                  colour="black", size = 7) +
  theme_classic() + 
  ggtitle("Fall 2021") +
  theme(plot.title = element_text(size = 20), axis.title.x =element_text(size = 20), axis.text.x = element_text(size = 18),
        axis.title.y =element_text(size = 20), axis.text.y  = element_text(size = 18),
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  annotate("text", x = 0.5, y = -0.7, label = "Stress = 0.22", size = 5)
plot

png("SpeciesNMDSALLS.png", width = 650, height = 550)
plot(plot)
dev.off()


###### NMDS by Stream ########

bugs <- read_csv("bugid.csv") %>%
  filter(Season == "Spring") %>%
  filter(Stream == "CRO (R)")
bugs <- bugs[, colSums(bugs != 0, na.rm = TRUE) > 0]

data_1 <- bugs[,4:67] # Fall: EAS 58, CRO 68, FRY 49, SPC 58, ROL and LLW 50
# Spring: LLW ROL 39, FRY

nmds <- metaMDS(data_1, distance = "bray", k=2, autotransform=TRUE)

limited <- ordiselect(data_1, nmds, ablim = 0.15) # ablim = proportion of species with highest abundances
limited

scoresspecies <- as.data.frame(scores(nmds)$species) 
scoresspecies <- subset(scoresspecies, rownames(scoresspecies) %in% limited)
scoresspecies$taxa <- row.names(scoresspecies)

datascores <- as.data.frame(scores(nmds)$sites)
data_2 <- bugs[,70] 
scores <- cbind(datascores, Site = data_2$Site.ord)

stress <- nmds$stress
stress

plot <- ggplot(scores, aes(x= NMDS1, y = NMDS2, label = Site)) +
  geom_point(color = "blue", size = 4) +
  geom_text(nudge_y = 0.06, size =8, color = "blue") + #chartreuse3, blue, chocolate2 
  geom_segment(data=scoresspecies, 
               aes(x=0, y=0, xend=NMDS1*2.5, yend=NMDS2*2.5), 
               colour="black", arrow=arrow(),inherit.aes = FALSE) +
  geom_text_repel(data=scoresspecies, 
                  aes(x=NMDS1*3 ,y=NMDS2*3, label=taxa),
                  colour="black", size = 7) +
  theme_classic() + 
  ggtitle("CRO (R) Spring 2022") +
  theme(plot.title = element_text(size = 20), axis.title.x =element_text(size = 20), axis.text.x = element_text(size = 18),
        axis.title.y =element_text(size = 20), axis.text.y  = element_text(size = 18),
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  annotate("text", x = -0.8, y = -0.25, label = "Stress = 0.05", size = 5)
plot

png("SpeciesNMDSALLS.png", width = 500, height = 460)
plot(plot)
dev.off()

################################
###### Means Comparisons #######
################################

library("multcompView")
library(ggpubr) #stat_compare_means

metrics <- read.csv("metrics.f21-s22.csv") 

metrics.s <- filter(metrics, Season == "Spring")
metrics.s$group <- factor(metrics.s$Stream, levels = c("EAS (R)", "CRO (R)", 
                                                       "FRY (L)", "SPC (L)", 
                                                       "ROL (H)", "LLW (H)"))
metrics.f <- filter(metrics, Season == "Fall")
metrics.f$group <- factor(metrics.f$Stream, levels = c("EAS (R)", "CRO (R)", 
                                                       "FRY (L)", "SPC (L)", 
                                                       "ROL (H)", "LLW (H)"))

longterm.f <- metrics.f %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))

longterm.s <- metrics.s %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))

# pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom

#
# Spring ANOVA
## Spring Anova: rich.SC, p5dom, rich.Cling; no transformations

set.seed(1045)
model <- aov(p5dom~Stream, data=metrics.s)

letters.df <- data.frame(multcompLetters(TukeyHSD(
  aov(p5dom~Stream, data=metrics.s))$Stream[,4])$Letters)

colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Stream <- rownames(letters.df) #Create column based on rownames

placement <- metrics.s %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(p5dom)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.s <- left_join(letters.df, placement, by = "Stream") #Merge dataframes

plot <- metrics.s %>% 
  ggplot(aes(x = group, y = p5dom, color = Impact)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) +
  geom_point(size = 4) +
  geom_point(data = longterm.s, aes(Stream, p5dom), color = "black") + 
  theme_classic() + 
  xlab("Stream (Spring)") +
  ylab("% 5 Dominant Taxa") +
  stat_compare_means(method = "anova", size = 6) +
  geom_text(data = letters.df.s, 
            aes(x = Stream, y = Placement.Value, label = Letter), 
            color = "black", hjust = -0.5,
            vjust = -0.8, fontface = "bold", size = 6) +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))

plot

png("box.png",width = 650, height = 550)
plot(plot)
dev.off()

###
##### Fall Anova
#### Fall Anova: rich.SC, p5dom; no transformations

set.seed(1045)
model <- aov(p5dom~Stream, data=metrics.f)

letters.df <- data.frame(multcompLetters(TukeyHSD(
  aov(p5dom~Stream, data=metrics.f))$Stream[,4])$Letters)

colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Stream <- rownames(letters.df) #Create column based on rownames

placement <- metrics.f %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(p5dom)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.s <- left_join(letters.df, placement, by = "Stream") 

plot <- metrics.f %>% 
  ggplot(aes(x = group, y = p5dom, color = Impact)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) +
  geom_point(size = 4) +
  geom_point(data = longterm.f, aes(Stream, p5dom), color = "black") + 
  theme_classic() + 
  xlab("Stream (Fall)") +
  ylab("% 5 Dominant Taxa") +
  stat_compare_means(method = "anova", size = 6) +
  geom_text(data = letters.df.s, 
            aes(x = Stream, y = Placement.Value, label = Letter), 
            color = "black", hjust = -0.5,
            vjust = -0.8, fontface = "bold", size = 6) +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))
plot

png("box.png", width = 650, height = 550)
plot(plot)
dev.off()

###
##### Spring Kruskall
### pE, rich.E.less.B, rich.EPT, VASCI

Result = dunnTest(VASCI~Stream, data=metrics.s,
                  method="bonferroni")$res

X = Result$P.adj <= 0.05
names(X) = gsub(" ",  "", Result$Comparison)
X <- multcompLetters(X)
X <- as.data.frame(X$Letters) 
colnames(X) <- c("Letters")

placement <- metrics.s %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(VASCI)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.s <- cbind(X, placement)
unique(letters.df.s)

plot <- metrics.s %>% 
  ggplot(aes(x = group, y = VASCI, color = Impact)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) +
  geom_point(size = 4) +
  geom_point(data = longterm.s, aes(Stream, VASCI), color = "black") + 
  theme_classic() + #+ #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Spring)") +
  ylab("VASCI") +
  stat_compare_means(method = "kruskal", size = 6) +
  geom_text(data = letters.df.s, 
            aes(x = Stream, y = Placement.Value, label = Letters), 
            color = "black", hjust = -0.5,
            vjust = -0.8, fontface = "bold", size = 6) +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))
plot

png("box.png", width = 650, height = 550)
plot(plot)
dev.off()

###
##### Fall Kruskall
### pE, rich.E.less.B, rich.EPT, rich.Cling, VASCI

Result = dunnTest(VASCI~Stream, data=metrics.f,
                  method="bonferroni")$res

X = Result$P.adj <= 0.05
names(X) = gsub(" ",  "", Result$Comparison)
X <- multcompLetters(X)
X <- as.data.frame(X$Letters) 
colnames(X) <- c("Letters")

placement <- metrics.f %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(VASCI)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.f <- cbind(X, placement)
unique(letters.df.f)

plot <- metrics.f %>% 
  ggplot(aes(x = group, y = VASCI, color = Impact)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) +
  geom_point(size = 4) +
  geom_point(data = longterm.f, aes(Stream, VASCI), color = "black") + 
  theme_classic() + #+ #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Fall)") +
  ylab("VASCI") +
  stat_compare_means(method = "kruskal", size = 6) +
  geom_text(data = letters.df.f, 
            aes(x = Stream, y = Placement.Value, label = Letters), 
            color = "black", hjust = -0.5,
            vjust = -0.8, fontface = "bold", size = 6) +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))
plot

png("box.png", width = 650, height = 550)
plot(plot)
dev.off()

###
##### SC Kruskall
### sc.uScm

chem <- read.csv("chem.f21-s22.notrib.reduced.csv") %>%
  filter(season == "Fall")
chem$group <- factor(chem$Stream, levels = c("EAS (R)", "CRO (R)", 
                                                       "FRY (L)", "SPC (L)", 
                                                       "ROL (H)", "LLW (H)"))
longterm <- chem %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))

Result = dunnTest(sc.uScm~Stream, data=chem,
                  method="bonferroni")$res

X = Result$P.adj <= 0.05
names(X) = gsub(" ",  "", Result$Comparison)
X <- multcompLetters(X)
X <- as.data.frame(X$Letters) 
colnames(X) <- c("Letters")

placement <- chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(sc.uScm)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df <- cbind(X, placement)
unique(letters.df)

plot <- chem %>% 
  ggplot(aes(x = group, y = sc.uScm, color = Impact)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) +
  geom_point(size = 4) +
  geom_point(data = longterm, aes(Stream, sc.uScm), color = "black") + 
  theme_classic() + #+ #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Fall)") +
  ylab("SC (uS/cm)") +
  stat_compare_means(method = "kruskal", size = 6) +
  geom_text(data = letters.df, 
            aes(x = Stream, y = Placement.Value, label = Letters), 
            color = "black", hjust = -0.5,
            vjust = -0.8, fontface = "bold", size = 6) +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))
plot

png("box.png", width = 650, height = 550)
plot(plot)
dev.off()

#################################
######   Correlations   #########
#################################
library(scales) # for percent

Similar <- read_csv("Similarity.csv") 
Similar$group <- factor(Similar$Stream, levels = c("EAS (R)", "CRO (R)", 
                                                       "FRY (L)", "SPC (L)", 
                                                       "ROL (H)", "LLW (H)"))

plot <- Similar %>% ggplot(aes(x=sc.delt , y=SIM.b, group = Season, color = Season))+
  geom_point()  +
  facet_wrap("group") +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman", label.y.npc = 0.35, size = 5) +
  stat_regline_equation( label.y.npc = 0.15, size =5 ) +
  theme_classic()+
  xlab("Change in SC (uS/cm) Between Samples") + # (uS/cm)
  ylab("Bray-Curtis Similarity") +
  #labs(title = "Bray-Curtis Community Similarity vs \n Change in Specific Conductance") +
  theme( axis.title.x =element_text(size = 20), axis.text.x = element_text(size = 18),
        axis.title.y =element_text(size = 20), axis.text.y  = element_text(size = 18),
        strip.text.x = element_text(size = 20), 
        #strip.background =element_rect(fill=c("blue","blue","chartreuse3",
                                              #"chartreuse3", "chocolate2", 
                                              #"chocolate2")),
        legend.text=element_text(size=20), legend.title=element_text(size=25)) 
 
plot

#set size and save plot as png
png("ex.png", width = 750, height = 650)
plot(plot)
dev.off()

###
### Taxa Metrics
#### pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom

metrics <- read.csv("metrics.f21-s22.csv")
metrics$group <- factor(metrics$Stream, levels = c("EAS (R)", "CRO (R)", 
                                                   "FRY (L)", "SPC (L)", 
                                                   "ROL (H)", "LLW (H)"))



plot <- metrics %>% ggplot(aes(x=sc.uScm , y=p5dom, group = Season, color = Season))+
  geom_point(size = 4)  +
  facet_wrap("group") +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman", label.y.npc = 1, label.x.npc = 0.06, size = 5) +
  stat_regline_equation( label.y.npc = 0.1, label.x.npc = 0, size =5 ) +
  theme_classic()+
  xlab("Specific Conductance (uS/cm)") + # (uS/cm)
  ylab("% 5 Dominant Taxa") +
  #ylim(NA, 105) +
  #labs(title = "Bray-Curtis Community Similarity vs \n Change in Specific Conductance") +
  theme( axis.title.x =element_text(size = 20), axis.text.x = element_text(size = 18),
         axis.title.y =element_text(size = 20), axis.text.y  = element_text(size = 18),
         strip.text.x = element_text(size = 20), 
         #strip.background =element_rect(fill=c("blue","blue","chartreuse3",
         #"chartreuse3", "chocolate2", 
         #"chocolate2")),
         legend.text=element_text(size=20), legend.title=element_text(size=25)) 

plot

#set size and save plot as png
png("ex.png", width = 750, height = 650) #wider?
plot(plot)
dev.off()

#Global
# pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom

plot <- metrics %>% ggplot(aes(x=sc.uScm , y=pE))+
  geom_point(size = 4, aes(shape = Season))  +
  facet_wrap("Season") +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman", label.y.npc = 1, label.x.npc = 0, size = 7.5) +
  stat_regline_equation( label.y.npc = 1, label.x.npc = 0.5, size =7.5) +
  theme_classic()+
  xlab("Specific Conductance (uS/cm)") + # (uS/cm)
  ylab("% E") +
  #ylim(NA, 105) +
  #labs(title = "Bray-Curtis Community Similarity vs \n Change in Specific Conductance") +
  theme( axis.title.x =element_text(size = 20), axis.text.x = element_text(size = 18),
         axis.title.y =element_text(size = 20), axis.text.y  = element_text(size = 18),
         strip.text.x = element_text(size = 20), 
         #strip.background =element_rect(fill=c("blue","blue","chartreuse3",
         #"chartreuse3", "chocolate2", 
         #"chocolate2")),
         legend.text=element_text(size=20), legend.title=element_text(size=25)) 

plot

#set size and save plot as png
png("ex.png", width = 690, height = 520) #750 and 650 for taxa
plot(plot)
dev.off()




