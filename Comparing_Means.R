# Creating boxplots
# Perform One-way Anova and find pairwise values
library(tidyverse)
library(vegan) #decostand(), 
library(ggpubr) #stat_compare_means
library("multcompView")
library(car) #leveneTest

metrics <- read.csv("metrics.f21-s22.csv") 
metrics <- metrics %>%
  mutate(n.pSC = sqrt(pSC)) %>%
  mutate(n.rich.PR = log(rich.PR))
metrics.s <- filter(metrics, Season == "Spring")
metrics.s$group <- factor(metrics.s$Stream, levels = c("EAS (R)", "CRO (R)", 
                                                       "FRY (MI)", "SPC (MI)", 
                                                       "ROL (MI)", "LLW (MI)"))
metrics.f <- filter(metrics, Season == "Fall")
metrics.f$group <- factor(metrics.f$Stream, levels = c("EAS (R)", "CRO (R)", 
                                                       "FRY (MI)", "SPC (MI)", 
                                                       "ROL (MI)", "LLW (MI)"))

longterm <- metrics.f %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))

###### Transform and Test normality (normality =shapiro and variance with levene) #####
library(MASS)

trans <- apply(all.s, 2, shapiro.test)
df<- data.frame(t(sapply(trans,c)))

master.s <- master.s %>%mutate(across(c("Hshannon.mac" , "rich"     ,     "rich.E"   ,     "rich.EPT" ,    
                    "pEPT.less.H" ,  "pEPT"  ,        "pE"   ,         "rich.E.less.B",
                     "pE.less.B"  ,   "pP"      ,      "pP.less.Amph" , "pT"     ,      
                     "pPT.H"    ,     "pChi"    ,      "pD"   ,         "pChiO"  ,      
                   "p1dom"    ,     "p5dom"   ,      "rich.SC"   ,    "pPR"     ,     
                     "rich.Cling" ,   "pCling"   ,     "pSprawl"   ,    "pSwimm"  ,     
                     "ph"        ,    "temp.C"   ,     "sc.uScm"   ,    "no2no3.n.mgl" ,
                     "npoc.mgl"   ,   "na.mgl"   ,     "ca.mgl"  ,      "fe.ugl"      , 
                    "ni.ugl"     ,   "se.ugl"    ,    "ba.ugl"   ,     "u.ugl"  ,      
                     "ca.mg"     ,    "pfines"     ,   "D16"       ,    "D90"   ,       
                    "Hsimpson.peb" , "avgwetwidth"  , "avgembedd"   ,  "avgvegprotecL",
                    "avgvegprotecR" ,"avg.slope"), function(x)  boxcox(lm(x ~ 1) )))
# sqrt transformation
master.s <- cbind(master.s[1], master.s[,1:46]^(1/2)) %>%
  select(c(2:47))
master.s$Stream <- springrow$Stream

master.f <- cbind(master.f[1], master.f[,1:46]^(1/2)) %>%
  select(c(2:47))
master.f$Stream <- fallrow$Stream

# log transformation
master.s <- master %>%
  filter(Season == "Spring") %>%
  select(-c(Stream:Season, as.ugl)) %>%
  decostand(method = "log",  na.rm =TRUE) %>%
  select_if(~ ! any(is.na(.)))

master.f <- master %>%
  filter(Season == "Fall") %>%
  select(-c(Stream:Season)) %>%
  decostand(method = "log",  na.rm =TRUE) %>%
  select_if(~ ! any(is.na(.)))
master.f$Stream <- fallrow$Stream

shapiro.test(metrics.f$n.pSC)
# greater than 0.05 = normal
#library(car)
leveneTest(n.pSC~Stream, data=metrics.f)
#less than 0.05 means non-parametric (heterogeneity of variance)

##########    ANOVA    ##########
# https://www.r-bloggers.com/2021/08/how-to-perform-tukey-hsd-test-in-r/

# pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom
# Spring Anova: rich.SC, p5dom, rich.Cling; no transformations
# Fall Anova: rich.SC, p5dom; no transformations
leveneTest(n.rich.PR~Stream, data=metrics.s)
#less than 0.05 means non-parametric (heterogeneity of variance)


longterm <- metrics.s %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))

# Creating a letters dataframe for ANOVA
#https://www.mathiasecology.com/code/add-tukeys-significant-letters-to-ggplots
#library("multcompView")

###
##### SPRING Anova
###
set.seed(1045)
model <- aov(p5dom~Stream, data=metrics.s)
summary(model)
#less than 0.05 = there are differences
TukeyHSD(model, conf.level=.95)
plot(TukeyHSD(model, conf.level=.95), las = 2)

letters.df <- data.frame(multcompLetters(TukeyHSD(
  aov(p5dom~Stream, data=metrics.s))$Stream[,4])$Letters)

colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Stream <- rownames(letters.df) #Create column based on rownames

placement <- metrics.s %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(p5dom)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.s <- left_join(letters.df, placement, by = "Stream") #Merge dataframes
# This is exactly what I need for the table except instead of leaders it would have 
# sc.uScm and columns for all other metrics as well (also using kruskall values as Anova is not valid test for this data)
#library(ggpubr)

p5dom <- metrics.s %>% 
  ggplot(aes(x = group, y = p5dom, color = Impact)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) +
  geom_point(size = 3) +
  geom_point(data = longterm, aes(Stream, p5dom), color = "black") + 
  theme_classic() + #+ #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Spring)") +
  ylab("% 5 Dominant Taxa") +
  stat_compare_means(method = "anova", size = 6) +
  geom_text(data = letters.df.s, 
            aes(x = Stream, y = Placement.Value, label = Letter), 
            color = "black", hjust = -1.25,
            vjust = -0.8, fontface = "bold") +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))
rich.SC
rich.Cling
p5dom

png("VASCI.png")
plot(box.s)
dev.off()

###
##### Fall Anova
###
set.seed(1045)
model <- aov(rich.SC~Stream, data=metrics.f)
summary(model)
#less than 0.05 = there are differences
TukeyHSD(model, conf.level=.95)
plot(TukeyHSD(model, conf.level=.95), las = 2)

longterm <- metrics.f %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))

letters.df <- data.frame(multcompLetters(TukeyHSD(
  aov(rich.SC~Stream, data=metrics.f))$Stream[,4])$Letters)

colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Stream <- rownames(letters.df) #Create column based on rownames

placement <- metrics.f %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(rich.SC)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.s <- left_join(letters.df, placement, by = "Stream") #Merge dataframes
# This is exactly what I need for the table except instead of leaders it would have 
# sc.uScm and columns for all other metrics as well (also using kruskall values as Anova is not valid test for this data)
#library(ggpubr)

rich.SC <- metrics.f %>% 
  ggplot(aes(x = group, y = rich.SC, color = Impact)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) +
  geom_point(size = 3) +
  geom_point(data = longterm, aes(Stream, rich.SC), color = "black") + 
  theme_classic() + #+ #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Fall)") +
  ylab("Scraper Richness") +
  stat_compare_means(method = "anova", size = 6) +
  geom_text(data = letters.df.s, 
            aes(x = Stream, y = Placement.Value, label = Letter), 
            color = "black", hjust = -1.25,
            vjust = -0.8, fontface = "bold") +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))

rich.SC
p5dom

library(ggpubr)
box <- ggarrange(box.f, box.s, ncol = 2, common.legend = TRUE, legend = "bottom" )
box
# This graph is what I need but for the Kruskall-Wallis since Anova is not valid

png("box.aov.fall.ppr.png")
plot(box.f)
dev.off()

################################
#### Kruskall-Wallis non-parametric statistically valid option
####### see VASCI.R
library(tidyr)
library(FSA) #dunnTest

###
##### Spring Kruskall
###

longterm <- metrics.s %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))

kruskal.test(rich.EPT~Stream, data=metrics.s)
# less than 0.05 means there are differences among groups

Result = dunnTest(rich.EPT~Stream, data=metrics.s,
                  method="bonferroni")$res

### Use multcompView
#library(multcompView)

X = Result$P.adj <= 0.05
names(X) = gsub(" ",  "", Result$Comparison)
X <- multcompLetters(X)
X <- as.data.frame(X$Letters) 
colnames(X) <- c("Letters")

placement <- metrics.s %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(rich.EPT)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.s <- cbind(X, placement)
unique(letters.df.s)

rich.EPT <- metrics.s %>% 
  ggplot(aes(x = group, y = rich.EPT, color = Impact)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) +
  geom_point(size = 3) +
  geom_point(data = longterm, aes(Stream, rich.EPT), color = "black") + 
  theme_classic() + #+ #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Spring)") +
  ylab("Richness EPT") +
  stat_compare_means(method = "kruskal", size = 6) +
  geom_text(data = letters.df.s, 
            aes(x = Stream, y = Placement.Value, label = Letters), 
            color = "black", hjust = -1.25,
            vjust = -0.8, fontface = "bold") +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))
pE
rich.E.less.B
rich.EPT

plot <- ggarrange(pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom,  
                  common.legend = TRUE, legend = "right")
plot

png("box.metric.s.png", width = 1200, height = 750)
plot(plot)
dev.off()
# less than 0.05 means they are different

###
#### Fall Kruskall
### 

longterm <- metrics.f %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))


kruskal.test(rich.Cling~Stream, data=metrics.f)
# less than 0.05 means there are differences among groups

Result = dunnTest(rich.Cling~Stream, data=metrics.f,
                  method="bonferroni")$res

### Use multcompView
#library(multcompView)

X = Result$P.adj <= 0.05
names(X) = gsub(" ",  "", Result$Comparison)
X <- multcompLetters(X)
X <- as.data.frame(X$Letters) 
colnames(X) <- c("Letters")

placement <- metrics.f %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(rich.Cling)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.s <- cbind(X, placement)
unique(letters.df.s)

rich.Cling <- metrics.f %>% 
  ggplot(aes(x = group, y = rich.Cling, color = Impact)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) +
  geom_point(size = 3) +
  geom_point(data = longterm, aes(Stream, rich.Cling), color = "black") + 
  theme_classic() + #+ #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Fall)") +
  ylab("RClinger Richness") +
  stat_compare_means(method = "kruskal", size = 6) +
  geom_text(data = letters.df.s, 
            aes(x = Stream, y = Placement.Value, label = Letters), 
            color = "black", hjust = -1.25,
            vjust = -0.8, fontface = "bold") +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))
pE
rich.E.less.B
rich.EPT
rich.Cling

plot <- ggarrange(pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom,  
                  common.legend = TRUE, legend = "right")
plot

png("box.metric.f.png", width = 1800, height = 1000)
plot(plot)
dev.off()

# Other
# With 2 way Anova I was able to to analyze seasons separately but could not 
# create the letters dataframe or find the non-parametric equivalent
library("ggpubr")
res.aov <- aov(sc.uScm ~ season * Stream, data = chem)

box <- ggboxplot(chem, x = "Stream", y = "sc.uScm", color = "season",
                 palette = c("#00AFBB", "#E7B800"))
box

require("dplyr")
group_by(master, Season, Stream) %>%
  summarise(
    count = n(),
    mean = mean(sc.uScm, na.rm = TRUE),
    sd = sd(sc.uScm, na.rm = TRUE)
  )

TukeyHSD(res.aov, which = "Stream")

#Check test validity
library(car)
leveneTest(sc.uScm ~ season*Stream, data = chem)
#p must be more than 0.05 to be valid


# boxplots of SC
library(tidyverse)
chem <- read.csv("chem.f21-S22.notrib.reduced.csv")
longterm <- master.s %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))


box <- master.s %>% ggplot(aes(Stream, sc.uScm )) +
  geom_boxplot(aes(colour = Stream), outlier.colour = "black") +
  geom_point(data = longterm, aes(Stream, sc.uScm, color = Stream)) +
  theme_classic() + #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Spring)") +
  ylab("SC (uS/cm)")

box

png("box.sc.png", width = 700, height = 450)
plot(box)
dev.off()
