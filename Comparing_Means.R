# Creating boxplots
# Perform One-way Anova and find pairwise values
library(tidyverse)
library(vegan) #decostand(), 
library(ggpubr) #stat_compare_means

master <- read.csv("reduced.all.data.csv")

chem <- read.csv("chem.f21-S22.notrib.reduced.csv") #%>%
  filter(season == "Spring") %>%
  select(-c( 35))
habitat <- read.csv("habitatmaster.csv")
metrics <- read.csv("metrics.f21-s22.csv")%>%
  filter(Season == "Spring")

all.s <- left_join(metrics, chem,  by = "Site") %>%
  left_join(habitat,  by = "Site") %>%
  select(-c(1:6, 66:69, 105)) #%>%
  #decostand(method = "log",  na.rm =TRUE) 
all.s <- cbind(all.s[1], all.s[,1:125]^(1/2))# %>%
  select(c(2:126))


master.s <- master %>%
  filter(Season == "Spring") #%>%
  #filter(Stream == "LLW (MI)") %>%
  select(-c(Stream:Season, as.ugl)) %>%
  select_if(~ ! any(is.na(.)))
springrow <-  read.csv("reduced.all.data.csv") %>%
  filter(Season == "Spring") %>%
  select(Stream) 
master.s$Stream <- springrow$Stream

master.f <- master %>%
  filter(Season == "Fall") %>%
  select(-c(Stream:Season)) %>%
  select_if(~ ! any(is.na(.)))
fallrow <-  read.csv("reduced.all.data.csv") %>%
  filter(Season == "Fall") %>%
  select(Stream) 
master.f$Stream <- fallrow$Stream

#test normaility with shapiro and variance with levene
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

shapiro.test(metrics$VASCI)
# greater than 0.05 = normal
library(car)
leveneTest(sc.uScm~Stream, data=master.s)
#less than 0.05 means non-parametric (heterogeneity of variance)

p<-ggplot(master.s, aes(x=pPR)) + 
  geom_histogram(bins = 15)
p

#ANOVA
# https://www.r-bloggers.com/2021/08/how-to-perform-tukey-hsd-test-in-r/
library(car)
leveneTest(pPR~Stream, data=master.s)
#less than 0.05 means non-parametric (heterogeneity of variance)
# Anova Test for SC not valid because data is not normal 

set.seed(1045)
model <- aov(VASCI~Stream, data=metrics)
summary(model)
#less than 0.05 = there are differences
TukeyHSD(model, conf.level=.95)
plot(TukeyHSD(model, conf.level=.95), las = 2)

# Creating a letters dataframe for ANOVA
#https://www.mathiasecology.com/code/add-tukeys-significant-letters-to-ggplots
library("multcompView")

# Spring Anova
letters.df <- data.frame(multcompLetters(TukeyHSD(
  aov(VASCI~Stream, data=metrics))$Stream[,4])$Letters)

colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Stream <- rownames(letters.df) #Create column based on rownames

placement <- metrics %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(VASCI)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.s <- left_join(letters.df, placement) #Merge dataframes
# This is exactly what I need for the table except instead of leaders it would have 
# sc.uScm and columns for all other metrics as well (also using kruskall values as Anova is not valid test for this data)
library(ggpubr)

box.s <- metrics %>% #Dataframe from which data will be drawn
  ggplot(aes(x = reorder(Stream, VASCI, median), y = VASCI, color = Stream)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) + #I like to set the color of boxplots to black with the alpha at 0 (fully transparent). I also like geom_jitter() but do not use it here for simplicity.
  theme_classic() + #Clean, minimal theme courtesy of the "egg" package
  geom_point(data = longterm, aes(Stream, VASCI, color = Stream)) +
  xlab("Stream (Spring)") +
  ylab("VASCI") +
  stat_compare_means(method = "anova")+
  geom_text(data = letters.df.s, aes(x = Stream, y = Placement.Value, label = 
                                       Letter), color = "black", hjust = -1.25,
            vjust = -0.8, fontface = "bold")
box.s


png("VASCI.png")
plot(box.s)
dev.off()

# Fall Anova
letters.df <- data.frame(multcompLetters(TukeyHSD(
  aov(pChiO~Stream, data=master.f))$Stream[,4])$Letters)

colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Stream <- rownames(letters.df) #Create column based on rownames

placement <- master.f %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
summarise(quantile(pChiO)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.f <- left_join(letters.df, placement) #Merge dataframes

box.f <- master.f %>% #Dataframe from which data will be drawn
  ggplot(aes( x = reorder(Stream, pChiO, median), y = pChiO, color = Stream)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) + #I like to set the color of boxplots to black with the alpha at 0 (fully transparent). I also like geom_jitter() but do not use it here for simplicity.
  theme_classic() + #Clean, minimal theme courtesy of the "egg" package
  stat_summary(fun.y="mean")+ 
  xlab("Stream (Fall)") +
  #ylab("Percent Composed of 5 Dominant Taxa") +
  stat_compare_means(method = "anova")+
  geom_text(data = letters.df.f, aes(x = Stream, y = Placement.Value, label = Letter), size = 4, color = "black", hjust = -1.25, vjust = -0.8, fontface = "bold")
box.f

library(ggpubr)
box <- ggarrange(box.f, box.s, ncol = 2, common.legend = TRUE, legend = "bottom" )
box
# This graph is what I need but for the Kruskall-Wallis since Anova is not valid

png("box.aov.fall.ppr.png")
plot(box.f)
dev.off()

# Kruskall-Wallis non-parametric statistically valid option

# Spring Kruskall
kruskal.test(sc.uScm ~ Stream, data = master)
# less than 0.05 means there are differences among groups
library(tidyr)
library(FSA)

Result = dunnTest(sc.uScm ~ Stream, data = master,
                  method="bonferroni")$res


### Use cldList()

library(rcompanion)

X <- cldList(sc.uScm ~ Stream, data = master.s)

### Use multcompView
#library(multcompView)
#X = Result$P.adj <= 0.05
#names(X) = gsub(" ",  "",  Result$Comparison)
#multcompLetters(X)

placement <- master.s %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(sc.uScm)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.s <- cbind(X, placement)

#p.s <- data.frame(pairwise.wilcox.test(master.s$sc.uScm, master.s$Stream,
                                    # p.adjust.method = "BH")$p.value)


#compare_means(sc.uScm ~ Stream, data = master.s)
# only LLW-ROL was ns so they are the only value I graphed although ideally I need a way to make all others a and these 2 b like with the anova

#https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
#my_comparisons <- list( c("EAS (R)", "CRO (R)"), c("EAS (R)", "FRY (MI)"), c("EAS (R)", "ROL (MI)"),
                        #c("EAS (R)", "SPC (MI)"), c("EAS (R)", "LLW (MI)"), c("CRO (R)", "FRY (MI)"),
                        #c("CRO (R)", "ROL (MI)"), c("CRO (R)", "SPC (MI)"), c("CRO (R)", "LLW (MI)"),
                        #c("FRY (MI)", "ROL (MI)"), c("FRY (MI)", "SPC (MI)"), c("FRY (MI)", "LLW (MI)"), 
                        #c("ROL (MI)", "SPC (MI)"), c("ROL (MI)", "LLW (MI)"), c("SPC (MI)", "LLW (MI)"))

box.s <- master.s %>% #Dataframe from which data will be drawn
  ggplot(aes(x = reorder(Stream, sc.uScm, median), y = sc.uScm, color = Stream)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) + #I like to set the color of boxplots to black with the alpha at 0 (fully transparent). I also like geom_jitter() but do not use it here for simplicity.
  theme_classic() #+ #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Spring)") +
  #ylab("Clinger Richness") +
  stat_compare_means(method = "kruskal")+
  geom_text(data = letters.df.s, aes(x = Stream, y = Placement.Value, label = 
                                       Letter), color = "black", hjust = -1.25,
            vjust = -0.8, fontface = "bold")
box.s

png("box.krus.peptlessh.xstreamxspring.png")
plot(box.s)
dev.off()
# less than 0.05 means they are different

#Fall Kruskall
# Spring
kruskal.test(rich.EPT ~ Stream, data = master.f)
# less than 0.05 means there are differences among groups
library(tidyr)


box.f <- master.f %>% #Dataframe from which data will be drawn
  ggplot(aes(x = Stream, y = sc.uScm, color = Stream)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) + #I like to set the color of boxplots to black with the alpha at 0 (fully transparent). I also like geom_jitter() but do not use it here for simplicity.
  theme_classic() + #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Fall)") +
  ylab("SC (uS/cm)") +
  stat_compare_means() +
  #stat_compare_means(comparisons = my_comparisons)
  stat_compare_means(comparisons = list
                     (c("ROL (MI)", "SPC (MI)")), label.y = 1500)
box.f

box <- ggarrange(box.f, box.s, ncol = 2, common.legend = TRUE, legend = "bottom" )

png("box.krus+wilcox.scxstreamxseason.png", width = 900, height = 450)
plot(box)
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
