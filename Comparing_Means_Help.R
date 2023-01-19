# Creating boxplots
# Perform One-way Anova and find pairwise values
library(tidyverse)
library(vegan) #decostand(), 

master <- read.csv("reduced.all.data.csv")
master.s <- master %>%
  filter(Season == "Spring") %>%
  select(-c(Stream:Season, as.ugl)) %>%
  #decostand(method = "log",  na.rm =TRUE) %>%
  select_if(~ ! any(is.na(.)))
springrow <-  read.csv("reduced.all.data.csv") %>%
  filter(Season == "Spring") %>%
  select(Stream) 
master.s$Stream <- springrow$Stream

master.f <- master %>%
  filter(Season == "Fall")


#transforming data
library(MASS)

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
df_cube_root <- cbind(master.s[1], master.s[,1:46]^(1/5))
boxcox(lm(master.s$sc.uScm ~ 1))
apply(df_cube_root,2,shapiro.test)

p<-ggplot(master.s, aes(x=sc.uScm)) + 
  geom_histogram(color="black", fill="white")
p

#ANOVA
# https://www.r-bloggers.com/2021/08/how-to-perform-tukey-hsd-test-in-r/
set.seed(1045)
model <- aov(rich.Cling~Stream, data=master.s)
summary(model)
#less than 0.05 = there are differences
TukeyHSD(model, conf.level=.95)
plot(TukeyHSD(model, conf.level=.95), las = 2)

library(car)
leveneTest(pPR~Stream, data=master.s)
#less than 0.05 means non-parametric (heterogeneity of variance)
# Anova Test for SC not valid because data is not normal 

# Creating a letters dataframe for ANOVA
#https://www.mathiasecology.com/code/add-tukeys-significant-letters-to-ggplots
library("multcompView")

# Spring Anova
letters.df <- data.frame(multcompLetters(TukeyHSD(
  aov(rich.Cling~Stream, data=master.s))$Stream[,4])$Letters)

colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Stream <- rownames(letters.df) #Create column based on rownames

placement <- master.s %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(rich.Cling)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.s <- left_join(letters.df, placement) #Merge dataframes
# This is exactly what I need for the table except instead of leaders it would have 
# sc.uScm and columns for all other metrics as well (also using kruskall values as Anova is not valid test for this data)
library(ggpubr)

box.s <- master.s %>% #Dataframe from which data will be drawn
  ggplot(aes(x = Stream, y = rich.Cling, color = Stream)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) + #I like to set the color of boxplots to black with the alpha at 0 (fully transparent). I also like geom_jitter() but do not use it here for simplicity.
  theme_classic() + #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Spring)") +
  ylab("Clinger Richness") +
  stat_compare_means(method = "anova")+
  geom_text(data = letters.df.s, aes(x = Stream, y = Placement.Value, label = 
                                       Letter), color = "black", hjust = -1.25,
            vjust = -0.8, fontface = "bold")
box.s

png("box.tukeyhsd+aov.clingrichxstreamxspring.png")
plot(box.s)
dev.off()

# Fall Anova
letters.df <- data.frame(multcompLetters(TukeyHSD(
  aov(sc.uScm~Stream, data=master.f))$Stream[,4])$Letters)

colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Stream <- rownames(letters.df) #Create column based on rownames

placement <- master.f %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
summarise(quantile(sc.uScm)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.f <- left_join(letters.df, placement) #Merge dataframes

box.f <- master.f %>% #Dataframe from which data will be drawn
  ggplot(aes(x = Stream, y = sc.uScm, color = Stream)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) + #I like to set the color of boxplots to black with the alpha at 0 (fully transparent). I also like geom_jitter() but do not use it here for simplicity.
  theme_classic() + #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Fall)") +
  ylab("SC (uS/cm)") +
  stat_compare_means(method = "anova")+
  geom_text(data = letters.df.f, aes(x = Stream, y = Placement.Value, label = Letter), size = 4, color = "black", hjust = -1.25, vjust = -0.8, fontface = "bold")
box.f

library(ggpubr)
box <- ggarrange(box.f, box.s, ncol = 2, common.legend = TRUE, legend = "bottom" )
box
# This graph is what I need but for the Kruskall-Wallis since Anova is not valid

png("box.tukeyhsd+aov.scxstreamxseason.png", width = 900, height = 450)
plot(box)
dev.off()

# Kruskall-Wallis non-parametric statistically valid option

# Spring Kruskall
kruskal.test(rich.E ~ Stream, data = master.s)
# less than 0.05 means there are differences among groups
library(tidyr)
library(FSA)

Result = dunnTest(rich.E ~ Stream,
                  data=master.s,
                  method="bonferroni")$res


### Use cldList()

library(rcompanion)

X <- cldList(P.adj ~ Comparison, data=Result)

### Use multcompView
#library(multcompView)
#X = Result$P.adj <= 0.05
#names(X) = gsub(" ",  "",  Result$Comparison)
#multcompLetters(X)

placement <- master.s %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(rich.E)[4])

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
  ggplot(aes(x = Stream, y = rich.E, color = Stream)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) + #I like to set the color of boxplots to black with the alpha at 0 (fully transparent). I also like geom_jitter() but do not use it here for simplicity.
  theme_classic() + #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Spring)") +
  ylab("Ephemeroptera Richness") +
  stat_compare_means(method = "kruskal")+
  geom_text(data = letters.df.s, aes(x = Stream, y = Placement.Value, label = 
                                       Letter), color = "black", hjust = -1.25,
            vjust = -0.8, fontface = "bold")
box.s

png("box.krus.shan.xstreamxspring.png")
plot(box.s)
dev.off()
# less than 0.05 means they are different

#Fall Kruskall
# Spring
kruskal.test(sc.uScm ~ Stream, data = master.f)
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
res.aov <- aov(sc.uScm ~ Season * Stream, data = master)

box <- ggboxplot(master, x = "Stream", y = "sc.uScm", color = "Season",
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
leveneTest(sc.uScm ~ Season*Stream, data = master)
#p must be more than 0.05 to be valid
