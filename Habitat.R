###Determining pebble count quantiles

#Background (Based on Damion Drover Thesis)
# 1) Variables chosen to represent sediment characteristics were quantiles of the sediment
#distributions (10, 16, 25, 50, 75, 84, and 90; expressed as D10, D16, etc.), ratio of D84 to D50, and interquartile range of sediment distributions
# 2) proportion of sediment that occur in each size class (e.g. % fines, % gravel, % cobble) based on the
# Udden-Wentworth grain size classification system (Wentworth 1922, Blair and McPherson 1999), 
#ratio of % small cobble (64-128 mm) 
#to % fines (<2 mm), 
#ratio of %large cobble (128-255 mm) to
#% fines (<2 mm) (LCF), 
# 3) Simpsons diversity of size class distribution

# Why did Damion choose his?

#### 1) intial data exploration and determining quantiles
# load packages and read in pebble count tally data
library(tidyverse)
library(matrixStats)
library(tidyr)
library("vegan")
library("ggpubr")

pebraw <- read_csv("pebblecounts.csv") 
pebraw <- select(pebraw, c(1:16))
peb <- pivot_longer(pebraw, c(2:16)) 
num <- as.numeric(peb$name)
peb$pebble_size <- num

# now you can look at the raw data
ggplot(peb, aes(x = name, y = percent, fill = Site, position = "dodge")) +
  geom_boxplot(stat = "identity") +
  facet_wrap("Stream") +
  ylab("Percent") +
  xlab("Pebble size")+        
  scale_x_discrete(limits = c("2", "2.8", "4", "5.6", "8", "11", 
                              "16", "22.6", "32", "45"
                              ,"64", "90", "128", "180"))

# creating quantiles
pebquartiles <- peb %>%
  drop_na() %>%
  uncount(value) %>% #function has each size pebble show up as many times as it occurs
  group_by(Site) %>%
  summarise(D10 = quantile( pebble_size, .10, na.rm = TRUE ),
            D16 = quantile( pebble_size, .16, na.rm = TRUE ),
            D25 = quantile( pebble_size, .25, na.rm = TRUE ),
            D50 = median( pebble_size, na.rm = TRUE ),
            D75 = quantile( pebble_size, .75, na.rm = TRUE ),
            D84 = quantile( pebble_size, .845, na.rm = TRUE ),
            D90 = quantile( pebble_size, .90, na.rm = TRUE ))

pebquartiles <- pebquartiles %>%
  mutate(D84to50 = D84/D50) %>%
  mutate(IQR = D75-D25)

####2) Determining what percent occur in each size class
# Start by dividing the data into categories

pebcat<- peb %>%
  mutate(category = case_when(pebble_size <= 4 ~ "fines",
                                pebble_size > 4 &
                              pebble_size < 64 ~ "pebbles", 
                                pebble_size >= 64 &
                              pebble_size < 128  ~ "smallcobble", 
                                pebble_size >= 128 &
                              pebble_size <= 4000 ~ "largecobble"))
         
# determine percents
peb2adj <- aggregate(value ~ category + Site, data = pebcat, FUN = sum, na.rm = TRUE) %>%
  pivot_wider(names_from = category, values_from = value) %>%
  mutate(pfines = (fines/(largecobble + pebbles + smallcobble +fines)) * 100,
        ppebbles = (pebbles/(largecobble + fines + smallcobble + pebbles )) * 100,
        psmallcobble = (smallcobble/(largecobble + pebbles + fines + smallcobble )) * 100,
        plargecobble = (largecobble/(fines + pebbles + smallcobble + largecobble)) * 100)
# next determine smallcobble:fines and largecobble:fines
  peb2adj <- peb2adj %>% 
    mutate(LCF = plargecobble/pfines,
       SCF = psmallcobble/pfines)
  
#####3) Find simpson diversity of size class distribution
  pebraw2 <- pebraw[,-1]
  Hsimpson<- diversity(pebraw2, index="simpson")

#join peb2, pebquartiles, and Hsimpson
  peb2adj$Hsimpson<- Hsimpson
  habitatmaster <- left_join(peb2adj, pebquartiles, by = "Site")


# transect <- read_csv("habtransect.csv")
habitatmaster <- read_csv("habitatmaster.csv")
distance <- read_csv("distance.csv")
habitatmaster <- left_join( habitatmaster, distance , by = c("Site" = "site"))
    
write.csv(habitatmaster, file="habitatmasteradj.csv", sep = ",")
 
ggplot(habitatmaster, aes(x = dist.d, y = D25, colour = Stream )) +
  geom_line() +
  facet_wrap("Stream") +
  stat_cor(method = "spearman") +
  theme_classic()

## using transect avarages to add to habitatmaster

library(dplyr)

transect <- read_csv("habtransect.csv") 
  #filter(transect != "EAS9B") 
  #select(-c(2, 10)) #- site, transect, and channel code; non-numerics

#Option to assign groups
#mydf <- data.frame("transect"=seq(1:138))
#GroupLabels <- 0:(nrow(mydf) - 1) %/%  3
#transect$Group <- GroupLabels

unique(transect)

Avgs <- transect %>% group_by(site) %>% summarize(avgwetwidth = mean(wetwid), 
  avgdepth = mean((lft+  lctr+ ctr+ rctr+ rht)/5), 
  avgembedd  = mean(embedd), avgbankstabL = mean(bankstabL), 
  avgbankstaR = mean(bankstaR), avgriparianwidL = mean(riparianwidL), 
  avgriparianwidR = mean(riparianwidR), avgvegprotecL = mean(vegprotecL), 
  avgvegprotecR = mean(vegprotecR), avgcancov = mean(cancov))
#`summarise()` ungrouping output (override with `.groups` argument)
Avgs 

#habitatmaster <- read_csv("habitatmaster.csv")
  
habitatmaster <- left_join(habitatmaster, Avgs, by = c("Site" = "site"))
 
habitatmaster  

#write.csv(habitatmaster, file="habitatmaster.csv", sep = ",")
