# OBJECTIVE 1: Determine with-in stream variability of AMI communities
# PART 1: Is it changing? A) Simple regression for each stream of community similarity vs 
# distance downstream. 
# PART 2: What is changing? Use step-wise multiple regression to determine which bug 
# metrics and/or taxa are (or are not) shifting in each stream.

library(dplyr)
#create table of means, sd, pairwise for each streams insect metrics 
library(tidyr)

metrics <- read.csv("metrics.f21-s22.csv")
habitat_raw <- read.csv("habitatmaster.csv")
habitat <- habitat_raw %>%
  select(c( 1, 7:34)) 
chem <- read.csv("chem.f21-s22.notrib.csv") 
  

# create table of means and standard deviations
table <- chem %>% group_by(interaction(Stream, season)) %>% 
  summarize_all(list(mean = ~mean(.), sd = ~sd(.)))
table <- pivot_longer(table, cols = c(pfines_mean:avgcancov_sd))
table <- pivot_wider(table, names_from = Stream)

library(ggpubr)
library(tidyverse)
library(broom)

#one.way <- aov(rich.E ~ Stream, data = metrics)
#summary(one.way)

#seperate metrics by site and join all together at end?

metrics.f.s.results <- NULL;
for(i in 7:ncol(metrics))
{
  
  column <- names(metrics[i])
  #tidy will summarise and return neat format
  avz <- broom::tidy(kruskal.test(metrics[, i] ~ Stream, data = metrics))
  
  # Add this condition if you only want aov with P < 0.05 printed
  #if(avz$p.value[1] < 0.05) {
  
  
  metrics.f.s.results <- rbind(metrics.f.s.results, avz) # results <- print(column) #metrics <- cbind(metrics, column)
  # results <- print(avz) #metrics <- cbind(metrics, avz) 
}

#a.results <- filter(a.results, term != "Residuals") only needed for anova

metrics.f.s.results$s.statistic <- metrics.s.results$statistic
metrics.f.s.results$s.p.value <- metrics.s.results$p.value

a.metrics <- pivot_longer(metrics, cols = Hshannon:pTOL)
a.metrics <- a.metrics[c(1:59),]

metrics.f.s.results$Metric <- a.metrics$name
write.csv(metrics.f.s.results, file="stream.kruskal.metrics.csv", sep = ",")

a.results <- filter(a.results, p.value < 0.05) 

#will compile files in excel (saved as Obj1.stream.means&anova)
#write.csv(table, file="stream.means.csv", sep = ",")
#write.csv(a.results, file="stream.anova.csv", sep = ",")

#correlation matrix of distances to bug metrics
dis.simm <- read.csv("dis.simm.f21.csv") %>%
  filter(START %in% c("CRO-1", "EAS-1", "FRY-1", "LLW-1", "SPC-1", "ROL-1"))
metrics <- read.csv("metrics.f21.csv")
bug.dist.cor <- left_join(metrics, dis.simm, by = c("Site" = "End"))
#bug.dist.cor <- left_join(bug.dist.cor, chem, by = "Site" ) # for sc to distance graph

library(dplyr)

# prepare data to be used (sc and all metrics)
c.df <- read.csv("metrics.f21-s22.csv") %>%
  filter( Season == "Spring")
chem <- read.csv("chem.f21-S22.notrib.csv") %>%
  filter( season == "Spring")
c.df <- left_join(c.df, chem, by = c("Site","Season" = "season", "Stream", "dist.d"))
c.df <- select(c.df, c(7:65, 70)) # 3 is dist, 70 is sc

# Split by stream and determine if data is normal use rcorr instead of cor to get p and cor value values

library(tibble)
library(Hmisc) #rcorr()

#CRO
CRO <- filter(c.df, Stream == "CRO (R)") %>%
  select(c(2:61))
cro <- rcorr(as.matrix(CRO), type = "spearman")
cro.r = data.frame(cro$r) %>%
  select("dist.d")
cro.p = data.frame(cro$P) %>%
  select("dist.d")
cro <- NULL
cro$p.spring <- cro.p$dist.d
cro$metrics <- cro.p$rownames
cro$r.spring <- cro.r$dist.d
cro <- as.data.frame(cro)
rownames <- rownames(cro.p)
cro <- add_column(cro, rownames)
cro <- filter(cro, p.spring <= 0.05)

#set.seed(123)
#norm.cro <- apply(cro,2,shapiro.test)
#norm.cro <- do.call(rbind.data.frame, norm.cro)
#norm.cro <- norm.cro %>%
#mutate("test" = case_when(p.value <= 0.05 ~ "Pearson", 
# p.value >0.05 ~ "Spearman"))

#filter(norm.cro, test == "Pearson")

#EAS
EAS <- filter(c.df, Stream == "EAS (R)") %>%
  select(c(2:61))
eas <- rcorr(as.matrix(EAS), type = "spearman")
eas.r = data.frame(eas$r) %>%
  select("dist.d")
eas.p = data.frame(eas$P) %>%
  select("dist.d")
eas <- NULL
eas$p.spring <- eas.p$dist.d
eas$metrics <- eas.p$rownames
eas$r.spring <- eas.r$dist.d
eas <- as.data.frame(eas)
rownames <- rownames(eas.p)
eas <- add_column(eas, rownames)
eas <- filter(eas, p.spring <= 0.05)

#FRY
FRY <- filter(c.df, Stream == "FRY (MI)") %>%
  select(c(2:61))
fry <- rcorr(as.matrix(FRY), type = "spearman")
fry.r = data.frame(fry$r) %>%
  select("dist.d")
fry.p = data.frame(fry$P) %>%
  select("dist.d")
fry <- NULL
fry$p.spring <- fry.p$dist.d
fry$metrics <- fry.p$rownames
fry$r.spring <- fry.r$dist.d
fry <- as.data.frame(fry)
rownames <- rownames(fry.p)
fry <- add_column(fry, rownames)
fry <- filter(fry, p.spring <= 0.05)

#LLW
LLW <- filter(c.df, Stream == "LLW (MI)") %>%
  select(c(2:61))
llw <- rcorr(as.matrix(LLW), type = "spearman")
llw.r = data.frame(llw$r) %>%
  select("dist.d")
llw.p = data.frame(llw$P) %>%
  select("dist.d")
llw <- NULL
llw$p.spring <- llw.p$dist.d
llw$metrics <- llw.p$rownames
llw$r.spring <- llw.r$dist.d
llw <- as.data.frame(llw)
rownames <- rownames(llw.p)
llw <- add_column(llw, rownames)
llw <- filter(llw, p.spring <= 0.05)

#ROL
ROL<- filter(c.df, Stream == "ROL (MI)") %>%
  select(c(2:61))
rol <- rcorr(as.matrix(ROL), type = "spearman")
rol.r = data.frame(rol$r) %>%
  select("dist.d")
rol.p = data.frame(rol$P) %>%
  select("dist.d")
rol <- NULL
rol$p.spring <- rol.p$dist.d
rol$metrics <- rol.p$rownames
rol$r.spring <- rol.r$dist.d
rol <- as.data.frame(rol)
rownames <- rownames(rol.p)
rol <- add_column(rol, rownames)
rol <- filter(rol, p.spring <= 0.05)

#SPC
SPC<- filter(c.df, Stream == "SPC (MI)") %>%
  select(c(2:61))
spc <- rcorr(as.matrix(SPC), type = "spearman")
spc.r = data.frame(spc$r) %>%
  select("dist.d")
spc.p = data.frame(spc$P) %>%
  select("dist.d")
spc <- NULL
spc$p.spring <- spc.p$dist.d
spc$metrics <- spc.p$rownames
spc$r.spring <- spc.r$dist.d
spc <- as.data.frame(spc)
rownames <- rownames(spc.p)
spc <- add_column(spc, rownames)
spc <- filter(spc, p.spring <= 0.05)

#Global
spc <- rcorr(as.matrix(c.df), type = "spearman")
spc.r = data.frame(spc$r) %>%
  select("sc.uScm")
spc.p = data.frame(spc$P) %>%
  select("sc.uScm")
spc <- NULL
spc$p.spring <- spc.p$sc.uScm
spc$metrics <- spc.p$rownames
spc$r.spring <- spc.r$sc.uScm
spc <- as.data.frame(spc)
rownames <- rownames(spc.p)
spc <- add_column(spc, rownames)
spc <- filter(spc, p.spring <= 0.05)

# join correlations
cor.scxstream <- full_join(eas, cro, by = "rownames") %>%
  full_join(llw, by = "rownames") %>%
  full_join(spc, by = "rownames") %>%
  full_join(rol, by = "rownames") %>%
  full_join(fry, by = "rownames") 
write.csv(cor.scxstream, file="cor.dist.metric.stream.<0.05.s22.csv", sep = ",")






