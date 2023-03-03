#######################################################################################
# VASCI
#**************************************************************************************
library(vegan)
library(tidyverse)

# Read and prep bug count and taxa list ----------------------
bugs <- read.csv("bugid.csv") %>%
  filter(Season == "Fall")
bugs <- bugs[, colSums(bugs != 0, na.rm = TRUE) > 0]

genera <- colnames(bugs[4:97])
#genera <- genera[-c(3,7,65,71)] #Spring
taxa <- read.csv("taxa.list.csv")
taxa <- taxa[taxa$taxon %in% genera ,]
gen.func <- read.csv("taxa.traits.csv") %>%
  drop_na()
gen.func <- gen.func[gen.func$taxon %in% genera ,]
fam.func <- read.csv('Family FFG and Tolerance.csv')

func <- left_join(taxa, fam.func, by = c("family" ="taxon"))
#write.csv(func, file = "bugidfam.csv")
# created new bug count df (bugs.fam) with family as column instead of taxa
# manually summed columns of same family

# Read in VASCI BSVs --------------------------
bsv <- read.csv('VASCI BSV.csv')

# Make table of FAMILY-wise ABUNDANCE -----------------------------------------------------
bugs.fam <- read.csv("bugidfam.csv") %>% 
  filter(Season == "Fall")
fams <- taxa$family

# Determine which taxa in samples belong to each order by referencing Taxa List ------------
E <- fams[which(fams %in% taxa$family[taxa$order=="Ephemeroptera"])]
P <- fams[which(fams %in% taxa$family[taxa$order=="Plecoptera"])]
Tr <- fams[which(fams %in% taxa$family[taxa$order=="Trichoptera"])]
C <- fams[which(fams %in% taxa$family[taxa$order=="Coleoptera"])]
D <- fams[which(fams %in% taxa$family[taxa$order=="Diptera"])]
# If genus in sample, make it 1, sum by row for each ORDER
rich.E <- rowSums(ifelse(bugs.fam[E]>0,1,0))
rich.P <- rowSums(ifelse(bugs.fam[P]>0,1,0))
rich.T <- rowSums(ifelse(bugs.fam[Tr]>0,1,0))
rich.C <- rowSums(ifelse(bugs.fam[C]>0,1,0))
rich.D <- rowSums(ifelse(bugs.fam[D]>0,1,0))

# Calc raw metrics -------------------------------------------------------------------------
# Total No. Individuals per sample for proportion calculations ----------------------
totind <- rowSums(bugs.fam[4:55])

# 1) FAMILY RICHNESS
tot.rich <- rowSums(ifelse(bugs.fam[8:length(bugs.fam)]>0,1,0))

# 2) EPT RICHNESS
EPT.rich <- rich.E + rich.P + rich.T

# 3) Percent E
pE <- 100*rowSums(bugs.fam[E])/totind

# 4) Percent PT-Hydropsychidae
PT <- c(P,Tr)
PT.H <- PT[-which(PT=="Hydropsychidae")]
pPT.H <- 100*rowSums(bugs.fam[PT.H])/totind

# 5) Percent Scrapers
scrapers <- fams[which(fams %in% fam.func$taxon[fam.func$ffg =='SC'])]
pScrap <- 100*rowSums(bugs.fam[scrapers])/totind

# 6) Percent Chironomidae
pChiron <- 100*rowSums(bugs.fam['Chironomidae'])/totind

# 7) Percent 2 Dominant Taxa
dom2 <- NULL
for(i in 1:length(bugs.fam)){
  a <- sort(bugs.fam[i,8:length(bugs.fam)],decreasing=TRUE)
  dom2[i] <- sum(a[1:2])
}
dom2
p2dom <- 100*dom2/totind
p2dom

# 8) HBI - (taxon abund x taxon tol val) / tot abund ------------------
#bugs.fam.n <- read.csv("bugidfam.s.csv")
bugs.fam.sum <- bugs.fam %>%
  mutate(bugs.sum = rowSums(bugs.fam[4:55])) %>%
  select(Site, Stream, Season, bugs.sum)

bugs.fam.h <- pivot_longer(bugs.fam, cols = !c(Site, Stream, Season)) %>%
  left_join(fam.func, by = c("name" = "taxon")) %>%
  mutate(h = value*tv ) %>%
  drop_na() %>%
  select(Site, Stream, Season, name, h) %>%
  pivot_wider(names_from = name, values_from = h) %>%
  as.data.frame()
 
bugs.fam.h
bugs.fam.hbi <- bugs.fam.h %>%
  mutate(hsum = rowSums(bugs.fam.h[4:50])) %>%
  left_join(bugs.fam.sum, by = c("Site", "Season", "Stream")) %>%
  mutate(hbi = hsum/bugs.sum) 

hbi <- bugs.fam.hbi$hbi

#na_rows <-bugs.fam.long[!complete.cases(bugs.fam.long), ]
#na_rows

# Calc VASCI metric scores by comparing raw metric to BSV ------------------------------
# 1) FAMILY RICHNESS
tot.rich.score <- 100*tot.rich/bsv$BSV[bsv$metric=='tot.rich']       # full value
tot.rich.score100 <- ifelse(tot.rich.score > 100,100,tot.rich.score) # 100% ceiling

# 2) EPT RICHNESS
EPT.rich.score <- 100*EPT.rich/bsv$BSV[bsv$metric=='EPT.rich']       # full value
EPT.rich.score100 <- ifelse(EPT.rich.score > 100,100,EPT.rich.score) # 100% ceiling

# 3) Percent E
pE.score <- 100*pE/bsv$BSV[bsv$metric=='pE']       # full value
pE.score100 <- ifelse(pE.score > 100,100,pE.score) # 100% ceiling

# 4) Percent PT-Hydropsychidae
pPT.H.score <- 100*pPT.H/bsv$BSV[bsv$metric=='pPT.H']       # full value
pPT.H.score100 <- ifelse(pPT.H.score > 100,100,pPT.H.score) # 100% ceiling

# 5) Percent Scrapers
pScrap.score <- 100*pScrap/bsv$BSV[bsv$metric=='pScrap']       # full value 
pScrap.score100 <- ifelse(pScrap.score > 100,100,pScrap.score) # 100% ceiling

# 6) Percent Chironomidae
pChiron.score <- 100*(100-pChiron)/(100-bsv$BSV[bsv$metric=='pChiron']) # full value 
pChiron.score100 <- ifelse(pChiron.score > 100,100,pChiron.score)       # 100% ceiling

# 7) Percent 2 Dominant Taxa
p2dom.score <- 100*(100-p2dom)/(100-bsv$BSV[bsv$metric=='p2dom']) # full value 
p2dom.score100 <- ifelse(p2dom.score > 100,100,p2dom.score)       # 100% ceiling

# 8) HBI Hilsenhoff biotix index - (taxon abund x taxon tol val) / tot abund
hbi.score <- 100*(10-hbi)/(10-bsv$BSV[bsv$metric=='HBI']) # full value
hbi.score100 <- ifelse(hbi.score > 100,100,hbi.score)     # 100% ceiling

# Calc VASCI score --------------------------------------------------------------------
vasci200 <- round((tot.rich.score100 + EPT.rich.score100 + pE.score100 + pPT.H.score100  + pChiron.score100 + p2dom.score100 + pScrap.score100 + hbi.score100 )/8,3)
#metrics$vasci200 <- vasci200
vasci200

write.table(vasci200, "vasci.csv", sep = ",")


######### Plotting VASCI ##################
metrics <- read.csv("metrics.f21-s22.csv") 
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


shapiro.test(metrics$VASCI)
# greater than 0.05 = normal
library(car)
leveneTest(VASCI~Stream, data=metrics)
#less than 0.05 means non-parametric (heterogeneity of variance)

# Spring Kruskall
kruskal.test(VASCI~Stream, data=metrics.s)
# less than 0.05 means there are differences among groups

library(tidyr)
library(FSA)

Result = dunnTest(VASCI~Stream, data=metrics.s,
                  method="bonferroni")$res

### Use multcompView
library(multcompView)
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

box.s <- metrics.s %>% 
  ggplot(aes(x = group, y = VASCI, color = Impact)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) +
  geom_point(size = 3) +
  geom_point(data = longterm, aes(Stream, VASCI), color = "black") + 
  theme_classic() + #+ #Clean, minimal theme courtesy of the "egg" package
xlab("Stream (Spring)") +
  ylab("VASCI") +
  #stat_compare_means(method = "kruskal") +
  geom_text(data = letters.df.s, 
            aes(x = Stream, y = Placement.Value, label = Letters), 
            color = "black", hjust = -1.25,
            vjust = -0.8, fontface = "bold")
box.s

png("vascis.png")
plot(box.s)
dev.off()
# less than 0.05 means they are different

####Fall Kruskall
Result = dunnTest(VASCI~Stream, data=metrics.f,
                  method="bonferroni")$res

### Use multcompView
library(multcompView)
X = Result$P.adj <= 0.05
names(X) = gsub(" ",  "", Result$Comparison)
X <- multcompLetters(X)
X <- as.data.frame(X$Letters) 
colnames(X) <- c("Letters")

placement <- metrics.f %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(VASCI)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df.s <- cbind(X, placement)
unique(letters.df.s)

box.f <- metrics.f %>% 
  ggplot(aes(x = group, y = VASCI, color = Impact)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) +
  geom_point(size = 3) +
  geom_point(data = longterm, aes(Stream, VASCI), color = "black") + 
  theme_classic() + #+ #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Fall)") +
  ylab("VASCI") +
  #stat_compare_means(method = "kruskal") +
  geom_text(data = letters.df.s, 
            aes(x = Stream, y = Placement.Value, label = Letters), 
            color = "black", hjust = -1.25,
            vjust = -0.8, fontface = "bold")
box.f

png("vascif.png")
plot(box.f)
dev.off()
