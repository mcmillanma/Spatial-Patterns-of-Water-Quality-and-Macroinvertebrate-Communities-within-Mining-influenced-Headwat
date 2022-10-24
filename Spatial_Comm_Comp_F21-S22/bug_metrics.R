# ***** BIOLOGICAL METRICS *****
# Anthony J. Timpano    March 19, 2018
# McMillan Melanie edited Jan 20, 2022
# needs Kelly and Jeff results added to wq
# tributaries removed

#clear console ctrl+L
#clear workspace
rm(list=ls())
library(tidyverse)
library(vegan)
library(dplyr)
library(ggpubr)
library(corrplot)
library(patchwork)
library(ggplot2)


# Input:        'bugs.temp'  Data frame; Bug data matrix, raw abundance only (sample rows x taxon columns) \
bugs <- read_csv("bugid.csv")
taxa <- read_csv("taxa.list.csv")
gen.func <- read_csv("taxa.traits.csv")
chem <- NULL
chem <- read_csv("chem.f21-s22.notrib.csv")
chem <- chem %>%
  mutate("ca..mg" = ca.mgl / mg.mgl) %>%
  mutate ("so4..hco3" = so4.mgl /hco3.mgl)
#STREAM <- select(bugs, Site, Stream)
#chem <- left_join(chem, STREAM, by = c( "site.id" = "Site"))


# Computations: Multiple biological community metrics
# Output:       'metrics.temp'  Data frame; need to merge with sample info from full data frame

end_col <- ncol(bugs)
bugs.temp<- (bugs[,c(3:end_col)])
genera <- names(bugs.temp) # List of taxa found in bugs.temp
metrics.temp <- NULL  # Initialize metrics.temp data frame

# Total No. Individuals per sample for proportion calculations ----------------------
totind <- rowSums(bugs.temp)
metrics.temp <- data.frame(totind)

#abundance
#abund<-rowSums(bugs[,c(7:end_col)])
#metrics.temp$abund<- abund

#shannon
Hshannon <- diversity(bugs.temp,index="shannon")
metrics.temp$Hshannon<-Hshannon 

#Hsimpson
Hsimpson<- diversity(bugs.temp,index="simpson")
metrics.temp$Hsimpson<-Hsimpson

# GENUS RICHNESS of ALL TAXA ---------------------------------------------------------
rich <- rowSums(ifelse(bugs.temp > 0,1,0))
metrics.temp$rich <- rich

#evenness
J <- Hshannon/(log(rich))
metrics.temp$J<-J

# GENUS RICHNESS by ORDER ------------------------------------------------------------
# Determine which taxa in samples belong to each order by referencing Taxa List
E <- genera[which(genera %in% taxa$taxon[taxa$order=="Ephemeroptera"])]
P <- genera[which(genera %in% taxa$taxon[taxa$order=="Plecoptera"])]
Tr <- genera[which(genera %in% taxa$taxon[taxa$order=="Trichoptera"])]
C <- genera[which(genera %in% taxa$taxon[taxa$order=="Coleoptera"])]
D <- genera[which(genera %in% taxa$taxon[taxa$order=="Diptera"])]
# If genus in sample, make it 1, sum by row for each ORDER
rich.E <- rowSums(ifelse(bugs.temp[E] > 0,1,0))
rich.P <- rowSums(ifelse(bugs.temp[P] > 0,1,0))
rich.T <- rowSums(ifelse(bugs.temp[Tr] > 0,1,0))
rich.C <- rowSums(ifelse(bugs.temp[C] > 0,1,0))
rich.D <- rowSums(ifelse(bugs.temp[D] > 0,1,0))
metrics.temp$rich.E <- rich.E
metrics.temp$rich.P <- rich.P
metrics.temp$rich.T <- rich.T
metrics.temp$rich.C <- rich.C
metrics.temp$rich.D <- rich.D

# GENUS RICHNESS excluding Ephemeroptera ------
#not.E <- genera[which(genera %in% taxa$taxon[taxa$order!="Ephemeroptera"])]
#rich.less.E <- rowSums(ifelse(bugs.temp[not.E] > 0,1,0))
rich.less.E <- rich-rich.E
metrics.temp$rich.less.E <- rich.less.E

# EPT GENUS RICHNESS by ORDER --------------------------------------------------------
rich.EPT <- rich.E + rich.P + rich.T
metrics.temp$rich.EPT <- rich.EPT


# Richness/Percent Trichoptera less Hydropsychidae--------------------------------------
T.less.H <- genera[which(genera%in% taxa$taxon[taxa$order=="Trichoptera"&taxa$family!="Hydropsychidae"])]
rich.T.less.H <- rowSums(ifelse(bugs.temp[T.less.H]>0,1,0))
pT.less.H <- 100*rowSums(bugs.temp[T.less.H])/totind
metrics.temp$rich.T.less.H <- rich.T.less.H
metrics.temp$pT.less.H <- pT.less.H


#EPT less Hydropsychidae----------------------------------------
rich.EPT.less.H <- rich.E + rich.P + rich.T.less.H
metrics.temp$rich.EPT.less.H <- rich.EPT.less.H
# percent EPT less Hydropsychidae
mEPT <- c(E,P,T.less.H)
pEPT.less.H <- 100*rowSums(bugs.temp[mEPT])/totind
metrics.temp$pEPT.less.H <- pEPT.less.H

# Percent EPT --------------------------------
pEPT <- 100*rowSums(bugs.temp[c(E,P,Tr)])/totind
metrics.temp$pEPT <- pEPT

# Percent EPT less Cheumatopsyche -----------------------------------------
Tr.less.Cheum <- genera[which(genera %in% taxa$taxon[taxa$order=="Trichoptera" & taxa$taxon!="Cheumatopsyche"])]
rich.T.less.Cheum <- rowSums(ifelse(bugs.temp[Tr.less.Cheum] > 0,1,0))
mEPT <- c(E,P,Tr.less.Cheum)
pEPT.less.Cheum <- 100*rowSums(bugs.temp[mEPT])/totind
metrics.temp$pEPT.less.Cheum <- pEPT.less.Cheum

# Percent EPT less HBL ---------------------------------
E.less.B <- genera[which(genera %in% taxa$taxon[taxa$order=="Ephemeroptera" & taxa$family!="Baetidae"])]
P.less.L <- genera[which(genera %in% taxa$taxon[taxa$order=="Plecoptera" & taxa$family!="Leuctridae"])]
Tr.less.H <- genera[which(genera %in% taxa$taxon[taxa$order=="Trichoptera" & taxa$family!="Hydropsychidae"])]
mEPT <- c(E.less.B, P.less.L, Tr.less.H)
pEPT.less.HBL <- 100*rowSums(bugs.temp[mEPT])/totind
metrics.temp$pEPT.less.HBL <- pEPT.less.HBL

# PT GENUS richness -------
rich.PT <- rich.P + rich.T
metrics.temp$rich.PT <- rich.PT

# Percent Ephemeroptera --------------------------------------------------------------
pE <- 100*rowSums(bugs.temp[E],na.rm=TRUE)/totind
metrics.temp$pE <- round(pE,1)

# Richness/Percent Ephemeroptera less Baetidae -----------------------------------------
E.less.B <- genera[which(genera %in% taxa$taxon[taxa$order=="Ephemeroptera" & taxa$family!="Baetidae"])]
rich.E.less.B <- rowSums(ifelse(bugs.temp[E.less.B] > 0,1,0))
pE.less.B <- 100*rowSums(bugs.temp[E.less.B])/totind
metrics.temp$rich.E.less.B <- rich.E.less.B
metrics.temp$pE.less.B <- round(pE.less.B,1)

# Percent Plecoptera -----------------------------------------------------------------
pP <- 100*rowSums(bugs.temp[P],na.rm=TRUE)/totind
metrics.temp$pP <- round(pP,1)

# Percent Plecoptera less Leuctra-----------------------------------------------------------------
P.less.Leuc <- P[P!='Leuctra']
pP.less.Leuc <- 100*rowSums(bugs.temp[P.less.Leuc],na.rm=TRUE)/totind
metrics.temp$pP.less.Leuc <- round(pP.less.Leuc,1)

# Percent Plecoptera less Amphinemura-----------------------------------------------------------------
P.less.Amph <- P[P!='Amphinemura']
pP.less.Amph <- 100*rowSums(bugs.temp[P.less.Amph],na.rm=TRUE)/totind
metrics.temp$pP.less.Amph <- round(pP.less.Amph,1)

# Percent Plecoptera less Allocapnia-----------------------------------------------------------------
P.less.Allo <- P[P!='Allocapnia']
pP.less.Allo <- 100*rowSums(bugs.temp[P.less.Allo],na.rm=TRUE)/totind
metrics.temp$pP.less.Allo <- round(pP.less.Allo,1)

# Percent Plecoptera less Leuctra, Amphinemura, & Allocapnia -----------------------------------------------------------------
P.less.LAA <- P[!(P %in% c('Leuctra','Amphinemura','Allocapnia'))]
pP.less.LAA <- 100*rowSums(bugs.temp[P.less.LAA],na.rm=TRUE)/totind
metrics.temp$pP.less.LAA <- round(pP.less.LAA,1)

# Percent Trichoptera ----------------------------------------------------------------
pT <- 100*rowSums(bugs.temp[Tr],na.rm=TRUE)/totind
metrics.temp$pT <- round(pT,1)

# Percent PT less Hydropsychidae ------------
PT <- c(P,Tr)
PT.H <- PT[-which(PT %in% taxa$taxon[taxa$family=="Hydropsychidae"])]
pPT.H <- 100*rowSums(bugs.temp[PT.H])/totind
metrics.temp$pPT.H <- round(pPT.H,1)

# Percent CHIRONOMIDAE --------------------
pChi <- 100*bugs.temp$Chironomidae/totind
metrics.temp$pChi <- round(pChi,1)

# Percent DIPTERA---------------------------
pD <- 100*rowSums(bugs.temp[D], na.rm = TRUE)/totind
metrics.temp$pD <- round(pD,1)

# Percent OLIGOCHAETA 
pOligo <- 100*bugs.temp$Oligochaeta/totind
metrics.temp$pOligo <- round(pOligo,1)

# Percent CHIRONOMIDAE + OLIGOCHAETA --------------------
ChiO <- c('Chironomidae', 'Oligochaeta')
ChiO.found <- ChiO[ChiO %in% genera]
pChiO <- 100*(rowSums(bugs.temp[ChiO.found]))/totind
metrics.temp$pChiO <- pChiO

# Percent 1 Dominant Taxa------------------------------------------------------------
dom1 <- NULL
for(i in 1:length(bugs.temp[,1])){
  a <- sort(bugs.temp[i,1:length(bugs.temp)],decreasing = TRUE)
  dom1[i] <- sum(a[1])
}
p1dom <- 100*dom1/totind
metrics.temp$p1dom <- round(p1dom,1)

# Percent 2 Dominant Taxa ------------------------------------------------------------
dom2 <- NULL
for(i in 1:length(bugs.temp[,1])){
  a <- sort(bugs.temp[i,1:length(bugs.temp)],decreasing=TRUE)
  dom2[i] <- sum(a[1:2])
}
p2dom <- 100*dom2/totind
metrics.temp$p2dom <- round(p2dom,1)

# Percent 5 Dominant Taxa ------------------------------------------------------------
dom5 <- NULL
for(i in 1:length(bugs.temp[,1])){
  a <- sort(bugs.temp[i,1:length(bugs.temp)],decreasing=TRUE)
  dom5[i] <- sum(a[1:5])
}
p5dom <- 100*dom5/totind
metrics.temp$p5dom <- round(p5dom,1)

# HBI - (taxon abund x taxon tol val) / tot abund ------------------
#hbi <- NULL
#for(i in 1:length(bugs.temp[,1])){                    # one iteration per sample
# h <- 0                                             # initialize h to store hbi contribution of each taxon in sample
#  for(j in 1:length(bugs.temp)){                          # one iteration per taxon
#    abund <- bugs.temp[i,j]                        # abundance of taxon j in sample i
#    tol <- gen.func$tv[gen.func$taxon==genera[j]]   # tol value of taxon j
#    h <- h + abund*tol/totind[i]           # add hbi contribution of taxon to h
#  }
#  hbi[i] <- h                                        # store total hbi for sample i in ith spot of hbi vector
#}
#metrics.temp$hbi <- round(hbi,2)

# FFG metrics ----------------
shredders <- genera[which(genera %in% gen.func$taxon[gen.func$ffg=='SH'])]
rich.SH <- rowSums(ifelse(bugs.temp[shredders] > 0,1,0))
pSH <- 100*rowSums(bugs.temp[shredders])/totind

gatherers <- genera[which(genera %in% gen.func$taxon[gen.func$ffg=='CG'])]
rich.CG <- rowSums(ifelse(bugs.temp[gatherers] > 0,1,0))
pCG <- 100*rowSums(bugs.temp[gatherers])/totind

scrapers <- genera[which(genera %in% gen.func$taxon[gen.func$ffg=='SC'])]
rich.SC <- rowSums(ifelse(bugs.temp[scrapers] > 0,1,0))
pSC <- 100*rowSums(bugs.temp[scrapers])/totind

scrapers.less.E <- scrapers[which(scrapers %in% taxa$taxon[taxa$order != 'Ephemeroptera'])]
rich.SC.less.E <- rowSums(ifelse(bugs.temp[scrapers.less.E] > 0,1,0))
pSC.less.E <- 100*rowSums(bugs.temp[scrapers.less.E])/totind

filterers <- genera[which(genera %in% gen.func$taxon[gen.func$ffg=='CF'])]
rich.CF <- rowSums(ifelse(bugs.temp[filterers] > 0,1,0))
pCF <- 100*rowSums(bugs.temp[filterers])/totind

predators <- genera[which(genera %in% gen.func$taxon[gen.func$ffg=='PR'])]
rich.PR <- rowSums(ifelse(bugs.temp[predators] > 0,1,0))
pPR <- 100*rowSums(bugs.temp[predators])/totind

metrics.temp$rich.SH <- rich.SH
metrics.temp$rich.CG <- rich.CG
metrics.temp$rich.SC <- rich.SC
metrics.temp$rich.SC.less.E <- rich.SC.less.E
metrics.temp$rich.CF <- rich.CF
metrics.temp$rich.PR <- rich.PR
metrics.temp$pSH <- round(pSH,1)
metrics.temp$pCG <- round(pCG,1)
metrics.temp$pSC <- round(pSC,1)
metrics.temp$pSC.less.E <- round(pSC.less.E,1)
metrics.temp$pCF <- round(pCF,1)
metrics.temp$pPR <- round(pPR,1)

# HABIT metrics ----------------
clingers <- genera[which(genera %in% gen.func$taxon[gen.func$habit=='CN'])]
rich.Cling <- rowSums(ifelse(bugs.temp[clingers] > 0,1,0))
pCling <- 100*rowSums(bugs.temp[clingers])/totind
metrics.temp$rich.Cling <- rich.Cling
metrics.temp$pCling <- round(pCling,1)

burrowers <- genera[which(genera %in% gen.func$taxon[gen.func$habit== 'BU'])]
rich.Burrow <- rowSums(ifelse(bugs.temp[burrowers] > 0,1,0))
pBurrow <- 100*rowSums(bugs.temp[burrowers])/totind
metrics.temp$rich.Burrow <- rich.Burrow
metrics.temp$pBurrow <- round(pBurrow,1)

climbers <- genera[which(genera %in% gen.func$taxon[gen.func$habit=='CB'])]
rich.Climb <- rowSums(ifelse(bugs.temp[climbers] > 0,1,0))
pClimb <- 100*rowSums(bugs.temp[climbers])/totind
metrics.temp$rich.Climb <- rich.Climb
metrics.temp$pClimb <- round(pClimb,1)

sprawlers <- genera[which(genera %in% gen.func$taxon[gen.func$habit=='SP'])]
rich.Sprawl <- rowSums(ifelse(bugs.temp[sprawlers] > 0,1,0))
pSprawl <- 100*rowSums(bugs.temp[sprawlers])/totind
metrics.temp$rich.Sprawl <- rich.Sprawl
metrics.temp$pSprawl <- round(pSprawl,1)

swimmers <- genera[which(genera %in% gen.func$taxon[gen.func$habit=='SW'])]
rich.Swimm <- rowSums(ifelse(bugs.temp[swimmers] > 0,1,0))
pSwimm <- 100*rowSums(bugs.temp[swimmers])/totind
metrics.temp$rich.Swimm <- rich.Swimm
metrics.temp$pSwimm <- round(pSwimm,1)

# Diversity metrics -------------------------------------
#div.mat <- bugs.temp
#div.shan <- diversity(div.mat, index='shannon')
#div.simp <- diversity(div.mat, index='simpson')
#metrics.temp$div.shan <- round(div.shan,3)
#metrics.temp$div.simp <- round(div.simp,3)


# Evenness metrics -----
# Pielou
#even.pielou <- div.shan / log(rich)
#metrics.temp$even.pielou <- round(even.pielou,4)

# INTOLERANT Richness -----------------
INT.taxa <- gen.func$taxon[gen.func$tv <= 3]
INT <- genera[genera %in% INT.taxa]
rich.INT <- rowSums(ifelse(bugs.temp[INT] > 0,1,0))
metrics.temp$rich.INT <- rich.INT

# TOLERANT Richness 
TOL.taxa <- gen.func$taxon[gen.func$tv >= 7]
TOL <- genera[genera %in% TOL.taxa]
rich.TOL <- rowSums(ifelse(bugs.temp[TOL] > 0,1,0))
metrics.temp$rich.TOL <- rich.TOL

# PERCENT TOLERANT ----------------------------
TOL.taxa <- gen.func$taxon[gen.func$tv >= 7]
TOL <- genera[genera %in% TOL.taxa]
pTOL <- 100*rowSums(bugs.temp[TOL])/totind
metrics.temp$pTOL <- round(pTOL,digits=5)

# CLEAN UP -----
#rm(list=names(metrics.temp))
#rm(list=c('E','P','Tr','C','D','E.less.B','PT','PT.H','P.less.L','P.less.Allo','P.less.Amph','P.less.Leuc','P.less.LAA','h','i','j','tol','abund','clingers','dom2','dom5','filterers','gatherers','predators','scrapers','scrapers.less.E','shredders','Tr.less.Cheum','Tr.less.H','INT','INT.taxa','ChiO','ChiO.found','mEPT'))

# remove duplicates samples and sites that have been discontinued (i.e., CLE, COA, DAV)
#badsites<-c("RIC-DUP","DAV","KEL-DUP", "HCN-DUP","COA", "CLE")
#metrics <- metrics[ ! metrics$Site %in% badsites, ]
#remove tributaries
tribs <- c("FRY-T3", "CRO-T1", "CRO-T2", "LLW-T1", "LLW-T2", "SPC-T2", "SPC-T3", "SPC-T4","ROL-T1", "ROL-T2", "LLW-T3", "FRY-T1", 
           "FRY-T2", "EAS-T1", "CRO-T3")
metrics.temp$Site<-bugs$Site
metrics.temp$Season<-bugs$Season
metrics.temp$Stream_Type <- bugs$type




#metrics.temp <- metrics.temp[ ! metrics.temp$Site %in% tribs, ] 
#left_join(chem, by = c("Site" = "site.id") )
  

#save metrics as .csv
write.csv(metrics.temp, file="shannon.f21-s22.csv", sep = ",")

# Experimenting with ANOVA
#install.packages("ggpubr")


