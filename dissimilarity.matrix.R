Bray-Curtis and compare with Jaccard (presence absence)
library(vegan)
library(tidyverse)
library(ggpubr)
## Example: 
#data(varespec) 
#vare.dist <- vegdist(varespec)
#vare.dist %>%
  #as.matrix() %>%
  #as.data.frame() 
  #tibble::rownames_to_column() 
  #pivot_longer(-rowname)

# Bryan Brown Script
# Compute similarites by looping to produce independent ID variables and thus ensure that 
# we get correct correspondence between Distance and Similarity
bugs <- read.csv("metrics.f21-s22.csv") %>%
  filter( Season == "Spring")

Similar <- c()
Stream.list <- unique(bugs$Stream)
for(k in 1:length(Stream.list)){
  Metricz <- bugs[bugs$Stream==Stream.list[k],]
  Site.list <- unique(Metricz$Site)                            
  Simile <- c()
  for(i in 1:length(Site.list)){
    a <- Metricz[Metricz$Site==Site.list[i],-c(1:6)]
    Site.list2 <- Site.list[(i-1):length(Site.list)]
    for(j in 1:length(Site.list2)){
      b <- Metricz[Metricz$Site==Site.list2[j],-c(1:6)]
      ab <- rbind(a,b)
      sim <- vegdist(ab, method= 'jaccard', binary=TRUE)  # the similarity function
      SIM <- 1-(as.numeric(sim))    # converting from default dissimilarity to similarity
      out <- data.frame(START=Site.list[i], END=Site.list2[j], Season=Metricz$Season[Metricz$Site==Site.list[i]], Stream=Stream.list[k], SIM)
      Simile <- rbind(Simile,out)  
    }
  }
  Similar <- rbind(Similar, Simile) 
}

#Similar <- pivot_wider(Similar,names_from = END, values_from = SIM) %>%
  #select(-c(Season, Stream))


#write.csv(similar, file="Similarity.csv", sep = ",")

#plotting it

Similar <- read.csv("Similarity.csv")

Similar

bray.spearman.dist.dxsimm.stream <- Similar %>% 
  filter(START %in% c("CRO1", "EAS1", "FRY1", "LLW1", "ROL1", "SPC1")) %>%
  ggplot(aes(x=dis.simm ,y=SIM, group = Stream))+
  geom_point(aes(color=Stream, pch = Stream)) +
  geom_line(aes(color=Stream)) +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman") +
  theme_classic()+
  xlab("Distance Downstream (m)") +
  ylab("Simmilarity to most downstream site") +
  ggtitle("Distance Downstream vs Similarity")+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  facet_wrap("Stream", scales = "free") +
  geom_vline(aes(xintercept = Distance), trib.loc.diluting, color = "green") +
  geom_vline(aes(xintercept = Distance), trib.loc.concentrating, color = "red") +
  geom_vline(aes(xintercept = Distance), trib.loc.noninfluencing, color = "yellow") 

  bray.spearman.dist.dxsimm.stream

#set size and save plot as png
png("bray.spearman.dist.dxsimm.stream.png", width = 900, height = 450)
plot(bray.spearman.dist.dxsimm.stream)
dev.off()

#Different way to plot all site relations
#plotting similarity by stream type
par(mfrow=c(1,3))
plot(Similar$dis.simm[Similar$Stream %in% c("EAS (R)", "CRO (R)")], Similar$SIM[Similar$Stream %in% c("EAS (R)", "CRO (R)")], xlab='Distance', ylab='Community Similarity', main='Reference')
reg.1 <- lm(SIM[Similar$Stream %in% c("EAS (R)", "CRO (R)")]~dis.simm[Similar$Stream %in% c("EAS (R)", "CRO (R)")], data=Similar)  # calculating a simple linear regression
abline(reg.1$coefficients)   # Using the regression to draw a slope on the figure
plot(Similar$dis.simm[Similar$Stream == "FRY (MI)"], Similar$SIM[Similar$Stream== "FRY (MI)"], xlab='Distance', ylab='', main='Mined Non-gradient')
reg.2 <- lm(SIM[Similar$Stream == "FRY (MI)"]~dis.simm[Similar$Stream== "FRY (MI)"], data=Similar)
abline(reg.2$coefficients)
plot(Similar$dis.simm[Similar$Stream %in% c("LLW (MI)", "ROL (MI)", "SPC (MI)")], Similar$SIM[Similar$Stream %in% c("LLW (MI)", "ROL (MI)", "SPC (MI)")], xlab='Distance', ylab='Community Similarity', main='Mined Gradient')
reg.3 <- lm(SIM[Similar$Stream %in% c("LLW (MI)", "ROL (MI)", "SPC (MI)")]~dis.simm[Similar$Stream %in% c("LLW (MI)", "ROL (MI)", "SPC (MI)")], data=Similar)
abline(reg.3$coefficients)

#plotting similarity by stream

Similar <- read_csv("Similarity.csv")

b.plot <- Similar %>% ggplot(aes(x=sc.delt , y=SIM.b, group = Season, color = Season))+
  geom_point()  +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman", label.y.npc = 0.3) +
  stat_regline_equation(label.x.npc = 0.5, label.y.npc = 0.15) +
  theme_classic()+
  xlab("Change in Specific Conductance (uS/cm)") +
  ylab("Bray-Curtis Simmilarity") +
  ggtitle("Bray-Curtis Community Simmilarity vs Distance Between Samples") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  #scale_y_continuous(limits = c(0.4, 1)) +
  facet_wrap("Stream")
j.plot
b.plot

sim <- ggarrange(b.plot, j.plot, ncol= 1)
sim

#set size and save plot as png
png("bray+jaccxsc.stream.season.png", width = 900, height = 1200)
plot(sim)
dev.off()

par(mfrow=c(2,3), mar = c(3,3,8,3)) 
plot(Similar$dis.simm[Similar$Stream =="EAS (R)"], Similar$SIM[Similar$Stream =="EAS (R)"], xlab='Distance', ylab='Community Similarity', main='EAS (R)')
reg.1 <- lm(SIM[Similar$Stream =="EAS (R)"]~dis.simm[Similar$Stream =="EAS (R)"], data=Similar)  # calculating a simple linear regression
abline(reg.1$coefficients)   # Using the regression to draw a slope on the figure
plot(Similar$dis.simm[Similar$Stream == "FRY (MI)"], Similar$SIM[Similar$Stream== "FRY (MI)"], xlab='Distance', ylab='', main='FRY (MI)')
reg.2 <- lm(SIM[Similar$Stream == "FRY (MI)"]~dis.simm[Similar$Stream== "FRY (MI)"], data=Similar)
abline(reg.2$coefficients)
plot(Similar$dis.simm[Similar$Stream == "CRO (R)"], Similar$SIM[Similar$Stream == "CRO (R)"], xlab='Distance', ylab='', main='CRO (R)')
reg.3 <- lm(SIM[Similar$Stream == "CRO (R)"]~dis.simm[Similar$Stream == "CRO (R)"], data=Similar)
abline(reg.3$coefficients)
plot(Similar$dis.simm[Similar$Stream == "LLW (MI)"], Similar$SIM[Similar$Stream == "LLW (MI)"], xlab='Distance', ylab='Community Similarity', main='LLW (MI)')
reg.4 <- lm(SIM[Similar$Stream == "LLW (MI)"]~dis.simm[Similar$Stream == "LLW (MI)"], data=Similar)
abline(reg.4$coefficients)
plot(Similar$dis.simm[Similar$Stream == "ROL (MI)"], Similar$SIM[Similar$Stream == "ROL (MI)"], xlab='Distance', ylab='', main='ROL (MI)')
reg.5 <- lm(SIM[Similar$Stream == "ROL (MI)"]~dis.simm[Similar$Stream == "ROL (MI)"], data=Similar)
abline(reg.5$coefficients)
plot(Similar$dis.simm[Similar$Stream == "SPC (MI)"], Similar$SIM[Similar$Stream == "SPC (MI)"], xlab='Distance', ylab='', main='SPC (MI)')
reg.6 <- lm(SIM[Similar$Stream == "SPC (MI)"]~dis.simm[Similar$Stream == "SPC (MI)"], data=Similar)
abline(reg.6$coefficients)
mtext("Community Similarity vs Distance Between Samples", side = 3, line = -1, outer = TRUE))


summary(reg.6)


(CRO +EAS +FRY) /(ROL + LLW+ SPC) +plot_annotation(title = "Aquatic Insect Dissimilarity vs Distance Downstream")


