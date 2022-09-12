# OBJECTIVE 1: Determine with-in stream variability of AMI communities
# PART 1: Is it changing? A) Simple regression for each stream of community similarity vs 
# distance downstream. 
# PART 2: What is changing? Use step-wise multiple regression to determine which bug 
# metrics and/or taxa are (or are not) shifting in each stream.

### PART 1
# CREATE INSECT DISSIMILARITY INDEX
# Bryan Brown (VT) quantitative analysis script 2021
# Compute similarites by looping to produce independent ID variables and thus 
# ensure that we get correct correspondence between Distance and Similarity

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
      sim <- vegdist(ab, method= 'bray', binary=FALSE)  # the similarity function
      SIM <- 1-(as.numeric(sim))    # converting from default dissimilarity to similarity
      out <- data.frame(START=Site.list[i], END=Site.list2[j], Season=Metricz$Season[Metricz$Site==Site.list[i]], Stream=Stream.list[k], SIM)
      Simile <- rbind(Simile,out)  
    }
  }
  Similar <- rbind(Similar, Simile) 
}

#Similar <- pivot_wider(Similar,names_from = END, values_from = SIM) %>%
#select(-c(Season, Stream))

#write.csv(Similar, file="Similarity.f21.csv", sep = ",")

#PLOT INSECT DISSIMILARITY INDEX VS EUCLIDIAN DISTANCE

dis.simm <- read_csv("dis.simm.f21.csv") 
Similar$dis.simm <- dis.simm$dis.simm
Similar <- filter(Similar, SIM != 1) 

Similar

bray.spearman.dist.dxsimm.stream <- Similar %>% 
  filter(START %in% c("CRO-1", "EAS-1", "FRY-1", "LLW-1", "ROL-1", "SPC-1")) %>%
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
plot(Similar$dis.simm[Similar$Stream == "FRY (MN-G)"], Similar$SIM[Similar$Stream== "FRY (MN-G)"], xlab='Distance', ylab='', main='Mined Non-gradient')
reg.2 <- lm(SIM[Similar$Stream == "FRY (MN-G)"]~dis.simm[Similar$Stream== "FRY (MN-G)"], data=Similar)
abline(reg.2$coefficients)
plot(Similar$dis.simm[Similar$Stream %in% c("LLW (MG)", "ROL (MG)", "SPC (MG)")], Similar$SIM[Similar$Stream %in% c("LLW (MG)", "ROL (MG)", "SPC (MG)")], xlab='Distance', ylab='Community Similarity', main='Mined Gradient')
reg.3 <- lm(SIM[Similar$Stream %in% c("LLW (MG)", "ROL (MG)", "SPC (MG)")]~dis.simm[Similar$Stream %in% c("LLW (MG)", "ROL (MG)", "SPC (MG)")], data=Similar)
abline(reg.3$coefficients)

#plotting similarity by stream
Similar <- read_csv("Similarity.f21.csv")

bray.spearman.distxsimm.stream <- Similar %>% ggplot(aes(x=dis.simm ,y=SIM, group = Stream_Type))+
  geom_point()  +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman",
           label.y = 0.64) +
  stat_regline_equation(label.y = 0.6) +
  theme_classic(base_size = 15)+
  xlab("Distance Between Samples (m)") +
  ylab("Simmilarity") +
  ggtitle("Distance vs Similarity")+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  facet_wrap("Stream_Type") 

bray.spearman.distxsimm.stream

#set size and save plot as png
png("bray.spearman.distxsimm.stream.png", width = 900, height = 450)
plot(bray.spearman.distxsimm.stream)
dev.off()

### PART 2
### CREATE NMDS with vectors of key metrics and taxa

library(tidyverse)
library(vegan)
library(ggplot2)

bugs <- read_csv("bugID.f21.csv")
tribs <- c("FRY-T3", "CRO-T1", "CRO-T2", "LLW-T1", "LLW-T2", "SPC-T2", "SPC-T3", "SPC-T4","ROL-T1", "ROL-T2", "LLW-T3", "FRY-T1", 
           "FRY-T2", "EAS-T1", "CRO-T3")
bugs <- bugs[ ! bugs$Site %in% tribs, ]
metrics <- metrics.temp
NMDS <- "null"
NMDS <- left_join(metrics, bugs, by = "Site")
#habitatmaster from Habitat.R

com = NMDS[,76:224]
env = NMDS[,3:66] #habitat and water chem

#convert com to a matrix
m_com = as.matrix(com)

#nmds code
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
scores <- scores(nmds)


en = envfit(nmds, env, p.max = 0.001, permutations = 999, na.rm = TRUE)
scores(en, "vectors")
en

plot(nmds)
NMDS_metrics_graph <- plot(en, p.max = 0.001)
png("NMDS_abund_metrics.png", width = 950, height = 450)
plot(en, p.max = 0.001)
dev.off()

#n add columns from your original data that contain information that you would like to include in your plot
data.scores = as.data.frame(scores(nmds, "sites"))
data.scores$Stream = NMDS$Stream.x

#Because my data contained only continuous environmental variables, I’m extracting the information from both separately using the “vectors” and not “factors” options.
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)


#plot it
gg <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Stream), size = 3, alpha = 0.5) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") + 
  geom_point(data = en_coord_cont, aes(x = NMDS1, y = NMDS2)) +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont), nudge_y = 0.1) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Stream Type") +
  ggtitle("Water Quality Controls on Community Composition") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))

gg

#set size and save plot as png
png("NMDS_abund_metrics.png", width = 950, height = 450)
plot(gg)
dev.off()


## source: https://www.youtube.com/watch?v=QljEeBei-JA; 
# source: https://rstatisticsandresearch.weebly.com/uploads/1/0/2/6/1026585/nmds_updatedscript.r
# subset the dataframe on which to base the ordination (dataframe 1)
data_1 <- bugs[,8:156]

#Identify the columns that contains the descriptive/environmental data (dataframe 2)
data_2 <- NMDS[,c(3:66, 67,68)]


#ordination by NMDS
NMDS_metrics <- metaMDS(data_1, distance = "bray", k = 2)

#########################
#Data visualisation (THIS IS AN UPDATED VERSION OF THE SCRIPT, NOW USING GGPLOT)

#Extract the axes scores
datascores <- as.data.frame(scores(NMDS_metrics, "sites"))  #extract the site scores

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Stream = data_2$Stream)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Stream, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Stream','oNMDS1','oNMDS2')),
             by = 'Stream', sort = FALSE)

#plot
ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Stream)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # add spiders
  geom_point(data = centroids, size = 4) +                    # add centroids
  geom_point() +                                              
  coord_fixed()+                                              
  theme_bw()+ 
  theme(legend.position="right",legend.text=element_text(size=10),legend.direction='vertical')

#####################
#Bootstrapping and testing for differences between the groups
fit <- adonis(data_1 ~ Stream, data=data_2, permutations=999, method="bray")
fit
?be#####################
#Check assumption of homogeneity of multivariate dispersion
distances_data <- vegdist(data_1)
anova(betadisper(distances_data, data_2$Stream))