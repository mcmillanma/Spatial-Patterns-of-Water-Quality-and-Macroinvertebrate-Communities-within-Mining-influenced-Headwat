# experimenting with Shupryt funstions
envfit()
metaMDS()
adonis()
install.packages(lme4) ?
  install.packages(ncf) ?
  correlog()

library(vegan)
library(ggplot2)
library(tidyverse)
library(ggpubr)

#### PREPARING DATA ######
bugs <- read_csv("bugid.csv") %>%
  filter(Season == "Spring") %>%
  filter(Stream == "CRO (R)")

#source: https://jkzorz.github.io/2020/04/04/NMDS-extras.html

bugs <- read_csv("bugid.csv")
bugs.s <- bugs %>%
  filter(Season == "Spring") %>%
  filter(Site != "EAS1" & Site != "EAS9" & Site !="LLW6")
bugs.f <- bugs %>%
  filter(Season == "Fall") %>%
  filter(Site != "EAS1" & Site != "EAS9" & Site !="LLW6")
  
chem <- read.csv("chem.f21-s22.notrib.reduced.csv")  
  
# fall: sc, ca, hardness, so4, mg, cu, so4/hcos, k, hco3, u
# spring: sc, ca, hardness, so4, mg, cu, so4/hcos, k, hco3, u + ba, na, ph
chem.s <- chem %>%
  filter(season == "Fall") %>%
  filter(Site != "SPC6") #%>%
  select( c(1, 5, 7, 9, 12, 20, 22, 23, 24, 25, 26, 36, 45, 47, 49 ))
chem.f <- chem %>%
  filter(season == "Fall") %>%
  filter(Site != "SPC6") %>%
  select( c(1, 5, 9, 12, 20, 22, 24, 25, 26, 36, 47, 49 ))
hab <- read.csv("habitatmasteradj.csv") #%>%# needs to match length of bugs
  filter(Site != "EAS1" & Site != "EAS9" & Site !="LLW6" ) %>% #& Site !="FRY8"& Site !="FRY9"& Site !="ROL7"& Site !="SPC6") %>%
  select( c(1, 7:9, 11,12, 16, 17, 20, 21, 22, 23, 29, 32))   # avg.slope, avgriparianwidR, avgriparianwidL, avgvegprotecR, D50, avgwetwidth, pfines, psmallcobble, ppebbles, avgembedd


## source: https://www.youtube.com/watch?v=QljEeBei-JA; 
# source: https://rstatisticsandresearch.weebly.com/uploads/1/0/2/6/1026585/nmds_updatedscript.r
# subset the dataframe on which to base the ordination (dataframe 1)

# axis chem
data_1 <- chem.s[, 11:15] # nutrients
data_1 <- chem.s[, 20:41] # trace elements
data_1 <- decostand(data_1, method = "standardize",  na.rm =TRUE)
data_1 <- data_1 %>%
  select_if(~ ! any(is.na(.)))
bugs.pca <- prcomp(data_1, center = TRUE,scale. = TRUE)
data.scores = as.data.frame(scores(bugs.pca))

chem.pca <- NULL
chem.pca <- chem.s[,1:9]
chem.pca$hardness.mgl <- chem.s[,10]
chem.pca$cl.mgl <- chem.s[,17]
chem.pca$npoc.mgl <- chem.s[,16]
chem.pca$nutrients.pca <- data.scores$PC1
chem.pca$elements.pca <- data.scores$PC1

write.csv(chem.pca, file = "chem.pca.fall.csv")

# axis for habitat
data_1 <- hab[, 13:22] # mean diameters/pebbles size
data_1 <- hab[, 26:31] # riparian characterization
data_1 <- decostand(data_1, method = "standardize",  na.rm =TRUE)
data_1 <- data_1 %>%
  select_if(~ ! any(is.na(.)))
bugs.pca <- prcomp(data_1, center = TRUE,scale. = TRUE)
data.scores = as.data.frame(scores(bugs.pca))

hab.pca <- NULL
hab.pca <- hab[,c(1:2 , 7:9 ,11:12 , 23:25, 32,33)]
hab.pca$D.pca <- data.scores$PC1
hab.pca$Rip.pca <- data.scores$PC1

write.csv(hab.pca, file = "hab.pca.csv")

#for bugs
data_1 <- bugs[,4:120]
data_1 <- data_1[, colSums(data_1 != 0, na.rm = TRUE) > 0]

#Identify the columns that contains the descriptive/environmental data (dataframe 2)
data_2 <- bugs[,1:3]

######PCA########
#ordination by PCA https://www.datacamp.com/tutorial/pca-analysis-r
bugs.pca <- prcomp(data_1, center = TRUE,scale. = TRUE)
summary(bugs.pca)
str(bugs.pca)

library(ggbiplot)
ggbiplot(bugs.pca, ellipse=TRUE,labels = bugs$Site, groups = bugs$Stream) + 
  theme_classic()

data.scores = as.data.frame(scores(bugs.pca))
env <- data_1 #%>%
  select(c(4:6 , 14, 35, 41, 49, 51, 65, 87, 99, 114, 116))
ef <- envfit(bugs.pca, env)
ef
vectors <- colnames(env)
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
vec.df$variables <- rownames(vec.df)
scores <- cbind(as.data.frame(data.scores), Stream = data_2$Stream)

xx = ggplot(scores, aes(x = PC1, y = PC2)) + 
  geom_point(size = 4, aes(shape = Stream, color = Stream)) +
  stat_ellipse(aes(color = Stream)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               colour="blue",
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = PC1, y=PC2, label = variables),
            size=4) +
  ggtitle(" PCA Spring ") +
theme_classic()

xx

png("SpeciesNMDSallstreams.png", width = 650, height = 550)
plot(xx)
dev.off()

#########NMDS############
nmds <- metaMDS(data_1, distance = "bray", k=2, trymax=1, 
                autotransform=TRUE, noshare=0.1, expand=TRUE, trace=1, plot=FALSE)
plot1 <- ordiplot(nmds, choices=c(1,2))
str(nmds)

# Note, I set trymax = 100 because at the default 20 it was not converging

# envfit. First make dataframe for environmental variables
env <- data_1 #%>%
  select(c(1,5,6,12,13,19,22,23,35,41,42))
ef <- envfit(nmds, env)
ef

vectors <- colnames(env)
# Note that salinity is NOT significantly correlated with community structure

# Default vegan plot with vector
plot(nmds)
plot(ef)

# Edit dataframe
datascores <- as.data.frame(scores(nmds)$sites)
nmds$stress
stressplot(nmds)

# Make dataframe with vector to add to ggplot
# See https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2 for more discussion
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
vec.df$variables <- rownames(vec.df)
scores <- cbind(as.data.frame(datascores), Site = data_2$Site)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Stream, data = scores, FUN = mean)

# ggplot

xx = ggplot(scores, aes(x = NMDS1, y = NMDS2, label = Site)) + 
  geom_point(size = 4) +
  #stat_ellipse(aes(color = Stream)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour="blue",
               inherit.aes = FALSE) + 
  geom_text(hjust=1, vjust=-1) +
  geom_text(data = vec.df,
            aes(x = NMDS1, y=NMDS2, label = variables),
            size=4, nudge_y = 0.05) +
  labs(x = "NMDS1", y = "NMDS2", shape = "Stream") +
  ggtitle("SPC (MI) Spring") +
  theme_classic()
  
  xx
  
  png("SpeciesNMDSFRY.png", width = 650, height = 550)
  plot(xx)
  dev.off()
#########################
#Data visualisation (THIS IS AN UPDATED VERSION OF THE SCRIPT, NOW USING GGPLOT)
#install.packages("ggplot2")
library("ggplot2")

#Extract the axes scores
datascores <- as.data.frame(scores(NMDS)$sites)  #extract the site scores

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

# ADD VECTORS

#####################
#Bootstrapping and testing for differences between the groups

fit <- adonis2(data_1 ~ Stream, data=data_2, permutations=999, method="bray")
fit
#####################
#Check assumption of homogeneity of multivariate dispersion

distances_data <- vegdist(data_1)
anova(betadisper(distances_data, data_2$Stream))

########################################
#subset the dataframe on which to base the ordination (dataframe 1)
chem_1 <- chem[,7:48]

#Identify the columns that contains the descriptive/environmental data (dataframe 2)
chem_2 <- chem[c(2, 49)]



#ordination by NMDS
NMDSchem <- metaMDS(chem_1, distance = "bray", k = 2, na.rm = TRUE)

#Extract the axes scores
chemscores <- as.data.frame(scores(NMDSchem))  #extract the site scores

#Add/calculate spider diagram
scores2 <- cbind(as.data.frame(chemscores), Stream = chem_2$Stream)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Stream, data = scores, FUN = mean)
seg <- merge(scores2, setNames(centroids, c('Stream','oNMDS1','oNMDS2')),
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
fit <- adonis(chem_1 ~ Stream, data=chem_2, permutations=999, method="bray")
fit
#####################
#Check assumption of homogeneity of multivariate dispersion
distances_data <- vegdist(chem_1)
anova(betadisper(distances_data, chem_2$Stream))


#PCA
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot) 


#Remove columns with only 0
metrics.temp <- read_csv("metrics.f21.csv")
metrics.temp <- metrics.temp[, colSums(metrics.temp != 0) > 0]
row.names(metrics.temp) <- metrics.temp$Site

#subset the dataframe on which to base the ordination (dataframe 1)
metrics_1 <- metrics.temp[, 1:64]

#Identify the columns that contains the descriptive/environmental data (dataframe 2)
metrics_2 <- metrics.temp[c(65:66)]

bug.pca <- prcomp(metrics_1, center = TRUE,scale. = TRUE)
# source: https://www.datacamp.com/community/tutorials/pca-analysis-r
ggbiplot(bug.pca, labels = rownames(metrics.temp))

#How to group by stream?
#metrics.streams <- c(Stream("CRO (R)", Stream("EAS (R)"), Stream("FRY (MN-G)"), Stream("SPC (MG)"), 

#plotting
plot(bug.pca$x[,1], bug.pca$x[,2])

## make a scree plot: https://github.com/StatQuest/pca_demo/blob/master/pca_demo.R
pca.var <- bug.pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

#Same thing in ggplot
pca.data <- data.frame(Sample = rownames(metrics.temp),
                       X=bug.pca$x[,1],
                       Y=bug.pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, group = )) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Aquatic Insect Metrics")
