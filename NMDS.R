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



#source: https://jkzorz.github.io/2020/04/04/NMDS-extras.html

bugs <- read_csv("bugid.csv")
bugs.s <- bugs %>%
  filter(Season == "Spring") %>%
  filter(Site != "EAS1" & Site != "EAS9" & Site !="LLW6")
bugs.f <- bugs %>%
  filter(Season == "Fall") %>%
  filter(Site != "EAS1" & Site != "EAS9" & Site !="LLW6")
  
chem <- read.csv("chem.f21-s22.notrib.csv")  
  
# fall: sc, ca, hardness, so4, mg, cu, so4/hcos, k, hco3, u
# spring: sc, ca, hardness, so4, mg, cu, so4/hcos, k, hco3, u + ba, na, ph
chem.s <- chem %>%
  filter(season == "Spring") %>%
  select( c(1, 5, 7, 9, 12, 20, 22, 23, 24, 25, 26, 36, 45, 47, 49 ))
chem.f <- chem %>%
  filter(season == "Fall") %>%
  filter(Site != "SPC6") %>%
  select( c(1, 5, 9, 12, 20, 22, 24, 25, 26, 36, 47, 49 ))
habitat <- read.csv("habitatmaster.csv") %>% # needs to match length of bugs
  filter(Site != "EAS1" & Site != "EAS9" & Site !="LLW6" ) %>% #& Site !="FRY8"& Site !="FRY9"& Site !="ROL7"& Site !="SPC6") %>%
  select( c(1, 7:9, 11,12, 16, 17, 20, 21, 22, 23, 29, 32))   # avg.slope, avgriparianwidR, avgriparianwidL, avgvegprotecR, D50, avgwetwidth, pfines, psmallcobble, ppebbles, avgembedd


#habitatmaster from Habitat.R

#subset your data so that you have a data frame with only environmental variables (env), and a data frame with only species abundance data (com)
com = bugs.s[,3:119]
env = habitat[,2:14] #habitat and water chem

#convert com to a matrix
m_com = as.matrix(com)

#nmds code
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
scores <- scores(nmds)


en = envfit(nmds, env, permutations = 999, na.rm = TRUE)
en

plot(nmds)
plot(en, p.max = 0.05)

#n add columns from your original data that contain information that you would like to include in your plot
data.scores = as.data.frame(scores(nmds)$sites)
data.scores$Stream = habitat$Stream


#Because my data contained only continuous environmental variables, I’m extracting the information from both separately using the “vectors” and not “factors” options.
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)


#plot it
nmds.s <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
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
  labs(colour = "Stream") +
  ggtitle("Environmental Controls on Community Composition (Spring 2022)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))

nmds.s
nmds.f

nmds <- ggarrange(nmds.f, nmds.s, ncol=1)
nmds
#set size and save plot as png
png("NMDS_abund_hab_season.png", width = 800, height = 1000)
plot(nmds)
dev.off()


## source: https://www.youtube.com/watch?v=QljEeBei-JA; 
# source: https://rstatisticsandresearch.weebly.com/uploads/1/0/2/6/1026585/nmds_updatedscript.r
# subset the dataframe on which to base the ordination (dataframe 1)
data_1 <- bugs[,7:155]

#Identify the columns that contains the descriptive/environmental data (dataframe 2)
data_2 <- bugs[,4:5]


#ordination by NMDS
NMDS <- metaMDS(data_1, distance = "bray", k = 2)

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

#####################
#Bootstrapping and testing for differences between the groups
fit <- adonis(data_1 ~ Stream, data=data_2, permutations=999, method="bray")
fit
?be#####################
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

#########################
#Data visualisation (THIS IS AN UPDATED VERSION OF THE SCRIPT, NOW USING GGPLOT)
install.packages("ggplot2")
library("ggplot2")

#Extract the axes scores
datascores <- as.data.frame(scores(NMDS))  #extract the site scores

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
#as.numeric(unlist(x)) 

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
