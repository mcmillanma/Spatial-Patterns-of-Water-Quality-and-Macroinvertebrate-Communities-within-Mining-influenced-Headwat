#PCoA
library(ape) # load the ape package with a nifty PCoA function
library(vegan) # load the vegan package for vegdist
library(dplyr)
library(tidyr)
library(ggPCoA)

#Spring

bugs <- read.csv("bugid.csv") 
bugs <- dplyr::filter(bugs, Season == "Spring")
bugs <- bugs[, colSums(bugs != 0) > 0]
bugs <- bugs[, 4:93]

Design <- read.csv("metrics.reduced.csv") %>%
  filter(Season == "Spring")
Design <- Design[, colSums(Design != 0) > 0]

bugs.dist <- vegdist(log(bugs[,-1]+1), method='bray', binary=FALSE) # producing a distance matrix using vegdist, notice I applied a log transformation in there like in Wepking et al. ####
bugs.dist = as.data.frame(as.table(bugs.dist))
histogram(bugs.dist$Freq)

#bugs.dist1 <- vegdist(bugs, method='bray', binary=FALSE)
#bugs.dist1 = as.data.frame(as.table(bugs.dist1))
#histogram(bugs.dist1$Freq)

bugs.dist <- vegdist(log(bugs[,-1]+1), method='bray', binary=FALSE)
bugs.pco <- pcoa(bugs.dist, correction='none') # running the PCoA ####
summary(bugs.pco)
bugs.pco
bugs.pco$values
bugs.pco$vectors


Z <- data.frame(Design, bugs.pco$vectors[,1:2]) # merging the Design information with the PCoA output for the 1st 2 axes

plot <- Z %>%
  ggplot( aes(x=Axis.1, y=Axis.2, 
           color= 'Stream', size= VASCI))+
geom_point(alpha=0.5,aes(color = Stream))+ 
theme_classic() +
scale_shape_manual( name="Impact") + 
  #scale_size_continuous(name="VASCI") + 
  #scale_color_discrete(name="Stream") + 
  stat_ellipse(aes(color = Stream))

plot




