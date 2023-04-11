# Final Figures

### NMDS ###
library(BiodiversityR) # also loads vegan
library(ggplot2)
library(goeveg)
library(ggord)
library(tidyverse)
library(ggrepel)


bugs <- read_csv("bugid.csv") %>%
  filter(Season == "Fall")
bugs <- bugs[, colSums(bugs != 0, na.rm = TRUE) > 0]
bugs$group <- factor(bugs$Stream, levels = c("EAS (R)", "CRO (R)", 
                                             "FRY (L)", "SPC (L)", 
                                             "ROL (H)", "LLW (H)"))

data_1 <- bugs[,4:97] # Spring 93, Fall 97

nmds <- metaMDS(data_1, distance = "bray", k=2, autotransform=TRUE)

limited <- ordiselect(data_1, nmds, ablim = 0.15) # ablim = proportion of species with highest abundances
limited

scoresspecies <- as.data.frame(scores(nmds)$species) 
scoresspecies <- subset(scoresspecies, rownames(scoresspecies) %in% limited)
scoresspecies$taxa <- row.names(scoresspecies)

datascores <- as.data.frame(scores(nmds)$sites)
data_2 <- bugs[,98:99]
scores <- cbind(datascores, Stream = data_2$group)
scores <- cbind(scores, Impact = data_2$Impact)

stress <- nmds$stress
stress

#ggord(nmds, data_2$group, txt = limited) +
  #theme_classic()

plot <- ggplot(scores, aes(x= NMDS1, y = NMDS2)) +
         geom_point(aes(color = Impact, shape = Stream, size = 4)) +
  stat_ellipse(aes(data = Impact, color = Impact)) +
  geom_segment(data=scoresspecies, 
             aes(x=0, y=0, xend=NMDS1*1, yend=NMDS2*1), 
             colour="black", size=0.7, arrow=arrow()) +
  geom_text_repel(data=scoresspecies, 
                  aes(x=NMDS1*1.5, y=NMDS2*1.5, label=taxa),
                  colour="black", size = 7) +
  theme_classic() + 
  ggtitle("Fall 2021") +
  theme(plot.title = element_text(size = 20), axis.title.x =element_text(size = 20), axis.text.x = element_text(size = 18),
        axis.title.y =element_text(size = 20), axis.text.y  = element_text(size = 18),
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  annotate("text", x = 0.5, y = -0.7, label = "Stress = 0.22", size = 5)
plot

png("SpeciesNMDSALLS.png", width = 650, height = 550)
plot(plot)
dev.off()


###### NMDS by Stream ########

bugs <- read_csv("bugid.csv") %>%
  filter(Season == "Fall") %>%
  filter(Stream == "LLW (H)")
bugs <- bugs[, colSums(bugs != 0, na.rm = TRUE) > 0]

data_1 <- bugs[,4:50] # Fall: EAS 58, CRO 68, FRY 49, SPC 58, ROL and LLW 50

nmds <- metaMDS(data_1, distance = "bray", k=2, autotransform=TRUE)

limited <- ordiselect(data_1, nmds, ablim = 0.15) # ablim = proportion of species with highest abundances
limited

scoresspecies <- as.data.frame(scores(nmds)$species) 
scoresspecies <- subset(scoresspecies, rownames(scoresspecies) %in% limited)
scoresspecies$taxa <- row.names(scoresspecies)

datascores <- as.data.frame(scores(nmds)$sites)
data_2 <- bugs[,53] 
scores <- cbind(datascores, Site = data_2$Site.ord)

stress <- nmds$stress
stress

plot <- ggplot(scores, aes(x= NMDS1, y = NMDS2, label = Site)) +
  geom_point(color = "chocolate2", size = 4) +
  geom_text_repel(nudge_y = 0.06, size =8, color = "chocolate2") + #chartreuse3, blue, chocolate2 
  geom_segment(data=scoresspecies, 
               aes(x=0, y=0, xend=NMDS1*2.5, yend=NMDS2*2.5), 
               colour="black", arrow=arrow(),inherit.aes = FALSE) +
  geom_text_repel(data=scoresspecies, 
                  aes(x=NMDS1*3 ,y=NMDS2*3, label=taxa),
                  colour="black", size = 7) +
  theme_classic() + 
  ggtitle("LLW (H) Fall 2021") +
  theme(plot.title = element_text(size = 20), axis.title.x =element_text(size = 20), axis.text.x = element_text(size = 18),
        axis.title.y =element_text(size = 20), axis.text.y  = element_text(size = 18),
        legend.text=element_text(size=20), legend.title=element_text(size=18)) +
  annotate("text", x = -1, y = -0.3, label = "Stress = 0", size = 5)
plot

png("SpeciesNMDSALLS.png", width = 570, height = 500)
plot(plot)
dev.off()
