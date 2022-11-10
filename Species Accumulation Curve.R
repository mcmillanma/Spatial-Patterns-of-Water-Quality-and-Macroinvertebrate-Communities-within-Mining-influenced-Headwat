#Species Accumulation Curve
#https://rpubs.com/Roeland-KINDT/694021
library("BiodiversityR") # also loads vegan
library(ggplot2)
library(ggsci)
library(readxl)
library(ggplot2)

Accum.1 <- accumcomp(warcom, y=warenv, factor='population', 
                     method='exact', conditioned=FALSE, plotit=FALSE)
accum.long1 <- accumcomp.long(Accum.1, ci=NA, label.freq=5)

com <- read_csv("bugid.csv") %>%
  filter(Season == "Spring") %>%
  select(-c(Site, Season))
chem <- read.csv("chem.f21-s22.notrib.csv") %>%
  filter(season == "Spring") #%>%
  filter(Site != "SPC6" ) %>%
  select(Stream, do.mgl:so4.hco3)
summary(chem)

hab <- read.csv("habitatmaster.csv")

Accum.1 <- accumcomp(com, y=com, factor='Stream', 
                     method='exact', conditioned=FALSE, plotit=FALSE)
Accum.1

accum.long1 <- accumcomp.long(Accum.1, ci=NA, label.freq=5)
head(accum.long1)

Spring <- ggplot(data=accum.long1, aes(x = Sites, y = Richness, color= Grouping)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  ggtitle("Species Accumulation Curve Spring 2022") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  labs(x = "Sample Number", y = "Species", colour = "Stream", shape = "Stream")
Fall

gg <- ggarrange(Fall, Spring, ncol= 1, common.legend = TRUE, legend="bottom")
gg

png("speciesaccum.png", width = 900, height = 1200)
plot(gg)
dev.off()
