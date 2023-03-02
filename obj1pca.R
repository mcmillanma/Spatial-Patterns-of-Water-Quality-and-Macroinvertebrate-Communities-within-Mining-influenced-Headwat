# RDA for metrics and Species
library(tidyverse)
library(vegan) #decostand
library (car) #vif
library(ggord)

com <- read.csv("bugid.csv") %>%
  filter(Season == "Spring") %>%
  #filter(Stream == "LLW (MI)") %>%
  dplyr::select(-c(Stream:Site))
com <- com[, colSums(com != 0, na.rm = TRUE) > 0]
spe.hel <- decostand(com, "hellinger")

met <- read.csv("metrics.reduced.csv")%>%
  filter(Season == "Spring") %>%
  #filter(Stream == "LLW (MI)") %>%
  dplyr::select(-c(Stream:Season, rich.E:rich.T)) 
met.z <- decostand(met, method = "standardize",  na.rm =TRUE)
met.z <- met.z %>%
  select_if(~ ! any(is.na(.)))

bugs.pca <- prcomp(spe.hel)
summary(bugs.pca)
str(bugs.pca)
ggord(bugs.pca, xlims= c(-1,1) , ylims = c(-1,1))
bugs.pca
bugscore <- as.data.frame(bugs.pca$rotation)

loadings=NULL
loadings$All.i <- bugscore$PC1
loadings <- as.data.frame(loadings)

sites = as.data.frame(scores(bugs.pca))
env <- met.z
ef <- envfit(bugs.pca, env)
ef
spec.data.envfit <- data.frame(r=ef$vectors$r, p=ef$vectors$pvals)
species.long3 <- spec.data.envfit[spec.data.envfit$r >= 0.6, ]
species.long3$variables <- rownames(species.long3)

ef <- left_join(species.long3, vec.df, by = "variables" )

ggord(bugs.pca)
llw = ggord(bugs.pca, exp = 0.05, repel = TRUE, ext = 1, max.overlaps = 100,
           size = 5, veclsz = 0.5, labcol = "red" ,obslab = TRUE, ptslab = TRUE) +
  geom_segment(data = ef,
                    aes(x = 0, xend = PC1*4, y = 0, yend = PC2*4),
                    colour="blue",
                  inherit.aes = FALSE, arrow = arrow()) + 
  geom_text(data = ef,
            aes(x = PC1*5, y=PC2*5.2, label = variables),
            size=4, color = "blue") +
  ggtitle("FRY (MI) Spring ") +
  theme_classic()

llw
fry

library(ggpubr)
  ggarrange(llw,fry)

png("PCACROS.png", width = 900, height = 700)
plot(xx)
dev.off()

