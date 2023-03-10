RDA Plotting
library(ggord)
library(vegan)
library(ggplot2)
library(tidyverse)

env <- read.csv("env.rda.csv") %>%
  filter( Season == "Spring")

d.S <- d.S[-c(8),]
com.hel <- com.hel[-c(8),]

#write.csv(d.S, file = "d.S.csv")

env.z <- env %>%
  select(-c(Stream, Site, Season)) %>%
  decostand( method = "standardize", na.rm =TRUE)

mod.best <- rda(formula = com.hel ~ pfines + avg.slope + D.pca + avgembedd + avgwetwidth +
                  sc.uScm + hardness.mgl + npoc.mgl + elements.pca + temp.C + do.mgl +       
                nutrients.pca + Rip.pca + PCNM1 + PCNM2 +PCNM4 + PCNM5, data = env.z)

R2adj <- RsquareAdj(mod.best)$adj.r.squared
R2adj

env$group <- factor(env$Stream, levels = c("EAS (R)", "CRO (R)", 
                                           "FRY (L)", "SPC (L)", 
                                           "ROL (H)", "LLW (H)"))

plot <- ggord(mod.best, env$Stream, exp = 0.07, repel = TRUE, 
              veclsz = 0.7, ptslab = TRUE, addsize = 4,cols = c("deepskyblue3",
                                                               "deepskyblue",
                                                               "seagreen3",
                                                               "red4",
                                                               "red2",
                                                               "seagreen")) + #facet = TRUE, xlims = -1, 1, ylims = -1,1
  labs(title = "All Streams Spring") +
  theme_classic() +
  annotate("text", x = -0.5, y = 0.8, label = "R2adj = 0.52") 
plot

png("RDAALL.png", width = 600, height = 600) # one faceted: width = 1000, height = 800
# 6 faceted width = 1800, height = 2000
plot(plot)
dev.off()
