### Ploting correlations ###
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)

# metrics of interest: percent E, E less B richness, EPT richness, 
# scraper richness, clinger richness, and percent 5 dominant taxa
# pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom

chem <- read_csv("chem.pca.csv") 
metrics <- read_csv("metrics.f21-s22.csv") 
all <- left_join(chem, metrics, by = c( "Stream", "Site", "dist.d", "Season"))
all$group <- factor(all$Stream, levels = c("EAS (R)", "CRO (R)", 
                                           "FRY (MI)", "SPC (MI)", 
                                           "ROL (MI)", "LLW (MI)"))
boxplot(all$sc.uScm ~ all$group)


longterm <- all %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))

# pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom

p5dom <- all %>% ggplot(aes(x=sc.uScm , y=p5dom, group = Season, color = Season))+
  geom_point(aes(color = Season, shape = Impact), size = 4)  +
  stat_smooth(method = "lm", linetype = 2, colour = "black", alpha = 0.3) +
  stat_cor(method = "spearman", label.y.npc = 1, size = 6) + # label.y.npc = 0.3
  stat_regline_equation( label.y.npc = 0, label.x.npc =0.5, size = 6) + # label.x.npc = 0.5, label.y.npc = 0.15
  geom_point(data = longterm, aes(x=sc.uScm , y=p5dom, shape = Season), 
             color = "black" , size = 4) +
  facet_wrap("group") +
  theme_classic()+
  xlab("Specific Conductance (uS/cm)") +
  ylab("% 5 Dominant Taxa") +
  scale_y_continuous() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))

pE
rich.E.less.B
rich.EPT
rich.SC
rich.Cling
p5dom

plot <- ggarrange(ggarrange(pE, rich.E.less.B, common.legend = TRUE, legend = "right"), 
ggarrange( rich.EPT, rich.SC,common.legend = TRUE, legend = "right"), 
ggarrange(rich.Cling, p5dom ,common.legend = TRUE, legend = "right"), common.legend = TRUE, 
legend = "right", ncol = 1)

plot

png("slope.png", width = 1800, height = 2000) # one faceted: width = 1000, height = 800
# 6 faceted width = 1800, height = 2000
plot(plot)
dev.off()

##### Global ####

# pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom

rich.Cling <- all %>% ggplot(aes(x=sc.uScm , y=rich.Cling, group = Season, color = Season))+
  geom_point(aes(color = Season, shape = Stream), size = 4)  +
  stat_smooth(method = "lm", linetype = 2, colour = "black", alpha = 0.3) +
  stat_cor(method = "spearman", label.y.npc = 1, size = 6) + # label.y.npc = 0.3
  stat_regline_equation( label.y.npc = 0.85, label.x.npc =0, size = 6) + # label.x.npc = 0.5, label.y.npc = 0.15
  geom_point(data = longterm, aes(x=sc.uScm , y=rich.Cling, shape = Season), 
             color = "black" , size = 4) +
  #facet_wrap("group") +
  theme_classic()+
  xlab("Specific Conductance (uS/cm)") +
  ylab("Clinger Richness") +
  #stat_ellipse(aes(color = Impact????)) +
  scale_y_continuous() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x =element_text(size = 20, color = "black"),
        axis.title.y =element_text(size = 20, color = "black"),
        legend.text=element_text(size=18, color = "black"), 
        legend.title=element_text(size=20, color = "black"),
        strip.text.x = element_text(size = 18))

pE
rich.E.less.B
rich.EPT
rich.SC
rich.Cling
p5dom

plot <- ggarrange(pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom,  
                  common.legend = TRUE, legend = "right")

plot

png("slope.png", width = 1800, height = 1200) # one faceted: width = 1000, height = 800
# 3 x2 faceted width = 1800, height = 2000
# 2 x 3 global
plot(plot)
dev.off()




