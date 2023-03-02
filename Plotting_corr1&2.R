### Ploting correlations ###

# metrics of interest: percent E, E less B richness, EPT richness, 
# scraper richness, clinger richness, and percent 5 dominant taxa
# pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom

chem <- read_csv("chem.pca.csv") 
metrics <- read_csv("metrics.f21-s22.csv") 
all <- left_join(chem, metrics, by = c( "Stream", "Site", "dist.d"))
all$group <- factor(all$Stream, levels = c("EAS (R)", "CRO (R)", 
                                           "FRY (MI)", "SPC (MI)", 
                                           "ROL (MI)", "LLW (MI)"))
boxplot(all$sc.uScm ~ all$group)


longterm <- all %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))

# pE, rich.E.less.B, rich.EPT, rich.SC, rich.Cling, p5dom

p5dom <- all %>% ggplot(aes(x=dist.d , y=p5dom, group = Season, color = Season))+
  geom_point(aes(color = Season, shape = Impact))  +
  stat_smooth(method = "lm", linetype = 2, colour = "black", alpha = 0.3) +
  stat_cor(method = "spearman", label.y.npc = 1, size = 5) + # label.y.npc = 0.3
  stat_regline_equation( label.y.npc = 0, label.x.npc =0.5, size = 5) + # label.x.npc = 0.5, label.y.npc = 0.15
  geom_point(data = longterm, aes(x=dist.d , y=p5dom, shape = Season), 
             color = "black") +
  facet_wrap("group") +
  theme_classic()+
  xlab("Distance Downstream (m)") +
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

png("slope.png",width = 900, height = 600) # one faceted: width = 900, height = 600
# 6 faceted
plot(p5dom)
dev.off()
