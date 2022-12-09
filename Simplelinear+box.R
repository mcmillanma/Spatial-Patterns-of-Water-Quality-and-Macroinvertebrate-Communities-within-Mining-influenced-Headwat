
library(ggpubr)

metrics <- read.csv("metrics.f21-s22.csv")
#distance.d <- read.csv("distance.csv")
chem <- read.csv("chem.f21-S22.notrib.csv") 
#chem <- left_join(distance.d, chem, by = c("site" = "Site")) 
hab <- read.csv("habitatmaster.csv")
reg <- left_join(metrics, chem, by = c("Site", "Season" = "season", "dist.d", "Stream")) 



library(rstatix)
library("tidyverse")
library("egg") #The egg package contains one of my favorite themes, theme_article.
library("multcompView")
library("sigminer")

# Table of metrics to SC and dist.d correlation
letters.df <- data.frame(multcompLetters(TukeyHSD(aov(rich.SC ~ Stream.y, 
                                                      data = reg))$Stream.y[,4])$Letters)

colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Stream.y <- rownames(letters.df) #Create column based on rownames

placement <- reg %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream.y) %>%
  summarise(quantile(rich.SC)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) #Merge dataframes


# Boxplots between stream comparisons for reference

box <- reg %>%
  ggplot(aes(x= Stream.y, y= rich.SC, group = Stream.y, color = Stream.y)) +
  geom_boxplot( alpha = 0) +
  facet_wrap("Season") + 
  geom_text(data = letters.df, aes(x = Stream.y, y = Placement.Value, label = Letter),
            size = 4, color = "black", hjust = -1.25, vjust = -0.8, fontface = "bold")+
  theme_classic() +
  ylab("Scraper Richness") +
  xlab("Stream") +
  ggtitle("Scraper Richness by Stream") +
  guides(color = guide_legend(title = "Stream")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  stat_compare_means(method = "kruskal", label.x.npc = "left", label.y.npc = "bottom") 
box

box <- reg %>%
  ggplot(aes(x= Stream.y, y= rich.SC, group = Stream.y, color = Stream.y)) +
  geom_boxplot( alpha = 0) +
  facet_wrap("Season") + 
  #geom_text(data = letters.df, aes(x = Stream.y, y = Placement.Value, label = Letter),
            #size = 4, color = "black", hjust = -1.25, vjust = -0.8, fontface = "bold")+
  theme_classic() +
  ylab("Scraper Richness") +
  xlab("Stream") +
  ggtitle("Scraper Richness by Stream") +
  guides(color = guide_legend(title = "Stream")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  stat_compare_means(method = "kruskal", label.x.npc = "left", label.y.npc = "bottom") +
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.")

box

pairs <- compare_means(rich.SC ~ Stream.y, data = reg, 
              group.by = "Season", method = "wilcox.test")

png("box.tukeyhsd+aov.richElessBxstreamxseason.png", width = 900, height = 450)
plot(box)
dev.off()

# Linear regressions SC and distance
# Top insect metrics vs SC Timpano: rich.E.less.B, pElessB, rich.E, richEPT, rich.SC, pE, rich, Hshannon, 5dom
# Top insect metrics vs SC Cianciolo: rich, rich.EPT, rich.E, rich.SC, Hshannon, pE
# Drover Multiple stressors: SC - rich.E, pE; Se - Predator density; LCF- rich.Cling; LRBS - rich.Cling

#My Metrics of interest: rich.SC, rich.Cling, rich.EPT, 5dom, Hshannon, rich.PR
# rich.E, rich.EPT, rich.SC, rich, Shannon, rich.Cling (habitat)

library (ggplot2)

plot3 <- ggplot(reg ,aes(x=dist.d, y=rich.Cling, group = Season, color = Season))+
  geom_point()+
  # geom_jitter() +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman", label.y.npc = 1) +
  stat_regline_equation(label.x.npc = 0.5, label.y.npc = 0) +
  #ylim(0,8) +
  #xlim(0, 2500) +
  theme_classic() +
  xlab("Distance Downstream") +
  ylab("Clinger Richness") +
  ggtitle("Clinger Richness vs Distance Downstream") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  labs(color="Season") +
  facet_wrap("Stream") 
plot3

plot <- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")


plot7 <- ggplot(reg ,aes(x=sc.uScm, y=rich.E.less.B, group = Season, color = Season))+
  geom_point()+
  # geom_jitter() +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman") +
  stat_regline_equation(label.x.npc = 0.5, label.y.npc = 0.95) +
  #ylim(0,8) +
  #xlim(0, 2500) +
  theme_classic() +
  xlab("Specific Conductance (uS/cm)") +
  ylab("EPT less B") +
  ggtitle("EPT less B vs Specific Conductance (uS/cm)") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) + 
  facet_wrap("Stream")  

plot4

plotg <- ggarrange(plot1g, plot2g, plot3g, plot4g, plot5g, plot6g, ncol = 3, nrow = 2, 
                  common.legend = TRUE, legend = "bottom")

plot <- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 3, nrow = 2, 
                  common.legend = TRUE, legend = "bottom")
ggarrange(plot2, plot7)

#set size and save plot as png
png("metricsxsc.season.png", width = 1700, height = 1000)
plot(plot)
dev.off()

Sumj <-Similar %>%
  group_by(Stream, Season) %>%
  slice_min(SIM.j) 
Sumb <-Similar %>%
  group_by(Stream, Season) %>%
  slice_min(SIM.b) 


write.csv(Sumb, file="sumb.csv", sep = ",")
