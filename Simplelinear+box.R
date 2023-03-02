
library(ggpubr)

metrics <- read.csv("metrics.f21-s22.csv")
#distance.d <- read.csv("distance.csv")
chem <- read.csv("chem.f21-S22.notrib.reduced.csv") %>%
  dplyr::filter(season == "Spring")


#chem <- left_join(distance.d, chem, by = c("site" = "Site")) 
hab <- read.csv("habitatmaster.csv")
reg <- left_join(metrics, chem, by = c("Site", "Season" = "season", "dist.d", "Stream")) # %>%
  filter(Season == "Fall")


library("tidyverse")
library("multcompView")

# Boxplots between stream comparisons for reference
#preform anova and find pairwise values

shapiro.test(reg$pE)
library(car)
qqPlot(reg$sc.uScm)
leveneTest(sc.uScm ~ Stream, data = reg)

# https://www.r-bloggers.com/2021/08/how-to-perform-tukey-hsd-test-in-r/
set.seed(1045)
model <- aov(sc.uScm~Stream, data=reg)
summary(model)
TukeyHSD(model, conf.level=.95)
plot(TukeyHSD(model, conf.level=.95), las = 2)

#https://www.mathiasecology.com/code/add-tukeys-significant-letters-to-ggplots
letters.df <- data.frame(multcompLetters(TukeyHSD(
  aov(sc.uScm~Stream, data=reg))$Stream[,4])$Letters)

colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$Stream <- rownames(letters.df) #Create column based on rownames

placement <- reg %>% #We want to create a dataframe to assign the letter position.
  group_by(Stream) %>%
  summarise(quantile(sc.uScm)[4])

colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) #Merge dataframes
library(ggplot2)
longterm <- chem %>%
  filter( Site %in% c("EAS1", "CRO2", "FRY1", "LLW3", "ROL2", "SPC1"))

box.s <- chem %>% #Dataframe from which data will be drawn
  ggplot(aes(x = reorder(Stream, sc.uScm, median), y = sc.uScm, color = Impact)) +#Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0, size = 1) + #I like to set the color of boxplots to black with the alpha at 0 (fully transparent). I also like geom_jitter() but do not use it here for simplicity.
  geom_point(size = 3) +
    geom_point(data = longterm, aes(Stream, sc.uScm), color = "black") +
  theme_classic() +
  #scale_colour_viridis_d()+ #viridis color blind friendly
  xlab("Stream (Spring)") +
  ylab("SC (uS/cm)") +
  theme(axis.text = element_text(size = 18, color = "black"))#+
  #geom_text(data = letters.df, aes(x = Stream, y = Placement.Value, label = Letter), size = 4, color = "black", hjust = -1.25, vjust = -0.8, fontface = "bold")
box.s

box.f <- reg %>% #Dataframe from which data will be drawn
  ggplot(aes(x = reorder(Stream, sc.uScm, median), y = sc.uScm, color = Stream)) + #Instead of hard-coding a factor reorder, you can call it within the plotting function
  geom_boxplot(alpha = 0) + #I like to set the color of boxplots to black with the alpha at 0 (fully transparent). I also like geom_jitter() but do not use it here for simplicity.
  theme_classic() + #Clean, minimal theme courtesy of the "egg" package
  xlab("Stream (Fall)") +
  ylab("SC (uS/cm)") +
  geom_text(data = letters.df, aes(x = Stream, y = Placement.Value, label = Letter), size = 4, color = "black", hjust = -1.25, vjust = -0.8, fontface = "bold")
box.f

library(ggpubr)
box <- ggarrange(box.f, box.s, ncol = 2, common.legend = TRUE, legend = "bottom" )

png("box.tukeyhsd+aov.scxstreamxseason.png", width = 900, height = 450)
plot(box)
dev.off()

# option 2
box <- ggboxplot(reg, x = "Stream", y = "sc.uScm", color = "Season",
          palette = c("#00AFBB", "#E7B800"))
box

# Are streams community composition different 


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
