# PART 1) Simple regression of main metrics known to differ between streams by stream
# Part 2) Determine which habitat and water chemistry metrics are correlated to which insect metrics using NMDS. One NMDS graph with vectors of water quality and habitat on-top of raw insect community data.
# PART 3) Step-wise multiple regression of habitat and water quality to community similarity

#PART 1: metrics of interest: bugs) E.rich, E.rich.less.B, %E less B, EPT.rich, SC.rich, %E
# shannon, total richness, %5 Dominant Taxa, predator richness, scraper richness; 
# Water) SC, Se, SO4:HCO3; Ca:Mg ; Habitat) watershed size, stream gradient, LCF

# determine if data is normal and then create correlation matrix of sc to all 
#64 insect metrics

library(dplyr)

# prepare data to be used (sc and all metrics)
metrics <- read.csv("metrics.f21.csv")
chem <- read.csv("chem.f21-S22.csv") %>%
  filter( season == "Fall")
c.df <- left_join(metrics, chem, by = "Site")
c.df <- select(c.df, c(3, 6:57, 60:69, 76))

# Split by stream and determin if data is normal use rcorr instead of cor to get p values
#CRO
library(tibble)
library(Hmisc)
cro <- filter(c.df, Stream.x == "CRO (R)") %>%
  select(-c(1))
cro <- cor(cro, method = c("pearson"), use = "complete.obs") #cro <- rcorr(as.matrix(cro), type = "spearman")
cro <- as.data.frame(cro) %>%
  select("sc.uScm")
rownames <- rownames(cro)
cro <- add_column(cro, rownames)

set.seed(123)
norm.cro <- apply(cro,2,shapiro.test)
norm.cro <- do.call(rbind.data.frame, norm.cro)
norm.cro <- norm.cro %>%
  mutate("test" = case_when(p.value <= 0.05 ~ "Pearson", 
                            p.value >0.05 ~ "Spearman"))

filter(norm.cro, test == "Pearson")

#EAS
eas <- filter(c.df, Stream.x == "EAS (R)") %>%
  select(-c(1))
eas <- cor(eas, method = c("pearson"), use = "complete.obs") #eas <- rcorr(as.matrix(eas), type = "spearman")
eas <- as.data.frame(eas) %>%
  select("sc.uScm")
rownames <- rownames(eas)
eas <- add_column(eas, rownames)

norm.eas <- apply(eas,2,shapiro.test)
norm.eas <- do.call(rbind.data.frame, norm.eas)
norm.eas <- norm.eas %>%
  mutate("test" = case_when(p.value <= 0.05 ~ "Pearson", 
                            p.value >0.05 ~ "Spearman"))
filter(norm.eas, test == "Pearson")

#FRY
fry <- filter(c.df, Stream.x == "FRY (MNG)") %>%
  select(-c(1, 62, 63))
fry <- cor(fry, method = c("pearson"), use = "complete.obs") # fry <- rcorr(as.matrix(fry), type = "spearman") 
fry <- as.data.frame(fry) %>%
  select("sc.uScm")
rownames <- rownames(fry)
fry <- add_column(fry, rownames)

norm.fry <- apply(fry,2,shapiro.test)
norm.fry <- do.call(rbind.data.frame, norm.fry)
norm.fry <- norm.fry %>%
  mutate("test" = case_when(p.value <= 0.05 ~ "Pearson", 
                            p.value >0.05 ~ "Spearman"))

#LLW
llw <- filter(c.df, Stream.x == "LLW (MG)") %>%
  select(-c(1))
llw <- cor(llw, method = c("pearson"), use = "complete.obs") #llw <- rcorr(as.matrix(llw), type = "spearman")
llw <- as.data.frame(llw) %>%
  select("sc.uScm")
rownames <- rownames(llw)
llw <- add_column(llw, rownames)

norm.llw <- apply(llw,2,shapiro.test)
norm.llw <- do.call(rbind.data.frame, norm.llw)
norm.llw <- norm.llw %>%
  mutate("test" = case_when(p.value <= 0.05 ~ "Pearson", 
                            p.value >0.05 ~ "Spearman"))
filter(norm.llw, test == "Pearson")

#ROL
rol <- filter(c.df, Stream.x == "ROL (MG)") %>%
  select(-c(1, 39))
rol <- cor(rol, method = c("spearman"), use = "complete.obs") # rol <- rcorr(as.matrix(rol), type = "spearman") 
rol <- as.data.frame(rol) %>%
  select("sc.uScm")
rownames <- rownames(rol)
rol <- add_column(rol, rownames)

norm.rol <- apply(rol,2,shapiro.test)
norm.rol <- do.call(rbind.data.frame, norm.rol)
norm.rol <- norm.rol %>%
  mutate("test" = case_when(p.value <= 0.05 ~ "Pearson", 
                            p.value >0.05 ~ "Spearman"))


#SPC
spc <- filter(c.df, Stream.x == "SPC (MG)") %>%
  select(-c(1))
spc <- cor(spc, method = c("spearman"), use = "complete.obs") #spc <- rcorr(as.matrix(spc), type = "spearman") 
spc <- as.data.frame(spc) %>%
  select("sc.uScm")
rownames <- rownames(spc)
spc <- add_column(spc, rownames)

cor.scxstream <- left_join(eas, cro, by = "rownames") %>%
  left_join(llw, by = "rownames") %>%
  left_join(spc, by = "rownames") %>%
  left_join(rol, by = "rownames") %>%
  left_join(fry, by = "rownames") 
write.csv(cor.scxstream, file="cor.scxstream.f21.csv", sep = ",")

norm.spc <- apply(spc,2,shapiro.test)
norm.spc <- do.call(rbind.data.frame, norm.spc)
norm.spc <- norm.spc %>%
  mutate("test" = case_when(p.value <= 0.05 ~ "Pearson", 
                            p.value >0.05 ~ "Spearman"))

filter(norm.spc, test == "Pearson")

#Determine if data is normal globally
set.seed(123)
norm.results <- apply(c.df,2,shapiro.test)

norm.results <- do.call(rbind.data.frame, norm.results)

# if data is normal use Pearson, or if non-normal, use Spearman rank correlations
norm.results <- norm.results %>%
mutate("test" = case_when(p.value <= 0.05 ~ "Pearson", 
                                        p.value >0.05 ~ "Spearman"))



#install.packages("lmerTest")
library(ggplot2)
library(ggpubr)
library(dplyr)

metrics <- read.csv("metrics.f21.csv")
chem <- read.csv("chem.f21-S22.csv")
chem <- chem %>%
  mutate("ratiocamg" = ca.mgl / mg.mgl) %>%
  mutate ("ratioso4hco3" = so4.mgl /hco3.mgl)
chem <- dplyr::filter(chem, season == "Fall") 
hab <- read.csv("habitatmaster.csv")
reg2 <- left_join(metrics, chem, by = "Site") 
reg2 <- left_join(reg2, hab, by = "Site")

# Boxplots. Top 3 that change between my streams pEPT.less.HBL, pE.less.B, and pT
# most expected to change with SC rich.SC, rich.EPT, rich.E.less.B, p5dom, Hshannon
my_comparisons <- list( c("Reference 1", "Reference 2") )
compare_means(rich.SC ~ Stream.Type, data = reg2)

box.hshanxstream <- ggplot(reg2, aes(x=Stream.Type, y=Hshannon, group = Stream.Type, color = Stream.Type))+
  stat_boxplot() +
  theme_classic() +
  ylab("Shannon Diversity") +
  xlab("Stream") +
  ggtitle("Shannon Diversity by Stream") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  stat_compare_means(method = "kruskal", label.x.npc = "left", label.y = 3) +
  stat_compare_means(label = "p.signif", method = "anova.test",
                     ref.group = ".all.") 
  #stat_compare_means(comparisons = my_comparisons)+
box.hshanxstream

png("box.hshanxstream.png", width = 900, height = 450)
plot(box.hshanxstream)
dev.off()

# Top insect vs sc metrics of interest: rich.SC, rich.EPT, rich.E.less.B, p5dom, Hshannon
# Top insect vs water quality metrics: sc.uScm, se.ugl, ca/mg, so4/hco3; Hshannon and rich.EPT
# Top insect vs habitat metrics: LCF, avg.slope, avgembedd, pfines; rich.PR and rich.SC

scxdistancexstream <- ggplot(bug.dist.cor,aes(x=dis.simm, y=sc.uScm, group = Stream.Type))+
  geom_point(aes(color=Stream.Type))+
  # geom_jitter() +
  stat_smooth(method = "lm", linetype = 2) +
  stat_cor(method = "spearman",
           label.y = 150) +
  stat_regline_equation(label.y = 320) + #label.y = 
  #ylim(0,8) +
  #xlim(0, 2500) +
  theme_classic(base_size = 15, )+
  xlab("Distance Downstream (m)") +
  ylab("Specific Conductance (uS/cm)") +
  ggtitle("Specific Conductance vs Distance Downstream")+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  facet_wrap("Stream.Type") +
  theme(legend.position = "none") +
  theme_classic()
scxdistancexstream


#set size and save plot as png
png("scxdistancexstream.png", width = 900, height = 450)
plot(scxdistancexstream)
dev.off()


#PART 2 
#install.packages("glmulti")
#library(glmulti)
library(dplyr)
library(car)

# prepare data to be used
sim <- read.csv("Similarity.f21.csv")
sim <- sim %>% filter(START %in% c("CRO-1", "EAS-1", "LLW-1", "ROL-1", "SPC-1", "FRY-2")) %>%
  filter(END != "FRY-1") %>%
  filter(END != "SPC-6") %>%
  filter(END != "FRY-8")
chem <- read.csv("chem.f21-S22.csv")
chem <- filter( chem, season == "Fall") 
multireg <- left_join(sim, chem, by = c("END" = "site.id")) 
multireg <- select(multireg, -c(1:5, 8:12, 28, 46:54))

which(is.na(multireg))

library(lme4)

#define intercept-only model
intercept_only <- lm( SIM ~ 1, data=multireg)
summary(intercept_only)

#define model with all predictors
colnames(multireg)
all <- lm(SIM ~ temp.C + hardness.mgl + no2no3.n.mgl + so4.mgl + k.mgl + ti.ugl +
            mn.ugl + zn.ugl + dis.simm + sc.uScm + tn.mgl + ortho.po4.p.mgl + hco3.mgl 
          + ca.mgl + v.ugl + co.ugl + as.ugl +do.mgl + ac.uScm + tp.mgl + npoc.mgl + 
            na.mgl + li.ugl + ph + alkalinity.mgl + nh3.n.mgl +
            cl.mgl + mg.mgl + al.ugl + se.ugl + fe.ugl + cu.ugl, data=multireg)
summary(all)
plot(all)

#perform forward stepwise regression
forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)
plot(forward)
#view results of forward stepwise regression
forward$anova
