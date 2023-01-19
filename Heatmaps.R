# Habitat
hab <- read.csv("habitatmaster.csv") %>%
  filter(Stream == "SPC (MI)") %>%
  select(-c(Stream:smallcobble)) #%>%
  slice(-(8), -(31))

habheat <- heatmap(abs(cor(hab, method = "spearman")), 
                   # Compute pearson correlation (note they are absolute values)
                   col = rev(heat.colors(6)), 
                   Colv = NA, Rowv = NA)
legend("topright", 
       title = "SPC",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
habheat

#Chem

chem <- read.csv("chem.f21-s22.notrib.reduced.csv")%>%
  filter(season == "Spring") %>%
  #filter(Stream == "ROL (MI)") %>%
  select(-c("Stream", Site:season)) %>%
  slice( -(32))
  

chemheat <- heatmap(abs(cor(chem)), 
                   # Compute pearson correlation (note they are absolute values)
                   col = rev(heat.colors(6)), 
                   Colv = NA, Rowv = NA) 
legend("topright", 
       title = "Chem Spring All",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

