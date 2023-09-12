require(tidyverse)
library(raster)
library(readxl)
library(hrbrthemes)
library(viridis)
library(gridExtra)

setwd("D:/PROJECTS/Rafflesia_SDM")

#read in envi data
# from points
envipts <- read_excel("Species_pts_30s_envdata.xlsx")
summary(envipts)
envipts$Species <- as.factor(envipts$Species)
envipts <- filter(envipts, Species == "Rafflesia lobata" | Species == "Rafflesia lagascae" |Species == "Rafflesia speciosa" )

# from SDM results
envisdm <- read_excel("All_Rafflesia_spp_suitable_masked_envdata.xlsx")
summary(envisdm)
envisdm$Species <- as.factor(envisdm$Species)
envisdm$Scenario <- as.factor(envisdm$Scenario)

# native ranges only
envintv <- read.csv("Island_analysis/All_Rafflesia_spp_a1_suitable_masked_native_ranges_envdata.csv")
summary(envintv)
envintv$Species <- as.factor(envintv$Species)
envintv$Scenario <- as.factor(envintv$Scenario)

#filter data to ind species
#R. speciosa
rspe.pts <- filter(envipts, Species == "Rafflesia speciosa")
rspe.pts$Coverage <- rep("Species points", dim(rspe.pts)[1])

rspe.ph <- filter(envisdm, Species == "Rspe") #entire PH
rspe.ph$Species <- rep("Rafflesia speciosa", dim(rspe.ph)[1])
rspe.ph$Coverage <- rep("Philippine extent SDM", dim(rspe.ph)[1])

rspe.nt <- filter(envintv, Species == "Rspe" ) #native range only
rspe.nt$Species <- rep("Rafflesia speciosa", dim(rspe.nt)[1])
rspe.nt$Coverage <- rep("Native range SDM", dim(rspe.nt)[1])


#R lobata
rlob.pts <- filter(envipts, Species == "Rafflesia lobata")
rlob.pts$Coverage <- rep("Species points", dim(rlob.pts)[1])

rlob.ph <- filter(envisdm, Species == "Rlob") #entire PH
rlob.ph$Species <- rep("Rafflesia lobata", dim(rlob.ph)[1])
rlob.ph$Coverage <- rep("Philippine extent SDM", dim(rlob.ph)[1])

rlob.nt <- filter(envintv, Species == "Rlob") #native range only
rlob.nt$Species <- rep("Rafflesia lobata", dim(rlob.nt)[1])
rlob.nt$Coverage <- rep("Native range SDM", dim(rlob.nt)[1])


#R lagascae 
rlag.pts <- filter(envipts, Species == "Rafflesia lagascae")
rlag.pts$Coverage <- rep("Species points", dim(rlag.pts)[1])

rlag.ph <- filter(envisdm, Species == "Rlag") #entire PH
rlag.ph$Species <- rep("Rafflesia lagascae", dim(rlag.ph)[1])
rlag.ph$Coverage <- rep("Philippine extent SDM", dim(rlag.ph)[1])

rlag.nt <- filter(envintv, Species == "Rlag") #native range only
rlag.nt$Species <- rep("Rafflesia lagascae", dim(rlag.nt)[1])
rlag.nt$Coverage <- rep("Native range SDM", dim(rlag.nt)[1])



#Combine points and SDM results 
# April 20 2023: removing occurrence points data from the sets

#altitude
#df1 <- dplyr::select(rspe.pts, Species, Coverage, altitude)
df2 <- dplyr::select(rspe.ph, Species, Coverage, altitude)
df3 <- dplyr::select(rspe.nt, Species, Coverage, altitude)
rspe.alt <- bind_rows(df2, df3)

#df1 <- dplyr::select(rlob.pts, Species, Coverage, altitude)
df2 <- dplyr::select(rlob.ph, Species, Coverage, altitude)
df3 <- dplyr::select(rlob.nt, Species, Coverage, altitude)
rlob.alt <- bind_rows(df2, df3)

#df1 <- dplyr::select(rlag.pts, Species, Coverage, altitude)
df2 <- dplyr::select(rlag.ph, Species, Coverage, altitude)
df3 <- dplyr::select(rlag.nt, Species, Coverage, altitude)
rlag.alt <- bind_rows(df2, df3)

alt.all <- bind_rows(rspe.alt, rlob.alt, rlag.alt)
alt.all$Coverage <- factor(alt.all$Coverage, levels = c("Native range SDM", "Philippine extent SDM"), ordered = T)

# temperature
#df1 <- dplyr::select(rspe.pts, Species, Coverage, bio1)
df2 <- dplyr::select(rspe.ph, Species, Coverage, bio1)
df3 <- dplyr::select(rspe.nt, Species, Coverage, bio1)
rspe.temp <- bind_rows(df2, df3)

#df1 <- dplyr::select(rlob.pts, Species, Coverage, bio1)
df2 <- dplyr::select(rlob.ph, Species, Coverage,  bio1)
df3 <- dplyr::select(rlob.nt, Species, Coverage, bio1)
rlob.temp <- bind_rows(df2, df3)

#df1 <- dplyr::select(rlag.pts, Species, Coverage, bio1)
df2 <- dplyr::select(rlag.ph, Species, Coverage, bio1)
df3 <- dplyr::select(rlag.nt, Species, Coverage, bio1)
rlag.temp <- bind_rows(df2, df3)

temp.all <- bind_rows(rspe.temp, rlob.temp, rlag.temp)
temp.all$Coverage <- factor(temp.all$Coverage, levels = c("Native range SDM", "Philippine extent SDM"), ordered = T)


# precipitation
#df1 <- dplyr::select(rspe.pts, Species, Coverage, bio12)
df2 <- dplyr::select(rspe.ph, Species, Coverage, bio12)
df3 <- dplyr::select(rspe.nt, Species, Coverage, bio12)
rspe.pre <- bind_rows(df2, df3)

#df1 <- dplyr::select(rlob.pts, Species, Coverage, bio12)
df2 <- dplyr::select(rlob.ph, Species, Coverage, bio12)
df3 <- dplyr::select(rlob.nt, Species, Coverage, bio12)
rlob.pre <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlag.pts, Species, Coverage,  bio12)
df2 <- dplyr::select(rlag.ph, Species, Coverage, bio12)
df3 <- dplyr::select(rlag.nt, Species, Coverage, bio12)
rlag.pre <- bind_rows( df2, df3)

pre.all <- bind_rows(rspe.pre, rlob.pre, rlag.pre)
pre.all$Coverage <- factor(pre.all$Coverage, levels = c("Native range SDM", "Philippine extent SDM"), ordered = T)
summary(pre.all)

# temp seasonality
#df1 <- dplyr::select(rspe.pts, Species, Coverage, bio4)
df2 <- dplyr::select(rspe.ph, Species, Coverage, bio4)
df3 <- dplyr::select(rspe.nt, Species, Coverage, bio4)
rspe.ts <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlob.pts, Species, Coverage, bio4)
df2 <- dplyr::select(rlob.ph, Species, Coverage, bio4)
df3 <- dplyr::select(rlob.nt, Species, Coverage, bio4)
rlob.ts <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlag.pts, Species, Coverage, bio4)
df2 <- dplyr::select(rlag.ph, Species, Coverage, bio4)
df3 <- dplyr::select(rlag.nt, Species, Coverage, bio4)
rlag.ts <- bind_rows( df2, df3)

ts.all <- bind_rows(rspe.ts, rlob.ts, rlag.ts)
ts.all$Coverage <- factor(ts.all$Coverage, levels = c("Native range SDM", "Philippine extent SDM"), ordered = T)


# precip seasonality
#df1 <- dplyr::select(rspe.pts, Species, Coverage, bio15)
df2 <- dplyr::select(rspe.ph, Species, Coverage, bio15)
df3 <- dplyr::select(rspe.nt, Species, Coverage, bio15)
rspe.ps <- bind_rows(df2, df3)

#df1 <- dplyr::select(rlob.pts, Species, Coverage, bio15)
df2 <- dplyr::select(rlob.ph, Species, Coverage, bio15)
df3 <- dplyr::select(rlob.nt, Species, Coverage, bio15)
rlob.ps <- bind_rows(df2, df3)

#df1 <- dplyr::select(rlag.pts, Species, Coverage, bio15)
df2 <- dplyr::select(rlag.ph, Species, Coverage, bio15)
df3 <- dplyr::select(rlag.nt, Species, Coverage, bio15)
rlag.ps <- bind_rows( df2, df3)

ps.all <- bind_rows(rspe.ps, rlob.ps, rlag.ps)
ps.all$Coverage <- factor(ps.all$Coverage, levels = c("Native range SDM", "Philippine extent SDM"), ordered = T)

# soils, adding climate scenarios here for the revision (9/11/2023)

#soil ph
#df1 <- dplyr::select(rspe.pts, Species, Coverage, soilph)
df2 <- dplyr::select(rspe.ph, Species, Coverage, Scenario, soilph)
df3 <- dplyr::select(rspe.nt, Species, Coverage, Scenario, soilph)
rspe.sph <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlob.pts, Species, Coverage, soilph)
df2 <- dplyr::select(rlob.ph, Species, Coverage, Scenario, soilph)
df3 <- dplyr::select(rlob.nt, Species, Coverage, Scenario, soilph)
rlob.sph <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlag.pts, Species, Coverage, soilph)
df2 <- dplyr::select(rlag.ph, Species, Coverage, Scenario, soilph)
df3 <- dplyr::select(rlag.nt, Species, Coverage, Scenario, soilph)
rlag.sph <- bind_rows( df2, df3)

sph.all <- bind_rows(rspe.sph, rlob.sph, rlag.sph)
sph.all$Coverage <- factor(sph.all$Coverage, levels = c("Native range SDM", "Philippine extent SDM"), ordered = T)


#soil CEC
#df1 <- dplyr::select(rspe.pts, Species, Coverage, cec)
df2 <- dplyr::select(rspe.ph, Species, Coverage, Scenario, cec)
df3 <- dplyr::select(rspe.nt, Species, Coverage, Scenario, cec)
rspe.cec <- bind_rows(df2, df3)

#df1 <- dplyr::select(rlob.pts, Species, Coverage, cec)
df2 <- dplyr::select(rlob.ph, Species, Coverage, Scenario, cec)
df3 <- dplyr::select(rlob.nt, Species, Coverage, Scenario, cec)
rlob.cec <- bind_rows(df2, df3)

#df1 <- dplyr::select(rlag.pts, Species, Coverage, cec)
df2 <- dplyr::select(rlag.ph, Species, Coverage, Scenario, cec)
df3 <- dplyr::select(rlag.nt, Species, Coverage, Scenario, cec)
rlag.cec <- bind_rows(df2, df3)

cec.all <- bind_rows(rspe.cec, rlob.cec, rlag.cec)
cec.all$Coverage <- factor(cec.all$Coverage, levels = c("Native range SDM", "Philippine extent SDM"), ordered = T)

#soil bulk density
#df1 <- dplyr::select(rspe.pts, Species, Coverage, bulk_density)
df2 <- dplyr::select(rspe.ph, Species, Coverage, Scenario, bulk_density)
df3 <- dplyr::select(rspe.nt, Species, Coverage, Scenario, bulk_density)
rspe.bd <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlob.pts, Species, Coverage, bulk_density)
df2 <- dplyr::select(rlob.ph, Species, Coverage, Scenario, bulk_density)
df3 <- dplyr::select(rlob.nt, Species, Coverage, Scenario, bulk_density)
rlob.bd <- bind_rows(df2, df3)

#df1 <- dplyr::select(rlag.pts, Species, Coverage, bulk_density)
df2 <- dplyr::select(rlag.ph, Species, Coverage, Scenario, bulk_density)
df3 <- dplyr::select(rlag.nt, Species, Coverage, Scenario, bulk_density)
rlag.bd <- bind_rows(df2, df3)

bd.all <- bind_rows(rspe.bd, rlob.bd, rlag.bd)
bd.all$Coverage <- factor(bd.all$Coverage, levels = c("Native range SDM", "Philippine extent SDM"), ordered = T)

#soil sand
#df1 <- dplyr::select(rspe.pts, Species, Coverage, sand)
df2 <- dplyr::select(rspe.ph, Species, Coverage, Scenario, sand)
df3 <- dplyr::select(rspe.nt, Species, Coverage, Scenario, sand)
rspe.sand <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlob.pts, Species, Coverage, sand)
df2 <- dplyr::select(rlob.ph, Species, Coverage, Scenario, sand)
df3 <- dplyr::select(rlob.nt, Species, Coverage, Scenario, sand)
rlob.sand <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlag.pts, Species, Coverage, sand)
df2 <- dplyr::select(rlag.ph, Species, Coverage, Scenario, sand)
df3 <- dplyr::select(rlag.nt, Species, Coverage, Scenario, sand)
rlag.sand <- bind_rows( df2, df3)

sand.all <- bind_rows(rspe.sand, rlob.sand, rlag.sand)
sand.all$Coverage <- factor(sand.all$Coverage, levels = c("Native range SDM", "Philippine extent SDM"), ordered = T)

#soil silt
#df1 <- dplyr::select(rspe.pts, Species, Coverage, silt)
df2 <- dplyr::select(rspe.ph, Species, Coverage, Scenario, silt)
df3 <- dplyr::select(rspe.nt, Species, Coverage, Scenario, silt)
rspe.silt <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlob.pts, Species, Coverage, silt)
df2 <- dplyr::select(rlob.ph, Species, Coverage, Scenario, silt)
df3 <- dplyr::select(rlob.nt, Species, Coverage, Scenario, silt)
rlob.silt <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlag.pts, Species, Coverage, silt)
df2 <- dplyr::select(rlag.ph, Species, Coverage, Scenario, silt)
df3 <- dplyr::select(rlag.nt, Species, Coverage, Scenario, silt)
rlag.silt <- bind_rows( df2, df3)

silt.all <- bind_rows(rspe.silt, rlob.silt, rlag.silt)
silt.all$Coverage <- factor(silt.all$Coverage, levels = c( "Native range SDM", "Philippine extent SDM"), ordered = T)

#soil clay
#df1 <- dplyr::select(rspe.pts, Species, Coverage, clay)
df2 <- dplyr::select(rspe.ph, Species, Coverage, Scenario, clay)
df3 <- dplyr::select(rspe.nt, Species, Coverage, Scenario, clay)
rspe.clay <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlob.pts, Species, Coverage, clay)
df2 <- dplyr::select(rlob.ph, Species, Coverage, Scenario, clay)
df3 <- dplyr::select(rlob.nt, Species, Coverage, Scenario, clay)
rlob.clay <- bind_rows( df2, df3)

#df1 <- dplyr::select(rlag.pts, Species, Coverage, clay)
df2 <- dplyr::select(rlag.ph, Species, Coverage, Scenario, clay)
df3 <- dplyr::select(rlag.nt, Species, Coverage, Scenario, clay)
rlag.clay <- bind_rows( df2, df3)

clay.all <- bind_rows(rspe.clay, rlob.clay, rlag.clay)
clay.all$Coverage <- factor(clay.all$Coverage, levels = c( "Native range SDM", "Philippine extent SDM"), ordered = T)



#### BOX PLOTS ####

#altitude + bioclim #

altboxgg <- ggplot() + 
  geom_boxplot(aes(y = altitude, x = Species, fill = Coverage), 
               data=alt.all , outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values=c("#f7f7f7", "#f1b6da", "#d01c8b")) +
  theme_ipsum() +
  ylab("Altitude (masl)")
ggsave(altboxgg, file="Plots/Boxplot-alt-Rafflesia_spp_masked_pts_currentSDM.png", width=19.89, height=15, units="cm", dpi=300)

tempboxgg <- ggplot() + 
  geom_boxplot(aes(y = bio1/10, x = Species, fill= Coverage), show.legend = FALSE, 
               data=temp.all , outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values=c("#f7f7f7", "#f1b6da", "#d01c8b")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  #theme_ipsum() +
  ylab("Mean Annual Temperature (°C)")
ggsave(tempboxgg, file="Plots/Boxplot-matemp-Rafflesia_spp_masked_pts_currentSDM.png", width=19.89, height=15, units="cm", dpi=300)

precboxgg <- ggplot() + 
  geom_boxplot(aes(y = bio12, x = Species, fill=Coverage), show.legend = FALSE,
               data=pre.all , outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values=c("#f7f7f7", "#f1b6da", "#d01c8b")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  #theme_ipsum() +
  ylab("Mean Annual Precipitation (mm/yr)")
ggsave(precboxgg, file="Plots/Boxplot-prec-Rafflesia_spp_masked_pts_currentSDM.png", width=19.89, height=15, units="cm", dpi=300)

tsboxgg <- ggplot() + 
  geom_boxplot(aes(y = bio4, x = Species, fill=Coverage), show.legend = FALSE,
               data=ts.all , outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values=c("#f7f7f7", "#f1b6da", "#d01c8b")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  #theme_ipsum() +
  ylab("Temperature Seasonality (s.d.)")
ggsave(precboxgg, file="Plots/Boxplot-tempsea-Rafflesia_spp_masked_pts_currentSDM.png", width=19.89, height=15, units="cm", dpi=300)

psboxgg <- ggplot() + 
  geom_boxplot(aes(y = bio15, x = Species, fill=Coverage), show.legend = FALSE,
               data=ps.all , outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values=c("#f7f7f7", "#f1b6da", "#d01c8b")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  #theme_ipsum() +
  ylab("Precipitation Seasonality (c.v.)")
ggsave(precboxgg, file="Plots/Boxplot-precsea-Rafflesia_spp_masked_pts_currentSDM.png", width=19.89, height=15, units="cm", dpi=300)


# soil parameters #
phboxgg <- ggplot() + 
  geom_boxplot(aes(y = soilph/10, x = Species, fill=Coverage), show.legend = FALSE,
               data=sph.all , outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values=c("#f7f7f7", "#f1b6da", "#d01c8b")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  #theme_ipsum() +
  ylab("soil pH")
ggsave(phboxgg, file="Plots/Boxplot-soilpH-Rafflesia_spp_masked_pts_currentSDM.png", width=19.89, height=15, units="cm", dpi=300)


cecboxgg <- ggplot() + 
  geom_boxplot(aes(y = cec, x = Species, fill=Coverage), 
               data=cec.all , outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values=c("#f7f7f7", "#f1b6da", "#d01c8b")) +
  theme_ipsum() +
  ylab("soil CEC (cmol/kg)")
ggsave(cecboxgg, file="Plots/Boxplot-soilCEC-Rafflesia_spp_masked_pts_currentSDM.png", width=19.89, height=15, units="cm", dpi=300)


bdboxgg <- ggplot() + 
  geom_boxplot(aes(y = bulk_density, x = Species, fill=Coverage), 
               data=bd.all , outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values=c("#f7f7f7", "#f1b6da", "#d01c8b")) +
  #theme_ipsum() +
  ylab("soil bulk density (kg/cubic-m)")
ggsave(bdboxgg, file="Plots/Boxplot-soilbulkdens-Rafflesia_spp_masked_pts_currentSDM.png", width=19.89, height=15, units="cm", dpi=300)


sandboxgg <- ggplot() + 
  geom_boxplot(aes(y = sand, x = Species, fill=Coverage), #show.legend = FALSE,
               data=sand.all , outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values=c("#f7f7f7", "#f1b6da", "#d01c8b")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  #theme_ipsum() +
  ylab("soil sand content (%)")
ggsave(sandboxgg, file="Plots/Boxplot-soilsand-Rafflesia_spp_masked_pts_currentSDM.png", width=19.89, height=15, units="cm", dpi=300)


siltboxgg <- ggplot() + 
  geom_boxplot(aes(y = silt, x = Species, fill=Coverage), show.legend = FALSE,
               data=silt.all , outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values=c("#f7f7f7", "#f1b6da", "#d01c8b")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  #theme_ipsum() +
  ylab("soil silt content (%)")
ggsave(siltboxgg, file="Plots/Boxplot-soilsilt-Rafflesia_spp_masked_pts_currentSDM.png", width=19.89, height=15, units="cm", dpi=300)


clayboxgg <- ggplot() + 
  geom_boxplot(aes(y = clay, x = Species, fill=Coverage), show.legend = FALSE,
               data=clay.all , outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values=c("#f7f7f7", "#f1b6da", "#d01c8b")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  ylab("soil clay content (%)")
ggsave(clayboxgg, file="Plots/Boxplot-soilclay-Rafflesia_spp_masked_pts_currentSDM.png", width=19.89, height=15, units="cm", dpi=300)



# plot together
grid.arrange(tempboxgg, precboxgg, ncol=2, nrow=1)
grid.arrange(tsboxgg, psboxgg, ncol=2, nrow=1)
grid.arrange(phboxgg, sandboxgg, ncol=2, nrow=1)
grid.arrange(siltboxgg, clayboxgg, ncol=2, nrow=1)

mylegend <- cowplot::get_legend(altboxgg)
grid.newpage()
grid.draw(mylegend)



#### BOXPLOTS FOR REVISION ####
### 9 Sep 2023 ###

# Check if soil conditions are changing under different climate scenarios #
# Only plotting soils since climate is the changing factor #

# Define the color gradient for species and climate scenarios
species_colors <- c("#d9d9d9", "#bdbdbd", "#4d4d4d") # Add more colors if needed
scenario_colors <- c("#f7f7f7",  "#f1b6da", "#d01c8b")

# create ggplots
# soil pH

# Split the Coverage column into Native island and entire PH
sph.all.nt <- filter(sph.all, Coverage=="Native range SDM")
sph.all.ph <- filter(sph.all, Coverage=="Philippine extent SDM")

# Generate the plots
phboxgg1 <- ggplot(data = sph.all.nt, aes(y = soilph/10, x = Species, fill = Scenario)) +
    geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
    scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    ylab("soil pH") +
    facet_wrap(~ Species, scales = "free_x")+
    coord_cartesian(ylim = c(4.2, 6.3))  # Set the desired y-axis limits
  
phboxgg2 <- ggplot(data = sph.all.ph, aes(y = soilph/10, x = Species, fill = Scenario)) +
    geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
    scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    ylab("soil pH") +
    facet_wrap(~ Species, scales = "free_x") +
    coord_cartesian(ylim = c(4.2, 6.3))  # Set the desired y-axis limits
  
# Add labels and combine the plots
combined_plot <- cowplot::plot_grid(phboxgg1 + ggtitle("Native Range"), phboxgg2 + ggtitle("Entire Philippines"), nrow = 1, align = "h", labels = c("A", "B"), hjust = -2, vjust = 1.5)
  

# Save the plot
ggsave(combined_plot, file = "Plots/Revisions/Boxplot-soilpH-Rafflesia_spp_masked_all_climates.png",
       width = 45, height = 15, units = "cm", dpi = 300)


# soil CEC
# Split the Coverage column into Native island and entire PH
cec.all.nt <- filter(cec.all, Coverage=="Native range SDM")
cec.all.ph <- filter(cec.all, Coverage=="Philippine extent SDM")

# Generate the plots
cecboxgg1 <- ggplot(data = cec.all.nt, aes(y = cec, x = Species, fill = Scenario)) +
  geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
 theme_minimal() +
  theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("soil CEC (cmol/kg)") +
  facet_wrap(~ Species, scales = "free_x") +
  coord_cartesian(ylim = c(5, 55))  # Set the desired y-axis limits

cecboxgg2 <- ggplot(data = cec.all.ph, aes(y = cec, x = Species, fill = Scenario)) +
  geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("soil CEC (cmol/kg)") +
  facet_wrap(~ Species, scales = "free_x") +
  coord_cartesian(ylim = c(5, 55))  # Set the desired y-axis limits


# Add labels and combine the plots
combined_plot <- cowplot::plot_grid(cecboxgg1 + ggtitle("Native Range"), cecboxgg2 + ggtitle("Entire Philippines"), nrow = 1, align = "h", labels = c("A", "B"), hjust = -2, vjust = 1.5)


# Save the plot
ggsave(combined_plot, file = "Plots/Revisions/Boxplot-soilCEC-Rafflesia_spp_masked_all_climates.png",
       width = 45, height = 15, units = "cm", dpi = 300)


# Bulk Density #
# Split the Coverage column into Native island and entire PH
bd.all.nt <- filter(bd.all, Coverage=="Native range SDM")
bd.all.ph <- filter(bd.all, Coverage=="Philippine extent SDM")

# Generate the plots
bdboxgg1 <- ggplot(data = bd.all.nt, aes(y = bulk_density, x = Species, fill = Scenario)) +
  geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("soil bulk density (kg/cubic-m)") +
  facet_wrap(~ Species, scales = "free_x")+
  coord_cartesian(ylim = c(550, 1400))  # Set the desired y-axis limits

bdboxgg2 <- ggplot(data = bd.all.ph, aes(y = bulk_density, x = Species, fill = Scenario)) +
  geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("soil bulk density (kg/cubic-m)") +
  facet_wrap(~ Species, scales = "free_x")+
  coord_cartesian(ylim = c(550, 1400))  # Set the desired y-axis limits


# Add labels and combine the plots
combined_plot <- cowplot::plot_grid(bdboxgg1 + ggtitle("Native Range"), bdboxgg2 + ggtitle("Entire Philippines"), nrow = 1, align = "h", labels = c("A", "B"), hjust = -2, vjust = 1.5)

# Save the plot
ggsave(combined_plot, file = "Plots/Revisions/Boxplot-soilbulkdens-Rafflesia_spp_masked_all_climates.png",
       width = 45, height = 15, units = "cm", dpi = 300)


#Sand 
# Split the Coverage column into Native island and entire PH
sand.all.nt <- filter(sand.all, Coverage=="Native range SDM")
sand.all.ph <- filter(sand.all, Coverage=="Philippine extent SDM")

# Generate the plots
sandboxgg1 <- ggplot(data = sand.all.nt, aes(y = sand, x = Species, fill = Scenario)) +
  geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("soil sand content (%)") +
  facet_wrap(~ Species, scales = "free_x")+
  coord_cartesian(ylim = c(20, 65))  # Set the desired y-axis limits

sandboxgg2 <- ggplot(data = sand.all.ph, aes(y = sand, x = Species, fill = Scenario)) +
  geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("soil sand content (%)") +
  facet_wrap(~ Species, scales = "free_x")+
  coord_cartesian(ylim = c(20, 65))  # Set the desired y-axis limits


# Add labels and combine the plots
combined_plot <- cowplot::plot_grid(sandboxgg1 + ggtitle("Native Range"), sandboxgg2 + ggtitle("Entire Philippines"), nrow = 1, align = "h", labels = c("A", "B"), hjust = -2, vjust = 1.5)

# Save the plot
ggsave(combined_plot, file = "Plots/Revisions/Boxplot-soilsand-Rafflesia_spp_masked_all_climates.png",
       width = 45, height = 15, units = "cm", dpi = 300)


#Silt 
# Split the Coverage column into Native island and entire PH
silt.all.nt <- filter(silt.all, Coverage=="Native range SDM")
silt.all.ph <- filter(silt.all, Coverage=="Philippine extent SDM")

# Generate the plots
siltboxgg1 <- ggplot(data = silt.all.nt, aes(y = silt, x = Species, fill = Scenario)) +
  geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("soil silt content (%)") +
  facet_wrap(~ Species, scales = "free_x")+
  coord_cartesian(ylim = c(15, 50))  # Set the desired y-axis limits

siltboxgg2 <- ggplot(data = silt.all.ph, aes(y = silt, x = Species, fill = Scenario)) +
  geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("soil silt content (%)") +
  facet_wrap(~ Species, scales = "free_x")+
  coord_cartesian(ylim = c(15, 50))  # Set the desired y-axis limits


# Add labels and combine the plots
combined_plot <- cowplot::plot_grid(siltboxgg1 + ggtitle("Native Range"), siltboxgg2 + ggtitle("Entire Philippines"), nrow = 1, align = "h", labels = c("A", "B"), hjust = -2, vjust = 1.5)

# Save the plot
ggsave(combined_plot, file = "Plots/Revisions/Boxplot-soilsilt-Rafflesia_spp_masked_all_climates.png",
       width = 45, height = 15, units = "cm", dpi = 300)

# Clay 
# Split the Coverage column into Native island and entire PH
clay.all.nt <- filter(clay.all, Coverage=="Native range SDM")
clay.all.ph <- filter(clay.all, Coverage=="Philippine extent SDM")

# Generate the plots
clayboxgg1 <- ggplot(data = clay.all.nt, aes(y = clay, x = Species, fill = Scenario)) +
  geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("soil clay content (%)") +
  facet_wrap(~ Species, scales = "free_x")+
  coord_cartesian(ylim = c(10, 50))  # Set the desired y-axis limits

clayboxgg2 <- ggplot(data = clay.all.ph, aes(y = clay, x = Species, fill = Scenario)) +
  geom_boxplot(show.legend = T, outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12, face = "italic"), # Add italic font to x-axis labels
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ylab("soil clay content (%)") +
  facet_wrap(~ Species, scales = "free_x")+
  coord_cartesian(ylim = c(10, 50))  # Set the desired y-axis limits


# Add labels and combine the plots
combined_plot <- cowplot::plot_grid(clayboxgg1 + ggtitle("Native Range"), clayboxgg2 + ggtitle("Entire Philippines"), nrow = 1, align = "h", labels = c("A", "B"), hjust = -2, vjust = 1.5)

# Save the plot
ggsave(combined_plot, file = "Plots/Revisions/Boxplot-soilclay-Rafflesia_spp_masked_all_climates.png",
       width = 45, height = 15, units = "cm", dpi = 300)



