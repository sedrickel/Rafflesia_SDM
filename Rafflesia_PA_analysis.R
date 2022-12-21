library(raster)
library(rgdal)
require(tidyverse)
library(xlsx)

setwd("D:/PROJECTS/Rafflesia_SDM/PA_analysis")

#get PH shapefile
phl.shp <- getData('GADM', country = 'PHL', level=0)
extent(phl.shp)

#very important! run and DO NOT change anything. this contains the raster specifications: CRS projection, extent, and resolution
x <-raster(xmn=116.9283  ,xmx=126.6053 ,ymn=4.58694 ,ymx=21.07014 ,res=0.0083333333,crs="+proj=longlat +datum=WGS84")


#read in dissolved PA boundaries 
wdpa.shp <- readOGR(dsn=getwd(), layer = "PH_WDPA_dissolved")
wdpa.shp@data$Id <- as.factor(1) 

#convert shp to raster
wdpa.ras <- rasterize(wdpa.shp, x, field = "Id", filename = "PH_WDPA_dissolved_30s.tif", datatype = "INT2U", format = "GTiff")

#read back in
wdpa.ras <- raster("PH_WDPA_dissolved_30s.tif")
plot(wdpa.ras)

#convert to data frame
wdpa.df <- as.data.frame(wdpa.ras, xy=T, na.rm = T)
summary(wdpa.df)
wdpa.xy <- wdpa.df[1:2]

#### read in current rasters #####
# entire PH
#masked, A1
rlag.cur.a1 <- raster("D:/PROJECTS/Rafflesia_SDM/A1/MASKED/current/R_lagascae_01_currentEnv_bin_ave_masked.tif")
rlob.cur.a1 <- raster("D:/PROJECTS/Rafflesia_SDM/A1/MASKED/current/R_lobata_01_currentEnv_bin_ave_masked.tif")
rspe.cur.a1 <- raster("D:/PROJECTS/Rafflesia_SDM/A1/MASKED/current/R_speciosa_01_currentEnv_bin_ave_masked.tif")

#masked, A2
rlag.cur.a2 <- raster("D:/PROJECTS/Rafflesia_SDM/A2/MASKED/current/R_lagascae_A2_currentEnv_masked.tif")
rlob.cur.a2 <- raster("D:/PROJECTS/Rafflesia_SDM/A2/MASKED/current/R_lobata_A2_currentEnv_masked.tif")
rspe.cur.a2 <- raster("D:/PROJECTS/Rafflesia_SDM/A2/MASKED/current/R_speciosa_A2_currentEnv_masked.tif")

#island-limited only



#convert to data frame
#A1
rlag.cur.a1.df <- as.data.frame(rlag.cur.a1, xy = F, na.rm = T)
colnames(rlag.cur.a1.df) <- c("Suitability")
rlag.cur.a1.df <- count(rlag.cur.a1.df, Suitability )

rlob.cur.a1.df <- as.data.frame(rlob.cur.a1, xy = F, na.rm = T)
colnames(rlob.cur.a1.df) <- c("Suitability")
rlob.cur.a1.df <- count(rlob.cur.a1.df, Suitability )

rspe.cur.a1.df <- as.data.frame(rspe.cur.a1, xy = F, na.rm = T)
colnames(rspe.cur.a1.df) <- c("Suitability")
rspe.cur.a1.df <- count(rspe.cur.a1.df, Suitability )

raff.cur.a1.df <- cbind(rlag.cur.a1.df, rlob.cur.a1.df[2], rspe.cur.a1.df[2])
raff.cur.a1.df <- as.data.frame(t(raff.cur.a1.df))
raff.cur.a1.df <- raff.cur.a1.df[-1,]
rownames(raff.cur.a1.df) <- c("Rlag", "Rlob", "Rspe")
colnames(raff.cur.a1.df) <- c("unsuitph", "suitph")

#A2
rlag.cur.a2.df <- as.data.frame(rlag.cur.a2, xy = F, na.rm = T)
colnames(rlag.cur.a2.df) <- c("Suitability")
rlag.cur.a2.df <- count(rlag.cur.a2.df, Suitability )

rlob.cur.a2.df <- as.data.frame(rlob.cur.a2, xy = F, na.rm = T)
colnames(rlob.cur.a2.df) <- c("Suitability")
rlob.cur.a2.df <- count(rlob.cur.a2.df, Suitability )

rspe.cur.a2.df <- as.data.frame(rspe.cur.a2, xy = F, na.rm = T)
colnames(rspe.cur.a2.df) <- c("Suitability")
rspe.cur.a2.df <- count(rspe.cur.a2.df, Suitability )

raff.cur.a2.df <- cbind(rlag.cur.a2.df, rlob.cur.a2.df[2], rspe.cur.a2.df[2])
raff.cur.a2.df <- as.data.frame(t(raff.cur.a2.df))
raff.cur.a2.df <- raff.cur.a2.df[-1,]
rownames(raff.cur.a2.df) <- c("Rlag", "Rlob", "Rspe")
colnames(raff.cur.a2.df) <- c("unsuitph", "suitph")


#### Extract information from binary rasters ####

#A1
#R lagascae
rlag.cur.a1.pa <- raster::extract(rlag.cur.a1, wdpa.xy, df = T, na.rm = T)
summary(rlag.cur.a1.pa)
rlag.cur.a1.pa <- na.omit(rlag.cur.a1.pa)
colnames(rlag.cur.a1.pa) <- c("ID", "Suitability")
rlag.cur.a1.pa2 <- count(rlag.cur.a1.pa, Suitability )

#R lobata
rlob.cur.a1.pa <- raster::extract(rlob.cur.a1, wdpa.xy, df = T, na.rm = T)
summary(rlob.cur.a1.pa)
rlob.cur.a1.pa <- na.omit(rlob.cur.a1.pa)
colnames(rlob.cur.a1.pa) <- c("ID", "Suitability")
rlob.cur.a1.pa2 <- count(rlob.cur.a1.pa, Suitability )

#R speciosa
rspe.cur.a1.pa <- raster::extract(rspe.cur.a1, wdpa.xy, df = T, na.rm = T)
summary(rspe.cur.a1.pa)
rspe.cur.a1.pa <- na.omit(rspe.cur.a1.pa)
colnames(rspe.cur.a1.pa) <- c("ID", "Suitability")
rspe.cur.a1.pa2 <- count(rspe.cur.a1.pa, Suitability )

#A2
#R lagascae
rlag.cur.a2.pa <- raster::extract(rlag.cur.a2, wdpa.xy, df = T, na.rm = T)
summary(rlag.cur.a2.pa)
rlag.cur.a2.pa <- na.omit(rlag.cur.a2.pa)
colnames(rlag.cur.a2.pa) <- c("ID", "Suitability")
rlag.cur.a2.pa2 <- count(rlag.cur.a2.pa, Suitability )

#R lobata
rlob.cur.a2.pa <- raster::extract(rlob.cur.a2, wdpa.xy, df = T, na.rm = T)
summary(rlob.cur.a2.pa)
rlob.cur.a2.pa <- na.omit(rlob.cur.a2.pa)
colnames(rlob.cur.a2.pa) <- c("ID", "Suitability")
rlob.cur.a2.pa2 <- count(rlob.cur.a2.pa, Suitability )

#R speciosa
rspe.cur.a2.pa <- raster::extract(rspe.cur.a2, wdpa.xy, df = T, na.rm = T)
summary(rspe.cur.a2.pa)
rspe.cur.a2.pa <- na.omit(rspe.cur.a2.pa)
colnames(rspe.cur.a2.pa) <- c("ID", "Suitability")
rspe.cur.a2.pa2 <- count(rspe.cur.a2.pa, Suitability )

#### Combine in one df and compute percentage protection ####
#A1
all.cur.pa2 <- cbind(rlag.cur.a1.pa2, rlob.cur.a1.pa2[2], rspe.cur.a1.pa2[2])
all.cur.pa2 <- as.data.frame(t(all.cur.pa2))
all.cur.pa2 <- all.cur.pa2[-1,]

rownames(all.cur.pa2) <- c("Rlag", "Rlob", "Rspe")
colnames(all.cur.pa2) <- c("unsuitpa", "suitpa")

all.cur.pa2 <- cbind(all.cur.pa2, raff.cur.a1.df[2])

all.cur.pa2 <- all.cur.pa2 %>% 
  mutate(percent = suitpa/(unsuitpa + suitpa)*100) %>%
  mutate(protectionsuit = (suitpa/suitph)*100)

write.xlsx(all.cur.pa2, file = "Pr_area_analysis_results.xlsx", sheetName = "A1", col.names = T, row.names = T, append = F)

#A2
all.cur.pa2 <- cbind(rlag.cur.a2.pa2, rlob.cur.a2.pa2[2], rspe.cur.a2.pa2[2])
all.cur.pa2 <- as.data.frame(t(all.cur.pa2))
all.cur.pa2 <- all.cur.pa2[-1,]

rownames(all.cur.pa2) <- c("Rlag", "Rlob", "Rspe")
colnames(all.cur.pa2) <- c("unsuitpa", "suitpa")

all.cur.pa2 <- cbind(all.cur.pa2, raff.cur.a2.df[2])

all.cur.pa2 <- all.cur.pa2 %>% 
  mutate(percent = suitpa/(unsuitpa + suitpa)*100) %>%
  mutate(protectionsuit = (suitpa/suitph)*100)

write.xlsx(all.cur.pa2, file = "Pr_area_analysis_results.xlsx", sheetName = "A2", col.names = T, row.names = T, append = T)


#### Counting number of points inside PA boundaries ####
datapoints <- read.csv("D:/PROJECTS/Rafflesia_SDM/maxent_clean2_Tetrastigma_Rafflesia_30s.csv")
summary(datapoints)
#datapoints$Species <- as.factor(datapoints$Species)
unique(datapoints$Species)

lonlat <- datapoints[,2:3]

spa.df <- raster::extract(wdpa.ras, lonlat, df = T, na.rm = T)
spa.df[is.na(spa.df)] = 0 # replace NA's with 0, meaning they are not inside PAs
summary(spa.df)

datapoints <- cbind(datapoints, spa.df$PH_WDPA_dissolved_30s)
colnames(datapoints)[4] <- "PA.in"
summary(datapoints)

sppcounts <- count(datapoints, Species)

datapoints.summ <- datapoints %>% 
  group_by(Species) %>%
  summarise(PA.ins = sum(PA.in))

datapoints.summ <- cbind(datapoints.summ, sppcounts[2])

datapoints.summ <- datapoints.summ %>% mutate(percentPA = PA.ins/n*100)


write.xlsx(datapoints.summ, file = "Pr_area_analysis_results.xlsx", sheetName = "Points", col.names = T, row.names = T, append = T)
