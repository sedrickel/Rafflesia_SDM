library(raster)
library(rgdal)
require(tidyverse)
library(xlsx)

# Island-level analysis of Rafflesia SDM results #

#### read in the SDM result rasters #####
#### current  
# entire PH

setwd("D:/PROJECTS/Rafflesia_SDM/A1/MASKED")
#masked, A1
rlag.cur.a1 <- raster("current/R_lagascae_01_currentEnv_bin_ave_masked.tif")
rlob.cur.a1 <- raster("current/R_lobata_01_currentEnv_bin_ave_masked.tif")
rspe.cur.a1 <- raster("current/R_speciosa_01_currentEnv_bin_ave_masked.tif")

#RCP 4.5
rlag.r45.a1 <- raster("future/rcp45/R_lagascae_01_2070_RCP45_ensem_bin_masked.tif")
rlob.r45.a1 <- raster("future/rcp45/R_lobata_01_2070_RCP45_ensem_bin_masked.tif")
rspe.r45.a1 <- raster("future/rcp45/R_speciosa_01_2070_RCP45_ensem_bin_masked.tif")

#RCP 8.5
rlag.r85.a1 <- raster("future/rcp85/R_lagascae_01_2070_RCP85_ensem_bin_masked.tif")
rlob.r85.a1 <- raster("future/rcp85/R_lobata_01_2070_RCP85_ensem_bin_masked.tif")
rspe.r85.a1 <- raster("future/rcp85/R_speciosa_01_2070_RCP85_ensem_bin_masked.tif")


#masked, A2
setwd("D:/PROJECTS/Rafflesia_SDM/A2/MASKED/")
rlag.cur.a2 <- raster("current/R_lagascae_A2_currentEnv_masked.tif")
rlob.cur.a2 <- raster("current/R_lobata_A2_currentEnv_masked.tif")
rspe.cur.a2 <- raster("current/R_speciosa_A2_currentEnv_masked.tif")

#RCP 4.5
rlag.r45.a2 <- raster("future/RCP4.5/R_lagascae_A2_RCP45_2070_masked.tif")
rlob.r45.a2 <- raster("future/RCP4.5/R_lobata_A2_RCP45_2070_masked.tif")
rspe.r45.a2 <- raster("future/RCP4.5/R_speciosa_A2_RCP45_2070_masked.tif")

#RCP 8.5
rlag.r85.a2 <- raster("future/RCP8.5/R_lagascae_A2_RCP85_2070_masked.tif")
rlob.r85.a2 <- raster("future/RCP8.5/R_lobata_A2_RCP85_2070_masked.tif")
rspe.r85.a2 <- raster("future/RCP8.5/R_speciosa_A2_RCP85_2070_masked.tif")

#island shapefiles 
setwd("D:/PROJECTS/Rafflesia_SDM/Island_analysis")

luzon <- readOGR(dsn=getwd(), layer = "Luzon_mainland")
panneg <- readOGR(dsn=getwd(), layer = "Panay_Negros_merge")
panay <- readOGR(dsn=getwd(), layer = "Panay_mainland")

windows()
par(mfrow=c(1,3))
plot(luzon)
plot(panneg)
plot(panay)

#### Crop the rasters as per species' native range ####

#Rlagascae - Mainland Luzon and Bicol
#A1
rlag.cur.a1 <- crop(rlag.cur.a1, luzon)
rlag.cur.a1 <- mask(rlag.cur.a1, luzon)

rlag.r45.a1 <- crop(rlag.r45.a1, luzon)
rlag.r45.a1 <- mask(rlag.r45.a1, luzon)

rlag.r85.a1 <- crop(rlag.r85.a1, luzon)
rlag.r85.a1 <- mask(rlag.r85.a1, luzon)

#A2
rlag.cur.a2 <- crop(rlag.cur.a2, luzon)
rlag.cur.a2 <- mask(rlag.cur.a2, luzon)

rlag.r45.a2 <- crop(rlag.r45.a2, luzon)
rlag.r45.a2 <- mask(rlag.r45.a2, luzon)

rlag.r85.a2 <- crop(rlag.r85.a2, luzon)
rlag.r85.a2 <- mask(rlag.r85.a2, luzon)

#Rlobata - Panay
#A1
rlob.cur.a1 <- crop(rlob.cur.a1, panay)
rlob.cur.a1 <- mask(rlob.cur.a1, panay)

rlob.r45.a1 <- crop(rlob.r45.a1, panay)
rlob.r45.a1 <- mask(rlob.r45.a1, panay)

rlob.r85.a1 <- crop(rlob.r85.a1, panay)
rlob.r85.a1 <- mask(rlob.r85.a1, panay)

#A2
rlob.cur.a2 <- crop(rlob.cur.a2, panay)
rlob.cur.a2 <- mask(rlob.cur.a2, panay)

rlob.r45.a2 <- crop(rlob.r45.a2, panay)
rlob.r45.a2 <- mask(rlob.r45.a2, panay)

rlob.r85.a2 <- crop(rlob.r85.a2, panay)
rlob.r85.a2 <- mask(rlob.r85.a2, panay)

#Rspeciosa - Panay and Negros
#A1
rspe.cur.a1 <- crop(rspe.cur.a1, panneg)
rspe.cur.a1 <- mask(rspe.cur.a1, panneg)

rspe.r45.a1 <- crop(rspe.r45.a1, panneg)
rspe.r45.a1 <- mask(rspe.r45.a1, panneg)

rspe.r85.a1 <- crop(rspe.r85.a1, panneg)
rspe.r85.a1 <- mask(rspe.r85.a1, panneg)

#A2
rspe.cur.a2 <- crop(rspe.cur.a2, panneg)
rspe.cur.a2 <- mask(rspe.cur.a2, panneg)

rspe.r45.a2 <- crop(rspe.r45.a2, panneg)
rspe.r45.a2 <- mask(rspe.r45.a2, panneg)

rspe.r85.a2 <- crop(rspe.r85.a2, panneg)
rspe.r85.a2 <- mask(rspe.r85.a2, panneg)

#plot to check
windows()
#R lagascae
par(mfrow=c(2,3))
plot(rlag.cur.a1)
plot(rlag.r45.a1)
plot(rlag.r85.a1)
plot(rlag.cur.a2)
plot(rlag.r45.a2)
plot(rlag.r85.a2)


#R lobata
par(mfrow=c(2,3))
plot(rlob.cur.a1)
plot(rlob.r45.a1)
plot(rlob.r85.a1)
plot(rlob.cur.a2)
plot(rlob.r45.a2)
plot(rlob.r85.a2)

#R speciosa
par(mfrow=c(2,3))
plot(rspe.cur.a1)
plot(rspe.r45.a1)
plot(rspe.r85.a1)
plot(rspe.cur.a2)
plot(rspe.r45.a2)
plot(rspe.r85.a2)


#Write it out 
setwd("D:/PROJECTS/Rafflesia_SDM/A1/MASKED")
#masked, A1
writeRaster(rlag.cur.a1, filename="current/R_lagascae_01_currentEnv_bin_ave_masked_island.tif", datatype="INT1U")
writeRaster(rlob.cur.a1, filename="current/R_lobata_01_currentEnv_bin_ave_masked_island.tif", datatype="INT1U")
writeRaster(rspe.cur.a1, filename="current/R_speciosa_01_currentEnv_bin_ave_masked_island.tif", datatype="INT1U")

#RCP 4.5
writeRaster(rlag.r45.a1, filename="future/rcp45/R_lagascae_01_2070_RCP45_ensem_bin_masked_island.tif", datatype="INT1U")
writeRaster(rlob.r45.a1, filename="future/rcp45/R_lobata_01_2070_RCP45_ensem_bin_masked_island.tif", datatype="INT1U")
writeRaster(rspe.r45.a1, filename="future/rcp45/R_speciosa_01_2070_RCP45_ensem_bin_masked_island.tif", datatype="INT1U")

#RCP 8.5
writeRaster(rlag.r85.a1, filename="future/rcp85/R_lagascae_01_2070_RCP85_ensem_bin_masked_island.tif", datatype="INT1U")
writeRaster(rlob.r85.a1, filename="future/rcp85/R_lobata_01_2070_RCP85_ensem_bin_masked_island.tif", datatype="INT1U")
writeRaster(rspe.r85.a1, filename="future/rcp85/R_speciosa_01_2070_RCP85_ensem_bin_masked_island.tif", datatype="INT1U")


#masked, A2
setwd("D:/PROJECTS/Rafflesia_SDM/A2/MASKED/")
writeRaster(rlag.cur.a2, filename="current/R_lagascae_A2_currentEnv_masked_island.tif", datatype="INT1U")
writeRaster(rlob.cur.a2, filename="current/R_lobata_A2_currentEnv_masked_island.tif", datatype="INT1U")
writeRaster(rspe.cur.a2, filename="current/R_speciosa_A2_currentEnv_masked_island.tif", datatype="INT1U")

#RCP 4.5
writeRaster(rlag.r45.a2, filename="future/RCP4.5/R_lagascae_A2_RCP45_2070_masked_island.tif", datatype="INT1U")
writeRaster(rlob.r45.a2, filename="future/RCP4.5/R_lobata_A2_RCP45_2070_masked_island.tif", datatype="INT1U")
writeRaster(rspe.r45.a2, filename="future/RCP4.5/R_speciosa_A2_RCP45_2070_masked_island.tif", datatype="INT1U")

#RCP 8.5
writeRaster(rlag.r85.a2, filename="future/RCP8.5/R_lagascae_A2_RCP85_2070_masked_island.tif", datatype="INT1U")
writeRaster(rlob.r85.a2, filename="future/RCP8.5/R_lobata_A2_RCP85_2070_masked_island.tif", datatype="INT1U")
writeRaster(rspe.r85.a2, filename="future/RCP8.5/R_speciosa_A2_RCP85_2070_masked_island.tif", datatype="INT1U")


####COUNT NO OF PIXELS CHANGE ####
setwd("D:/PROJECTS/Rafflesia_SDM/Island_analysis")

sp.name <- "rlob" #rlob, rspe

current.a1.df <- raster::as.data.frame(get(paste0(sp.name,".cur.a1")), na.rm = T)
colnames(current.a1.df) <- "x"
current.a1.sum <- count(current.a1.df, x)

rcp45.a1.df <- as.data.frame(get(paste0(sp.name,".r45.a1")), na.rm = T)
colnames(rcp45.a1.df) <- "x"
rcp45.a1.sum <- count(rcp45.a1.df, x)

rcp85.a1.df <- as.data.frame(get(paste0(sp.name,".r85.a1")), na.rm = T)
colnames(rcp85.a1.df) <- "x"
rcp85.a1.sum <- count(rcp85.a1.df, x)

current.a2.df <- raster::as.data.frame(get(paste0(sp.name,".cur.a1")), na.rm = T)
colnames(current.a2.df) <- "x"
current.a2.sum <- count(current.a2.df, x)

rcp45.a2.df <- as.data.frame(get(paste0(sp.name,".r45.a2")), na.rm = T)
colnames(rcp45.a2.df) <- "x"
rcp45.a2.sum <- count(rcp45.a2.df, x)

rcp85.a2.df <- as.data.frame(get(paste0(sp.name,".r85.a2")), na.rm = T)
colnames(rcp85.a2.df) <- "x"
rcp85.a2.sum <- count(rcp85.a2.df, x)

change.summ <- as.data.frame(cbind(current.a1.sum$x, current.a1.sum$n, current.a2.sum$n, rcp45.a1.sum$n, rcp85.a1.sum$n, rcp45.a2.sum$n, rcp85.a2.sum$n))
colnames(change.summ) <- c("State", "Current-A1", "Current-A2", "RCP4.5-2070-A1", "RCP8.5-2070-A1", "RCP4.5-2070-A2", "RCP8.5-2070-A2")

change.summ[3,1] <- "%Suitable Area"
change.summ[3,2] <- change.summ[2,2]/(change.summ[1,2] + change.summ[2,2])*100
change.summ[3,3] <- change.summ[2,3]/(change.summ[1,3] + change.summ[2,3])*100
change.summ[3,4] <- change.summ[2,4]/(change.summ[1,4] + change.summ[2,4])*100
change.summ[3,5] <- change.summ[2,5]/(change.summ[1,5] + change.summ[2,5])*100
change.summ[3,6] <- change.summ[2,6]/(change.summ[1,6] + change.summ[2,6])*100
change.summ[3,7] <- change.summ[2,7]/(change.summ[1,7] + change.summ[2,7])*100
change.summ

write.csv(change.summ, file= paste0(sp.name,"_pixel_changes.csv"), row.names = F)


# Extract data for boxplots #
# this assumes that envi data has been read in from main Rafflesia_SDM codes
# uses A1 - masked data only #

#speciosa
spe <- as.data.frame(rspe.cur.a1, na.rm = T, xy=T)
colnames(spe) <- c("lon", "lat", "bin")
spe <- filter(spe, bin==1)
spe <- cbind(spe$lon, spe$lat)

spe.df <- raster::extract(envdata, spe, method = "simple", na.rm = T)
spe.df <- as.data.frame(spe.df)
spe.df$Scenario <- rep("Current", dim(spe.df)[1])

#RCP4.5
spe45 <- as.data.frame(rspe.r45.a1, na.rm = T, xy=T)
colnames(spe45) <- c("lon", "lat", "bin")
spe45 <- filter(spe45, bin==1)
spe45 <- cbind(spe45$lon, spe45$lat)

spe45.df <- raster::extract(envdata, spe45, method = "simple", na.rm = T)
spe45.df <- as.data.frame(spe45.df)
spe45.df$Scenario <- rep("RCP4.5", dim(spe45.df)[1])


#RCP 8.5
spe85 <- as.data.frame(rspe.r85.a1, na.rm = T, xy=T)
colnames(spe85) <- c("lon", "lat", "bin")
spe85 <- filter(spe85, bin==1)
spe85 <- cbind(spe85$lon, spe85$lat)

spe85.df <- raster::extract(envdata, spe85, method = "simple", na.rm = T)
spe85.df <- as.data.frame(spe85.df)
spe85.df$Scenario <- rep("RCP8.5", dim(spe85.df)[1])


#combine all data from one species

allspe <- rbind(spe.df, spe45.df, spe85.df)
allspe$Species <- rep("Rspe", dim(allspe)[1]) #change for each species

summary(allspe)

## R. lagascae
#current
lag <- as.data.frame(rlag.cur.a1, na.rm = T, xy=T)
colnames(lag) <- c("lon", "lat", "bin")
lag <- filter(lag, bin==1)
lag <- cbind(lag$lon, lag$lat)

lag.df <- raster::extract(envdata, lag, method = "simple", na.rm = T)
lag.df <- as.data.frame(lag.df)
lag.df$Scenario <- rep("Current", dim(lag.df)[1])

#RCP4.5
lag45 <- as.data.frame(rlag.r45.a1, na.rm = T, xy=T)
colnames(lag45) <- c("lon", "lat", "bin")
lag45 <- filter(lag45, bin==1)
lag45 <- cbind(lag45$lon, lag45$lat)

lag45.df <- raster::extract(envdata, lag45, method = "simple", na.rm = T)
lag45.df <- as.data.frame(lag45.df)
lag45.df$Scenario <- rep("RCP4.5", dim(lag45.df)[1])


#RCP 8.5
lag85 <- as.data.frame(rlag.r85.a1, na.rm = T, xy=T)
colnames(lag85) <- c("lon", "lat", "bin")
lag85 <- filter(lag85, bin==1)
lag85 <- cbind(lag85$lon, lag85$lat)

lag85.df <- raster::extract(envdata, lag85, method = "simple", na.rm = T)
lag85.df <- as.data.frame(lag85.df)
lag85.df$Scenario <- rep("RCP8.5", dim(lag85.df)[1])


#combine all data from one species

alllag <- rbind(lag.df, lag45.df, lag85.df)
alllag$Species <- rep("Rlag", dim(alllag)[1]) #change for each species
summary(alllag)

## R. lobata
#current
lob <- as.data.frame(rlob.cur.a1, na.rm = T, xy=T)
colnames(lob) <- c("lon", "lat", "bin")
lob <- filter(lob, bin==1)
lob <- cbind(lob$lon, lob$lat)

lob.df <- raster::extract(envdata, lob, method = "simple", na.rm = T)
lob.df <- as.data.frame(lob.df)
lob.df$Scenario <- rep("Current", dim(lob.df)[1])

#RCP4.5
lob45 <- as.data.frame(rlob.r45.a1, na.rm = T, xy=T)
colnames(lob45) <- c("lon", "lat", "bin")
lob45 <- filter(lob45, bin==1)
lob45 <- cbind(lob45$lon, lob45$lat)

lob45.df <- raster::extract(envdata, lob45, method = "simple", na.rm = T)
lob45.df <- as.data.frame(lob45.df)
lob45.df$Scenario <- rep("RCP4.5", dim(lob45.df)[1])


#RCP 8.5
lob85 <- as.data.frame(rlob.r85.a1, na.rm = T, xy=T)
colnames(lob85) <- c("lon", "lat", "bin")
lob85 <- filter(lob85, bin==1)
lob85 <- cbind(lob85$lon, lob85$lat)

lob85.df <- raster::extract(envdata, lob85, method = "simple", na.rm = T)
lob85.df <- as.data.frame(lob85.df)
lob85.df$Scenario <- rep("RCP8.5", dim(lob85.df)[1])

#combine all data from one species

alllob <- rbind(lob.df, lob45.df, lob85.df)
alllob$Species <- rep("Rlob", dim(alllob)[1]) #change for each species

summary(alllob)

#remove NAs
allspe <- na.omit(allspe)
alllag <- na.omit(alllag)
alllob <- na.omit(alllob)

#combine all species df
all.raf <- rbind(allspe, alllag, alllob)
all.raf$Species <- as.factor(all.raf$Species)
all.raf$Scenario <- as.factor(all.raf$Scenario)
summary(all.raf)

write.csv(all.raf, file = "All_Rafflesia_spp_a1_suitable_masked_native_ranges_envdata.csv", row.names = F)

