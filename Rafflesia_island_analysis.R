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


#### Compute suitable areas percentage for the current era per species ####