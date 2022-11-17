library(raster)
require(tidyverse)

setwd("D:PROJECTS/Rafflesia_SDM")

#read in species binary rasters

a1ras <- raster("spec_bin_a1.tif")
a2ras <- raster("spec_a2.tif")

windows()
plot(a1ras)
plot(a2ras)

#convert A2 to binary
a2ras.b <- a2ras * 0
plot(a2ras.b)

a2ras.b[a2ras == 3] <- 1

#plot together

par(mfrow=c(1,3))
plot(a1ras, main = "A1 binary")
plot(a2ras.b, main = "A2 binary")
plot(a2ras, main = "A2 summed")

#check overlaps visually
sumras <- a1ras + a2ras.b
par(mfrow=c(1,1))
plot(sumras)

#convert to dataframe for quantitative comparison
a1ras.df <- as.data.frame(a1ras, xy = T, na.rm = T)
summary(a1ras.df)
colnames(a1ras.df)[3] <- "A1"

a2ras.b.df <- as.data.frame(a2ras.b, xy = T, na.rm = T)
summary(a2ras.b.df)
colnames(a2ras.b.df)[3] <- "A2"

#combine both df's
comb.df <- left_join(a1ras.df, a2ras.b.df, by = c("x", "y"))

#create new column that compares A1 vs A2
comb.df <- comb.df  %>% 
  mutate(
    T = case_when(
      A1 == A2 ~ 1,
      A1 != A2 ~ 0)
  )

#sum the columns of the data frame
#if T = total number of pixels, then the two species binary rasters completely overlap
#if T < total number of pixels, then the two species binary rasters do NOT completely overlap
colSums(comb.df)


#if there is no complete overlap (<100%), you can divide the sum of column T over the total number of pixels in the raster (in this case, 347947)
