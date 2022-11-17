require(tidyverse)
library(raster)
library(rasterVis)
library(rJava)
library(usdm)
library(ENMeval)
library(maps)
library(mapdata)
library(dismo)  
library(maptools)
library(jsonlite)
library(rworldmap)
library(RColorBrewer)
library(doParallel)
library(readxl)
library(tmap)
library(ENMTools)
library(forcats)
library(hrbrthemes)
library(viridis)
library(multcompView)
library(visreg)
library(jtools)

registerDoParallel()
getDoParWorkers() 

setwd("D:/PROJECTS/Rafflesia_SDM")

# hello again

### TRUE SKILL STATISTIC TEST FUNCTION #####
#to be able to run this script you need to have told the Maxent model to produce background predictions. If you are running MaxEnt in R this means putting the argument (after "args") "writebackgroundpredictions=true" as true not false. 

#FUNCTION CODE
TSS_calculations <- function (sample_clog, prediction_clog, n, th) {
  
  xx <- sum(sample_clog > th)
  yy <- sum(prediction_clog > th)
  xxx <- sum(sample_clog < th)
  yyy <- sum(prediction_clog < th)
  
  ncount <- sum(xx,yy,xxx,yyy)
  
  overallaccuracy <- (xx + yyy)/ncount 
  sensitivity <- xx / (xx + xxx)
  specificity <- yyy / (yy + yyy)
  tss <- sensitivity + specificity - 1
  
  #kappa calculations
  a <- xx + xxx
  b <- xx + yy
  c <- yy + yyy
  d <- xxx + yyy
  e <- a * b
  f <- c * d
  g <- e + f
  h <- g / (ncount * ncount)
  hup <- overallaccuracy - h
  hdown <- 1 - h
  
  kappa <- hup/hdown
  Po <- (xx + yyy) / ncount
  Pe <- ((b/ncount) * (a/ncount)) + ((d/ncount) * (c/ncount))
  Px1 <- Po - Pe
  Px2 <- 1 - Pe
  Px3 <- Px1/Px2
  
  tx1 <- xx + yyy
  tx2 <- 2 * a * c
  tx3 <- a - c
  tx4 <- xx - yyy
  tx5 <- ncount * (( 2 * a ) - tx4)
  tx6 <- ncount * tx1
  
  kappamax <- (tx6 - tx2 - (tx3 * tx4)) / ((tx5 - tx3) - (tx3 * tx4))
  
  return(tss) #returns the TSS value alone to an object file
  
  #cat(" Maxent results for model with\n",a,"training sample predictions\n",c ,"background predictions\n\n TSS value:        ", tss,"\n Overall accuracy: ",overallaccuracy,"\n Sensitivity:      ",sensitivity,"\n Specificity:      ",specificity,"\n Kappa:            ",kappa,"\n Kappa max:        ",kappamax)
  
}

### READ INPUTS ####

#read in predictor layers
#current Env
modelclimEnv <- brick("modelclimEnv.tif")
plot(modelclimEnv)
names(modelclimEnv) <- c('bio1', 'bio4', 'bio12', 'bio15', 'bio19', 'bulk_density', 'cec', 'clay', 'soilph', 'silt', 'sand')

soils.ph <- dropLayer(modelclimEnv, c('bio1', 'bio4', 'bio12', 'bio15', 'bio19'))
plot(soils.ph)

#future Env layers

#2061-2080
#CNRM-CM5 
#RCP 4.5
setwd("D:/DATA/GIS_Layers/CHELSA_Future/CNRM-CM5/2061-2080/RCP4.5")
Files1 <- list.files(pattern = ".tif$")
Files1
futureEnvCC.c <- stack(Files1)

futureEnvCC.c <- crop(futureEnvCC.c, modelclimEnv)
futureEnvCC.c <- mask(futureEnvCC.c, modelclimEnv$bio1) 

futureEnvCC.c <- stack(futureEnvCC.c, soils.ph)
names(futureEnvCC.c) <- names(modelclimEnv)

plot(futureEnvCC.c)
futureEnvCC.c


#RCP 8.5
setwd("D:/DATA/GIS_Layers/CHELSA_Future/CNRM-CM5/2061-2080/RCP8.5")
Files1 <- list.files(pattern = ".tif$")
Files1
futureEnvCC.d <- stack(Files1)

futureEnvCC.d <- crop(futureEnvCC.d, modelclimEnv)
futureEnvCC.d <- mask(futureEnvCC.d, modelclimEnv$bio1) 

futureEnvCC.d <- stack(futureEnvCC.d, soils.ph)
names(futureEnvCC.d) <- names(modelclimEnv)

plot(futureEnvCC.d)
futureEnvCC.d

#GFDL-CM3
#RCP 4.5
setwd("D:/DATA/GIS_Layers/CHELSA_Future/GFDL-CM3/2061-2080/RCP4.5")
Files1 <- list.files(pattern = ".tif$")
Files1
futureEnvGF.c <- stack(Files1)

futureEnvGF.c <- crop(futureEnvGF.c, modelclimEnv)
futureEnvGF.c <- mask(futureEnvGF.c, modelclimEnv$bio1) 

futureEnvGF.c <- stack(futureEnvGF.c, soils.ph)
names(futureEnvGF.c) <- names(modelclimEnv)

plot(futureEnvGF.c)
futureEnvGF.c

#RCP 8.5
setwd("D:/DATA/GIS_Layers/CHELSA_Future/GFDL-CM3/2061-2080/RCP8.5")
Files1 <- list.files(pattern = ".tif$")
Files1
futureEnvGF.d <- stack(Files1)

futureEnvGF.d <- crop(futureEnvGF.d, modelclimEnv)
futureEnvGF.d <- mask(futureEnvGF.d, modelclimEnv$bio1) 

futureEnvGF.d <- stack(futureEnvGF.d, soils.ph)
names(futureEnvGF.d) <- names(modelclimEnv)

plot(futureEnvGF.d)
futureEnvGF.d


#MPI-ESM-LR
#RCP 4.5
setwd("D:/DATA/GIS_Layers/CHELSA_Future/MPI-ESM-LR/2061-2080/RCP4.5")
Files1 <- list.files(pattern = ".tif$")
Files1
futureEnvMP.c <- stack(Files1)

futureEnvMP.c <- crop(futureEnvMP.c, modelclimEnv)
futureEnvMP.c <- mask(futureEnvMP.c, modelclimEnv$bio1) 

futureEnvMP.c <- stack(futureEnvMP.c, soils.ph)
names(futureEnvMP.c) <- names(modelclimEnv)

plot(futureEnvMP.c)
futureEnvMP.c

#RCP 8.5
setwd("D:/DATA/GIS_Layers/CHELSA_Future/MPI-ESM-LR/2061-2080/RCP8.5")
Files1 <- list.files(pattern = ".tif$")
Files1
futureEnvMP.d <- stack(Files1)

futureEnvMP.d <- crop(futureEnvMP.d, modelclimEnv)
futureEnvMP.d <- mask(futureEnvMP.d, modelclimEnv$bio1) 

futureEnvMP.d <- stack(futureEnvMP.d, soils.ph)
names(futureEnvMP.d) <- names(modelclimEnv)

plot(futureEnvMP.d)
futureEnvMP.d


#read in species points
setwd("D:/PROJECTS/Rafflesia_SDM")
datapoints <- read.csv("maxent_clean2_Tetrastigma_Rafflesia_30s.csv")
summary(datapoints)
datapoints$Species <- as.factor(datapoints$Species)
unique(datapoints$Species)

### TETRASTIGMA and RAFFLESIA A1 RUNS ####
sp.name <- 'T_harmandii_ss01'

predictors <- modelclimEnv
predictors <- dropLayer(modelclimEnv,c('bio4', 'bio15', 'bio19', 'bulk_density',  'cec', 'clay', 'soilph', 'silt',  'sand'))
predictors <- dropLayer(modelclimEnv, c('silt', 'cec'))
predictors #to check 

sp.occ <- filter(datapoints, Species == "Tetrastigma harmandii") #Tetrastigma ellipticum s.l., Tetrastigma loheri, Tetrastigma sp. A, Tetrastigma cf. magnum, Tetrastigma harmandii, Rafflesia lagascae, Rafflesia speciosa, Rafflesia lobata     
sp.occ <- as.data.frame(sp.occ)
head(sp.occ)
summary(sp.occ)
sp.occ <- sp.occ[,-1]
sp.occ

# ENMEval 
set.seed(123)
rm<-seq(from =1, to = 5, by = 0.5)
tune.args<- list(fc=c("L", "LQ", "LQH", "LQHP", "LQHPT", "H"),rm=rm)

#if occurrence points are more than 25 (see ENMEval vignette)
#enmeval_results <- ENMevaluate(sp.occ, predictors, method="block", n.bg=10000,  algorithm='maxent.jar',tune.args =tune.args)

#if occurrence points are less than or equal to 25
enmeval_results <- ENMevaluate(sp.occ, predictors, method="jackknife", n.bg=10000,  algorithm='maxent.jar',tune.args =tune.args)
#?ENMeval


#enmeval_results@results 
#find the model that has the lowest delta AICc

write.csv(enmeval_results@results, "Tharmandii_6reps_13sep2022ss1_enmeval_results_30s.csv", row.names = F)


# MAXENT run using the results from ENMEval
# NOTE: CHANGE MAXENT SETTINGS FOR EACH SPECIES ACCORDING TO ENMEVAL RESULTS!!!
getwd()
sp.model <- maxent(
  x=predictors,
  p=sp.occ,
  J = TRUE,
  path=paste0(sp.name),
  burnin=5000,
  nbg=10000,
  args=c(
    'betamultiplier=1',
    'noautofeature',
    'linear',
    'quadratic',
    'nohinge',
    'noproduct',
    'nothreshold',
    'writebackgroundpredictions=true',
    'doclamp',
    'threads=3',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=true',
    'removeduplicates=true',
    'replicatetype=bootstrap',
    'outputformat=Cloglog',
    'replicates=6',
    'maximumiterations=500000'
  )
)

# Compile AUC and TSS results

for (i in 0:5) { 
  
  #read in the background predictions
  backgroundpredictions <- read.csv(paste0(sp.name, "/species_", i,"_backgroundPredictions.csv"))
  
  #we need the last column so will set the number as x
  x <- length(backgroundpredictions)
  
  #extract the cloglog/logistic results
  backgroundclog <- backgroundpredictions[,x]
  
  #now read in the sample predictions for testing
  samplepredictions <- read.csv(paste0(sp.name,"/species_",i,"_samplePredictions.csv"))
  
  #we need the last column again of logistic or cloglog predictions so set a second x
  x2 <- length(samplepredictions)
  
  #extract the cloglog/logistic results for sample
  sampleclog <- samplepredictions[,x2]
  
  #set n the number of pseudoabsences used for background predictions by MaxEnt
  n <- 10000
  
  #read in maxent results
  maxres <- read.csv(paste0(sp.name,"/maxentResults.csv"))
  
  
  ### Set the threshold rule here
  th <- maxres[i+1,"Minimum.training.presence.Cloglog.threshold"]
  
  #run the function, the input values are the sampleclog values, then the background clog values, the sample number for the pseudo absences and then threshold value
  ts <- TSS_calculations(sampleclog,backgroundclog,n,th)
  
  if (i==0) {
    tss.df <- ts
  } else {
    tss.df <- rbind(tss.df, ts)
  }
  
}

tss.df <- as.data.frame(tss.df)
tss.df[7,1] <- summarise(tss.df, avg = mean(V1))

results.df <- cbind(maxres[,1], maxres[,6], tss.df)
colnames(results.df) <- c("Replicate", "Training.AUC", "TSS")
results.df

#setwd("D:/PROJECTS/Rafflesia_SDM")
write.csv(results.df, file = paste0(sp.name,"/",sp.name,"_Maxent_AUC_TSS_results.csv"), row.names = F)


getwd()
# predict to entire dataset
sp.pred.map <- predict(object=sp.model,
                       x=predictors,
                       filename=paste0(sp.name,"/",sp.name,"_R"),
                       na.rm=TRUE,
                       format='GTiff',
                       overwrite=TRUE,
                       doclamp= TRUE
)

sp.pred.map.ave <- calc(sp.pred.map,
                        fun = mean,
                        na.rm = T,
                        filename=paste0(sp.name,"/",sp.name,"_ave"),
                        format='GTiff',
                        overwrite=T)


#create binary maps
#extract the minimum training presence threshold
d <- as.data.frame(sp.model@results)
a <- d["Minimum.training.presence.Cloglog.threshold", "species (average)"]

#create matrix for reclassification
n <- c(0,a,0,a,1,1)
binmat <- matrix(n, ncol=3, byrow=T)
binmat

#reclassify current map
sp.pred.map.bin <- reclassify(x=sp.pred.map.ave, 
                              rcl=binmat, 
                              filename=paste0(sp.name,"/",sp.name,"_currentEnv_bin_ave"), 
                              format="GTiff",
                              overwrite=T)

#windows()
par(mfrow=c(1,2))
plot(sp.pred.map.ave, main="Predicted Suitability")
points(sp.occ$Longitude, sp.occ$Latitude, pch="+", cex=0.5)
plot(sp.pred.map.bin, main="Binary", legend=F)
#points(sp.occ$lon, sp.occ$lat, pch="+", cex=0.5)


### 
#predict into future climate

#2061-2080
#rcp45
future.mod4a <- predict(object=sp.model, 
                        x=futureEnvCC.c, 
                        filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_CC_45"),
                        na.rm=TRUE,
                        format='GTiff',
                        overwrite=TRUE,
                        doclamp= TRUE)
future.mod4a.ave <- calc(future.mod4a,
                         fun = mean,
                         na.rm = T,
                         filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_CC_45_ave"),
                         format='GTiff',
                         overwrite=T)

future.mod5a <- predict(object=sp.model, 
                        x=futureEnvGF.c,
                        filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_GF_45"),
                        na.rm=TRUE,
                        format='GTiff',
                        overwrite=TRUE,
                        doclamp= TRUE)
future.mod5a.ave <- calc(future.mod5a,
                         fun = mean,
                         na.rm = T,
                         filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_GF_45_ave"),
                         format='GTiff',
                         overwrite=T)

future.mod6a <- predict(object=sp.model, 
                        x=futureEnvMP.c,
                        filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_MP_45"),
                        na.rm=TRUE,
                        format='GTiff',
                        overwrite=TRUE,
                        doclamp= TRUE)
future.mod6a.ave <- calc(future.mod6a,
                         fun = mean,
                         na.rm = T,
                         filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_MP_45_ave"),
                         format='GTiff',
                         overwrite=T)

#rcp85
future.mod4b <- predict(object=sp.model, 
                        x=futureEnvCC.d,
                        filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_CC_85"),
                        na.rm=TRUE,
                        format='GTiff',
                        overwrite=TRUE,
                        doclamp= TRUE)
future.mod4b.ave <- calc(future.mod4b,
                         fun = mean,
                         na.rm = T,
                         filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_CC_85_ave"),
                         format='GTiff',
                         overwrite=T)

future.mod5b <- predict(object=sp.model, 
                        x=futureEnvGF.d,
                        filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_GF_85"),
                        na.rm=TRUE,
                        format='GTiff',
                        overwrite=TRUE,
                        doclamp= TRUE)
future.mod5b.ave <- calc(future.mod5b,
                         fun = mean,
                         na.rm = T,
                         filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_GF_85_ave"),
                         format='GTiff',
                         overwrite=T)

future.mod6b <- predict(object=sp.model, 
                        x=futureEnvMP.d,
                        filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_MP_85"),
                        na.rm=TRUE,
                        format='GTiff',
                        overwrite=TRUE,
                        doclamp= TRUE)
future.mod6b.ave <- calc(future.mod6b,
                         fun = mean,
                         na.rm = T,
                         filename=paste0(sp.name,"/",sp.name,"_MaxEnt_2070_MP_85_ave"),
                         format='GTiff',
                         overwrite=T)



#average model predictions
#2070
fut.2070.45.ave <- (future.mod4a.ave + future.mod5a.ave + future.mod6a.ave) / 3
fut.2070.85.ave <- (future.mod4b.ave + future.mod5b.ave + future.mod6b.ave) / 3


#write it out
#future
writeRaster(fut.2070.45.ave, filename=paste0(sp.name,"/",sp.name,"_2070_rcp45_ensemble"), format='GTiff', overwrite=TRUE)
writeRaster(fut.2070.85.ave, filename=paste0(sp.name,"/",sp.name,"_2070_rcp85_ensemble"), format='GTiff', overwrite=TRUE)


#convert to binary
#averages
#2070
rcp45.70.ave.bin <- reclassify(x=fut.2070.45.ave, 
                               rcl=binmat, 
                               filename=paste0(sp.name,"/",sp.name,"_2070_RCP45_ensem_bin"),                               format="GTiff",
                               overwrite=T)
rcp85.70.ave.bin <- reclassify(x=fut.2070.85.ave, 
                               rcl=binmat, 
                               filename=paste0(sp.name,"/",sp.name,"_2070_RCP85_ensem_bin"),                                 format="GTiff",
                               overwrite=T)

#single models
#2070
#RCP4.5
future.mod4a.bin <- reclassify(x=future.mod4a.ave, 
                               rcl=binmat, 
                               filename=paste0(sp.name,"/",sp.name,"_2070_CC_RCP45_ave_bin"),                                 format="GTiff",
                               overwrite=T)

future.mod5a.bin <- reclassify(x=future.mod5a.ave, 
                               rcl=binmat, 
                               filename=paste0(sp.name,"/",sp.name,"_2070_GF_RCP45_ave_bin"),                                 format="GTiff",
                               overwrite=T)

future.mod6a.bin <- reclassify(x=future.mod6a.ave, 
                               rcl=binmat, 
                               filename=paste0(sp.name,"/",sp.name,"_2070_MP_RCP45_ave_bin"),                                 format="GTiff",
                               overwrite=T)

#RCP 8.5
future.mod4b.bin <- reclassify(x=future.mod4b.ave, 
                               rcl=binmat, 
                               filename=paste0(sp.name,"/",sp.name,"_2070_CC_RCP85_ave_bin"),                                 format="GTiff",
                               overwrite=T)

future.mod5b.bin <- reclassify(x=future.mod5b.ave, 
                               rcl=binmat, 
                               filename=paste0(sp.name,"/",sp.name,"_2070_GF_RCP85_ave_bin"),                                 format="GTiff",
                               overwrite=T)

future.mod6b.bin <- reclassify(x=future.mod6b.ave, 
                               rcl=binmat, 
                               filename=paste0(sp.name,"/",sp.name,"_2070_MP_RCP85_ave_bin"),                                 format="GTiff",
                               overwrite=T)


####COUNT NO OF PIXELS CHANGE ####

current.df <- as.data.frame(sp.pred.map.bin, na.rm = T)
colnames(current.df) <- "x"
current.sum <- count(current.df, x)

rcp45.70.df <- as.data.frame(rcp45.70.ave.bin, na.rm = T)
colnames(rcp45.70.df) <- "x"
rcp45.70.sum <- count(rcp45.70.df, x)

rcp85.70.df <- as.data.frame(rcp85.70.ave.bin, na.rm = T)
colnames(rcp85.70.df) <- "x"
rcp85.70.sum <- count(rcp85.70.df, x)

change.summ <- as.data.frame(cbind(current.sum$x, current.sum$n, rcp45.70.sum$n, rcp85.70.sum$n))
colnames(change.summ) <- c("State", "Current", "RCP4.5-2070", "RCP8.5-2070")

change.summ[3,1] <- "%Suitable Area"
change.summ[3,2] <- change.summ[2,2]/(change.summ[1,2] + change.summ[2,2])*100
change.summ[3,3] <- change.summ[2,3]/(change.summ[1,3] + change.summ[2,3])*100
change.summ[3,4] <- change.summ[2,4]/(change.summ[1,4] + change.summ[2,4])*100
change.summ

write.csv(change.summ, file= paste0(sp.name,"/",sp.name,"_pixel_changes.csv"), row.names = F)


####PLOTS #### 
sp.name <- "Tetrastigma harmandii"

#Current suitability
#windows()
par(mfrow=c(1,2))
plot(sp.pred.map.ave, main = paste0(sp.name, " Predicted Suitability"))
points(sp.occ$Longitude, sp.occ$Latitude, pch="+", cex=0.5)
plot(sp.pred.map.bin, main = paste0(sp.name, " Binary Suitability"))

#### SUMMARY MAPS ####
#windows()
#continuous
par(mfrow=c(1,3))
plot(sp.pred.map.ave, main = paste0(sp.name, " current"))
plot(fut.2070.45.ave, main = paste0(sp.name, " 2070 RCP 4.5"))
plot(fut.2070.85.ave, main = paste0(sp.name, " 2070 RCP 8.5"))


#binary
par(mfrow=c(1,3))
plot(sp.pred.map.bin, main = paste0(sp.name, " current"))
plot(rcp45.70.ave.bin, main = paste0(sp.name, " 2070 RCP 4.5"))
plot(rcp85.70.ave.bin, main = paste0(sp.name, " 2070 RCP 8.5"))

#FUTURE suitabilities
#CONTINUOUS
#individual GCMs
#RCP 2.6
par(mfrow=c(2,4))
plot(future.mod4a.ave, main = paste0(sp.name, " CNRM-CM5 2070 RCP 4.5"))
plot(future.mod5a.ave, main = paste0(sp.name, " GFDL-CM3 2070 RCP 4.5"))
plot(future.mod6a.ave, main = paste0(sp.name, " MPI-ESM-LR 2070 RCP 4.5"))
plot(fut.2070.45.ave, main = paste0(sp.name, " ensemble 2070 RCP 4.5"))

#RCP 8.5
#par(mfrow=c(2,4))
plot(future.mod4b.ave, main = paste0(sp.name, " CNRM-CM5 2070 RCP 8.5"))
plot(future.mod5b.ave, main = paste0(sp.name, " GFDL-CM3 2070 RCP 8.5"))
plot(future.mod6b.ave, main = paste0(sp.name, " MPI-ESM-LR 2070 RCP 8.5"))
plot(fut.2070.85.ave, main = paste0(sp.name, " ensemble 2070 RCP 8.5"))

#BINARY MAPS

#individual GCMs

#rcp4.5
par(mfrow=c(2,3))
#2070
plot(future.mod4a.bin, main = paste0(sp.name, " CNRM-CM5 RCP 4.5 2070"))
plot(future.mod5a.bin, main = paste0(sp.name, " GFDL-CM3 RCP 4.5 2070"))
plot(future.mod6a.bin, main = paste0(sp.name, " MPI-ESM-LR RCP 4.5 2070"))
#rcp8.5
plot(future.mod4b.bin, main = paste0(sp.name, " CNRM-CM5 RCP 8.5 2070"))
plot(future.mod5b.bin, main = paste0(sp.name, " GFDL-CM3 RCP 8.5 2070"))
plot(future.mod6b.bin, main = paste0(sp.name, " MPI-ESM-LR RCP 8.5 2070"))


### RAFFLESIA A3 RUNS ####

#### PREDICTORS ####

# Read in Tetrastigma continuous SDM results 

tetloh <- raster("D:/PROJECTS/Rafflesia_SDM/T_loheri_01/T_loheri_01_ave.tif")
tetmag <- raster("D:/PROJECTS/Rafflesia_SDM/T_magnum_01/T_magnum_01_ave.tif")
tetspa <- raster("D:/PROJECTS/Rafflesia_SDM/T_spa_01ss/T_spa_01ss_ave.tif")
tethar <- raster("D:/PROJECTS/Rafflesia_SDM/T_harmandii_ss01/T_harmandii_ss01_ave.tif")
tetell <- raster("D:/PROJECTS/Rafflesia_SDM/T_ellipticum_01/T_ellipticum_01_ave.tif")

# Stack with the 
sp.name <- 'R_lobata_03A'

predictors <- modelclimEnv

#add Tetrastigma host SDM of each Rafflesia species
#predictors <- stack(predictors, tetmag, tethar) #speciosa
predictors <- stack(predictors, tetspa, tetloh, tetell) #lobata
#predictors <- stack(predictors, tetell, tetloh, tetspa) #lagascae


#remove layers with 0 contribution after initial run
#predictors <- dropLayer(predictors, c('cec', 'bulk_density', 'silt', 'bio12', 'bio19', 'T_harmandii_ss01_ave')) #speciosa
#predictors <- dropLayer(modelclimEnv, c('clay', 'bulk_density')) #lagascae
predictors <- dropLayer(predictors,c('bio12', 'bio1', 'silt', 'cec', 'clay', 'soilph', 'sand')) #lobata 

predictors #to check 
plot(predictors)

#filter data points for each species
sp.occ <- filter(datapoints, Species == "Rafflesia lobata") # Rafflesia lagascae, Rafflesia speciosa, Rafflesia lobata     
sp.occ <- as.data.frame(sp.occ)
head(sp.occ)
summary(sp.occ)
sp.occ <- sp.occ[,-1]
sp.occ


# ENMEval 

rm<-seq(from =1, to = 5, by = 0.5)
tune.args<- list(fc=c("L", "LQ", "LQH", "LQHP", "LQHPT", "H"),rm=rm)

#if occurrence points are more than 25 (see ENMEval vignette)
#enmeval_results <- ENMevaluate(sp.occ, predictors, method="block", n.bg=10000,  algorithm='maxent.jar',tune.args =tune.args)

#if occurrence points are less than or equal to 25
enmeval_results <- ENMevaluate(sp.occ, predictors, method="jackknife", n.bg=10000,  algorithm='maxent.jar',tune.args =tune.args)
#?ENMeval


#enmeval_results@results 
#find the model that has the lowest delta AICc

write.csv(enmeval_results@results, "Rlobata03ssB_w_Tetras_6reps_ENMeval_results_30s.csv", row.names = F)




# MAXENT run using the results from ENMEval
# NOTE: CHANGE MAXENT SETTINGS FOR EACH SPECIES ACCORDING TO ENMEVAL RESULTS!!!
getwd()
sp.model <- maxent(
  x=predictors,
  p=sp.occ,
  J = TRUE,
  path=paste0(sp.name),
  burnin=5000,
  nbg=10000,
  args=c(
    'betamultiplier=3',
    'noautofeature',
    'nolinear',
    'noquadratic',
    'hinge',
    'noproduct',
    'nothreshold',
    'writebackgroundpredictions=true',
    'doclamp',
    'threads=3',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=true',
    'removeduplicates=true',
    'replicatetype=bootstrap',
    'outputformat=Cloglog',
    'replicates=6',
    'maximumiterations=500000'
  )
)

# Compile AUC and TSS results

for (i in 0:5) { 
  
  #read in the background predictions
  backgroundpredictions <- read.csv(paste0(sp.name, "/species_", i,"_backgroundPredictions.csv"))
  
  #we need the last column so will set the number as x
  x <- length(backgroundpredictions)
  
  #extract the cloglog/logistic results
  backgroundclog <- backgroundpredictions[,x]
  
  #now read in the sample predictions for testing
  samplepredictions <- read.csv(paste0(sp.name,"/species_",i,"_samplePredictions.csv"))
  
  #we need the last column again of logistic or cloglog predictions so set a second x
  x2 <- length(samplepredictions)
  
  #extract the cloglog/logistic results for sample
  sampleclog <- samplepredictions[,x2]
  
  #set n the number of pseudoabsences used for background predictions by MaxEnt
  n <- 10000
  
  #read in maxent results
  maxres <- read.csv(paste0(sp.name,"/maxentResults.csv"))
  
  
  ### Set the threshold rule here
  th <- maxres[i+1,"Minimum.training.presence.Cloglog.threshold"]
  
  #run the function, the input values are the sampleclog values, then the background clog values, the sample number for the pseudo absences and then threshold value
  ts <- TSS_calculations(sampleclog,backgroundclog,n,th)
  
  if (i==0) {
    tss.df <- ts
  } else {
    tss.df <- rbind(tss.df, ts)
  }
  
}

tss.df <- as.data.frame(tss.df)
tss.df[7,1] <- summarise(tss.df, avg = mean(V1))

results.df <- cbind(maxres[,1], maxres[,6], tss.df)
colnames(results.df) <- c("Replicate", "Training.AUC", "TSS")
results.df

#setwd("D:/PROJECTS/Rafflesia_SDM")
write.csv(results.df, file = paste0(sp.name,"/",sp.name,"_Maxent_AUC_TSS_results.csv"), row.names = F)


getwd()
# predict to entire dataset
sp.pred.map <- predict(object=sp.model,
                       x=predictors,
                       filename=paste0(sp.name,"/",sp.name,"_R"),
                       na.rm=TRUE,
                       format='GTiff',
                       overwrite=TRUE,
                       doclamp= TRUE
)

sp.pred.map.ave <- calc(sp.pred.map,
                        fun = mean,
                        na.rm = T,
                        filename=paste0(sp.name,"/",sp.name,"_ave"),
                        format='GTiff',
                        overwrite=T)


#create binary maps
#extract the minimum training presence threshold
d <- as.data.frame(sp.model@results)
a <- d["Minimum.training.presence.Cloglog.threshold", "species (average)"]

#create matrix for reclassification
n <- c(0,a,0,a,1,1)
binmat <- matrix(n, ncol=3, byrow=T)
binmat

#reclassify current map
sp.pred.map.bin <- reclassify(x=sp.pred.map.ave, 
                              rcl=binmat, 
                              filename=paste0(sp.name,"/",sp.name,"_currentEnv_bin_ave"), 
                              format="GTiff",
                              overwrite=T)

#windows()
sp.name <- "Rafflesia lobata"
par(mfrow=c(1,2))
plot(sp.pred.map.ave, main=paste0(sp.name, " Predicted Suitability"))
points(sp.occ$Longitude, sp.occ$Latitude, pch="+", cex=0.5)
plot(sp.pred.map.bin, main=paste0(sp.name, " Binary"), legend=F)
#points(sp.occ$lon, sp.occ$lat, pch="+", cex=0.5)


### RAFFLESIA A2 ####

#### Current #### 
#read in binary rasters
#tetrastigma binary rasters
tetloh.b <- raster("D:/PROJECTS/Rafflesia_SDM/T_loheri_01/T_loheri_01_currentEnv_bin_ave.tif")
tethar.b <- raster("D:/PROJECTS/Rafflesia_SDM/T_harmandii_ss01/T_harmandii_ss01_currentEnv_bin_ave.tif")
tetell.b <- raster("D:/PROJECTS/Rafflesia_SDM/T_ellipticum_01/T_ellipticum_01_currentEnv_bin_ave.tif")
tetmag.b <- raster("D:/PROJECTS/Rafflesia_SDM/T_magnum_01/T_magnum_01_currentEnv_bin_ave.tif") 
tetspa.b <- raster("D:/PROJECTS/Rafflesia_SDM/T_spa_01ss/T_spa_01ss_currentEnv_bin_ave.tif")

#Rafflesias
rspe.b <- raster("D:/PROJECTS/Rafflesia_SDM/R_speciosa_01/R_speciosa_01_currentEnv_bin_ave.tif")
rlob.b <- raster("D:/PROJECTS/Rafflesia_SDM/R_lobata_01/R_lobata_01_currentEnv_bin_ave.tif")
rlag.b <- raster("D:/PROJECTS/Rafflesia_SDM/R_lagascae_01/R_lagascae_01_currentEnv_bin_ave.tif")

# approach 2: raster addition of rafflesia + tetrastigma hosts binary rasters
#predictors <- stack(predictors, tetmag, tethar) #speciosa
#predictors <- stack(predictors, tetspa, tetloh, tetell) #lobata
#predictors <- stack(predictors, tetell, tetloh, tetspa) #lagascae

# R. speciosa (tetmag + tethar)
rspe.a2 <- (rspe.b*10) + tetmag.b + tethar.b

rspe.a21 <- rspe.a2*0
rspe.a21[rspe.a2 >= 11] <- 3
rspe.a21[rspe.a2 <= 2 & rspe.a2 > 0] <- 2

rspe.a22 <- rspe.a2*0
rspe.a22[rspe.a2>=11] <- 1


writeRaster(rspe.a2, filename = "R_speciosa_01/R_speciosa_A2_currentEnv_w_hosts.tif", datatype='INT2S')
writeRaster(rspe.a21, filename = "R_speciosa_01/R_speciosa_A2_currentEnv_w_mergedhosts.tif", datatype='INT2S')
writeRaster(rspe.a22, filename = "R_speciosa_01/R_speciosa_A2_currentEnv.tif", datatype='INT2S')


#R. lobata (tetspa + tetloh+ tetell)
rlob.a2 <- (rlob.b*10) + tetspa.b + tetloh.b + tetell.b

rlob.a21 <- rlob.a2*0
rlob.a21[rlob.a2 >= 11] <- 3
rlob.a21[rlob.a2 <= 3 & rlob.a2 > 0] <- 2

rlob.a22 <- rlob.a2*0
rlob.a22[rlob.a2>=11] <- 1

writeRaster(rlob.a2, filename = "R_lobata_01/R_lobata_A2_currentEnv_w_hosts.tif", datatype='INT2S')
writeRaster(rlob.a21, filename = "R_lobata_01/R_lobata_A2_currentEnv_w_mergedhosts.tif", datatype='INT2S')
writeRaster(rlob.a22, filename = "R_lobata_01/R_lobata_A2_currentEnv.tif", datatype='INT2S')


#R. lagascae (tetell + tetloh + tetspa)
rlag.a2 <- (rlag.b*10) + tetell.b + tetloh.b + tetspa.b

rlag.a21 <- rlag.a2*0
rlag.a21[rlag.a2 >= 11] <- 3
rlag.a21[rlag.a2 <= 3 & rlag.a2 > 0] <- 2

rlag.a22 <- rlag.a2*0
rlag.a22[rlag.a2>=11] <- 1

writeRaster(rlag.a2, filename = "R_lagascae_01/R_lagascae_A2_currentEnv_w_hosts.tif", datatype='INT2S')
writeRaster(rlag.a21, filename = "R_lagascae_01/R_lagascae_A2_currentEnv_w_mergedhosts.tif", datatype='INT2S')
writeRaster(rlag.a22, filename = "R_lagascae_01/R_lagascae_A2_currentEnv.tif", datatype='INT2S')


#### RCP 4.5 ####
#read in binary rasters
#tetrastigma binary rasters
tetloh.b45 <- raster("D:/PROJECTS/Rafflesia_SDM/T_loheri_01/T_loheri_01_2070_RCP45_ensem_bin.tif")
tethar.b45 <- raster("D:/PROJECTS/Rafflesia_SDM/T_harmandii_ss01/T_harmandii_ss01_2070_RCP45_ensem_bin.tif")
tetell.b45 <- raster("D:/PROJECTS/Rafflesia_SDM/T_ellipticum_01/T_ellipticum_01_2070_RCP45_ensem_bin.tif")
tetmag.b45 <- raster("D:/PROJECTS/Rafflesia_SDM/T_magnum_01/T_magnum_01_2070_RCP45_ensem_bin.tif") 
tetspa.b45 <- raster("D:/PROJECTS/Rafflesia_SDM/T_spa_01ss/T_spa_01ss_2070_RCP45_ensem_bin.tif")

#Rafflesias
rspe.b45 <- raster("D:/PROJECTS/Rafflesia_SDM/R_speciosa_01/R_speciosa_01_2070_RCP45_ensem_bin.tif")
rlob.b45 <- raster("D:/PROJECTS/Rafflesia_SDM/R_lobata_01/R_lobata_01_2070_RCP45_ensem_bin.tif")
rlag.b45 <- raster("D:/PROJECTS/Rafflesia_SDM/R_lagascae_01/R_lagascae_01_2070_RCP45_ensem_bin.tif")


# R. speciosa (tetmag + tethar)
rspe.a2.45 <- (rspe.b45*10) + tetmag.b45 + tethar.b45

rspe.a21.45 <- rspe.a2.45*0
rspe.a21.45[rspe.a2.45 >= 11] <- 3
rspe.a21.45[rspe.a2.45 <= 2 & rspe.a2.45 > 0] <- 2

rspe.a22.45 <- rspe.a2.45*0
rspe.a22.45[rspe.a2.45>=11] <- 1


writeRaster(rspe.a2.45, filename = "R_speciosa_01/R_speciosa_A2_RCP45_w_hosts.tif", datatype='INT2S')
writeRaster(rspe.a21.45, filename = "R_speciosa_01/R_speciosa_A2_RCP45_w_mergedhosts.tif", datatype='INT2S')
writeRaster(rspe.a22.45, filename = "R_speciosa_01/R_speciosa_A2_RCP45.tif", datatype='INT2S')


#R. lobata (tetspa + tetloh+ tetell)
rlob.a2.45 <- (rlob.b45*10) + tetspa.b45 + tetloh.b45 + tetell.b45

rlob.a21.45 <- rlob.a2.45*0
rlob.a21.45[rlob.a2.45 >= 11] <- 3
rlob.a21.45[rlob.a2.45 <= 3 & rlob.a2.45 > 0] <- 2

rlob.a22.45 <- rlob.a2.45*0
rlob.a22.45[rlob.a2.45>=11] <- 1

writeRaster(rlob.a2.45, filename = "R_lobata_01/R_lobata_A2_RCP45_w_hosts.tif", datatype='INT2S')
writeRaster(rlob.a21.45, filename = "R_lobata_01/R_lobata_A2_RCP45_w_mergedhosts.tif", datatype='INT2S')
writeRaster(rlob.a22.45, filename = "R_lobata_01/R_lobata_A2_RCP45.tif", datatype='INT2S')


#R. lagascae (tetell + tetloh + tetspa)
rlag.a2.45 <- (rlag.b45*10) + tetell.b45 + tetloh.b45 + tetspa.b45

rlag.a21.45 <- rlag.a2.45*0
rlag.a21.45[rlag.a2.45 >= 11] <- 3
rlag.a21.45[rlag.a2.45 <= 3 & rlag.a2.45 > 0] <- 2

rlag.a22.45 <- rlag.a2.45*0
rlag.a22.45[rlag.a2.45>=11] <- 1

writeRaster(rlag.a2.45, filename = "R_lagascae_01/R_lagascae_A2_RCP45_w_hosts.tif", datatype='INT2S')
writeRaster(rlag.a21.45, filename = "R_lagascae_01/R_lagascae_A2_RCP45_w_mergedhosts.tif", datatype='INT2S')
writeRaster(rlag.a22.45, filename = "R_lagascae_01/R_lagascae_A2_RCP45.tif", datatype='INT2S')


#### RCP 8.5 ####
#read in binary rasters
#tetrastigma binary rasters
tetloh.b85 <- raster("D:/PROJECTS/Rafflesia_SDM/T_loheri_01/T_loheri_01_2070_RCP85_ensem_bin.tif")
tethar.b85 <- raster("D:/PROJECTS/Rafflesia_SDM/T_harmandii_ss01/T_harmandii_ss01_2070_RCP85_ensem_bin.tif")
tetell.b85 <- raster("D:/PROJECTS/Rafflesia_SDM/T_ellipticum_01/T_ellipticum_01_2070_RCP85_ensem_bin.tif")
tetmag.b85 <- raster("D:/PROJECTS/Rafflesia_SDM/T_magnum_01/T_magnum_01_2070_RCP85_ensem_bin.tif") 
tetspa.b85 <- raster("D:/PROJECTS/Rafflesia_SDM/T_spa_01ss/T_spa_01ss_2070_RCP85_ensem_bin.tif")

#Rafflesias
rspe.b85 <- raster("D:/PROJECTS/Rafflesia_SDM/R_speciosa_01/R_speciosa_01_2070_RCP85_ensem_bin.tif")
rlob.b85 <- raster("D:/PROJECTS/Rafflesia_SDM/R_lobata_01/R_lobata_01_2070_RCP85_ensem_bin.tif")
rlag.b85 <- raster("D:/PROJECTS/Rafflesia_SDM/R_lagascae_01/R_lagascae_01_2070_RCP85_ensem_bin.tif")


# R. speciosa (tetmag + tethar)
rspe.a2.85 <- (rspe.b85*10) + tetmag.b85 + tethar.b85

rspe.a21.85 <- rspe.a2.85*0
rspe.a21.85[rspe.a2.85 >= 11] <- 3
rspe.a21.85[rspe.a2.85 <= 2 & rspe.a2.85 > 0] <- 2

rspe.a22.85 <- rspe.a2.85*0
rspe.a22.85[rspe.a2.85>=11] <- 1


writeRaster(rspe.a2.85, filename = "R_speciosa_01/R_speciosa_A2_RCP85_w_hosts.tif", datatype='INT2S')
writeRaster(rspe.a21.85, filename = "R_speciosa_01/R_speciosa_A2_RCP85_w_mergedhosts.tif", datatype='INT2S')
writeRaster(rspe.a22.85, filename = "R_speciosa_01/R_speciosa_A2_RCP85.tif", datatype='INT2S')


#R. lobata (tetspa + tetloh+ tetell)
rlob.a2.85 <- (rlob.b85*10) + tetspa.b85 + tetloh.b85 + tetell.b85

rlob.a21.85 <- rlob.a2.85*0
rlob.a21.85[rlob.a2.85 >= 11] <- 3
rlob.a21.85[rlob.a2.85 <= 3 & rlob.a2.85 > 0] <- 2

rlob.a22.85 <- rlob.a2.85*0
rlob.a22.85[rlob.a2.85>=11] <- 1

writeRaster(rlob.a2.85, filename = "R_lobata_01/R_lobata_A2_RCP85_w_hosts.tif", datatype='INT2S')
writeRaster(rlob.a21.85, filename = "R_lobata_01/R_lobata_A2_RCP85_w_mergedhosts.tif", datatype='INT2S')
writeRaster(rlob.a22.85, filename = "R_lobata_01/R_lobata_A2_RCP85.tif", datatype='INT2S')


#R. lagascae (tetell + tetloh + tetspa)
rlag.a2.85 <- (rlag.b85*10) + tetell.b85 + tetloh.b85 + tetspa.b85

rlag.a21.85 <- rlag.a2.85*0
rlag.a21.85[rlag.a2.85 >= 11] <- 3
rlag.a21.85[rlag.a2.85 <= 3 & rlag.a2.85 > 0] <- 2

rlag.a22.85 <- rlag.a2.85*0
rlag.a22.85[rlag.a2.85>=11] <- 1

writeRaster(rlag.a2.85, filename = "R_lagascae_01/R_lagascae_A2_RCP85_w_hosts.tif", datatype='INT2S')
writeRaster(rlag.a21.85, filename = "R_lagascae_01/R_lagascae_A2_RCP85_w_mergedhosts.tif", datatype='INT2S')
writeRaster(rlag.a22.85, filename = "R_lagascae_01/R_lagascae_A2_RCP85.tif", datatype='INT2S')


### TETRASTIGMA COMBINATIONS ####

tet1 <- tetloh.b + tetspa.b + tetell.b
tet2 <- tetmag.b + tethar.b

par(mfrow=c(1,2))
plot(tet1, main = "T.loheri + T.sp.A + T.ellipticum")
plot(tet2, main = "T.magnum + T.harmandii")

writeRaster(tet1, filename = "Tetrastigma_SetA_LAE.tif", datatype='INT2S')
writeRaster(tet2, filename = "Tetrastigma_SetB_MH.tif", datatype='INT2S')


#### COMPARISON OF APPROACHES: USING CONTINUOUS MAXENT RESULTS, UNMASKED ####

#Rafflesia speciosa
rs1 <- enmtools.species(raster("R_speciosa_01/R_speciosa_01_currentEnv_bin_ave.tif"), species.name = "R.speciosa A1")
rs3 <- enmtools.species(raster("R_speciosa_03B/R_speciosa_03B_currentEnv_bin_ave.tif"), species.name = "R.speciosa A3")
rs2 <- enmtools.species(rspe.a22, species.name = "R.speciosa A2")

geog.range.overlap(rs1, rs3)
geog.range.overlap(rs1, rs2)
geog.range.overlap(rs2, rs3)

rs1c <- raster("R_speciosa_01/R_speciosa_01_ave.tif")
rs3c <- raster("R_speciosa_03B/R_speciosa_03B_ave.tif")

raster.overlap(rs1c, rs3c)

#Rafflesia lagascae
rlag1 <- enmtools.species(raster("R_lagascae_01/R_lagascae_01_currentEnv_bin_ave.tif"), species.name = "R.lagascae A1")
rlag3 <- enmtools.species(raster("R_lagascae_03A/R_lagascae_03A_currentEnv_bin_ave.tif"), species.name = "R.lagascae A3")
rlag2 <- enmtools.species(rlag.a22, species.name = "R.lagascae A2")

geog.range.overlap(rlag1, rlag3)
geog.range.overlap(rlag1, rlag2)
geog.range.overlap(rlag2, rlag3)

rlag1c <- raster("R_lagascae_01/R_lagascae_01_ave.tif")
rlag3c <- raster("R_lagascae_03A/R_lagascae_03A_ave.tif")

raster.overlap(rlag1c, rlag3c)

#Rafflesia lobata
rlob1 <- enmtools.species(raster("R_lobata_01/R_lobata_01_currentEnv_bin_ave.tif"), species.name = "R.lobata A1")
rlob3 <- enmtools.species(raster("R_lobata_03A/R_lobata_03A_currentEnv_bin_ave.tif"), species.name = "R.lobata A3")
rlob2 <- enmtools.species(rlob.a22, species.name = "R.lobata A2")

geog.range.overlap(rlob1, rlob3)
geog.range.overlap(rlob1, rlob2)
geog.range.overlap(rlob2, rlob3)

rlob1c <- raster("R_lobata_01/R_lobata_01_ave.tif")
rlob3c <- raster("R_lobata_03A/R_lobata_03A_ave.tif")

raster.overlap(rlob1c, rlob3c)

### FILTERING UNSUITABLE AREAS & COMPARISON OF APPROACHES, UNMASKED AND MASKED, BINARY RESULTS ####

unsuit <- raster("unsuit_mask.tif")
plot(unsuit)

#Convert to df
unsuit.df <- as.data.frame(unsuit, na.rm = T)
colnames(unsuit.df) <- "x"
unsuit.df.sum <- count(unsuit.df, x)


#### R. speciosa ####
sp.name <- "R_speciosa"
#A1
rspe1.cur <- raster("R_speciosa_01/R_speciosa_01_currentEnv_bin_ave.tif")
rspe1.45 <- raster("R_speciosa_01/R_speciosa_01_2070_RCP45_ensem_bin.tif")
rspe1.85 <- raster("R_speciosa_01/R_speciosa_01_2070_RCP85_ensem_bin.tif")

rspe1.cur.m <- (rspe1.cur*10) - unsuit
rspe1.cur.m[rspe1.cur.m < 10] <- 0
rspe1.cur.m[rspe1.cur.m == 10] <- 1

rspe1.45.m <- (rspe1.45*10) - unsuit
rspe1.45.m[rspe1.45.m < 10] <- 0
rspe1.45.m[rspe1.45.m == 10] <- 1

rspe1.85.m <- (rspe1.85*10) - unsuit
rspe1.85.m[rspe1.85.m < 10] <- 0
rspe1.85.m[rspe1.85.m == 10] <- 1

writeRaster(rspe1.cur.m, filename = "R_speciosa_01/R_speciosa_01_currentEnv_bin_ave_masked.tif", datatype='INT2S')
writeRaster(rspe1.45.m, filename = "R_speciosa_01/R_speciosa_01_2070_RCP45_ensem_bin_masked.tif", datatype='INT2S')
writeRaster(rspe1.85.m, filename = "R_speciosa_01/R_speciosa_01_2070_RCP85_ensem_bin_masked.tif", datatype='INT2S')

#A2
rspe2.cur <- raster("R_speciosa_01/R_speciosa_A2_currentEnv.tif")
rspe2.45.a2 <- raster("R_speciosa_01/R_speciosa_A2_RCP45.tif")
rspe2.85.a2 <- raster("R_speciosa_01/R_speciosa_A2_RCP85.tif")

rspe2.cur.m <- (rspe2.cur*10) - unsuit
rspe2.cur.m[rspe2.cur.m < 10] <- 0
rspe2.cur.m[rspe2.cur.m == 10] <- 1

rspe2.45.a2.m <- (rspe2.45.a2*10) - unsuit
rspe2.45.a2.m[rspe2.45.a2.m < 10] <- 0
rspe2.45.a2.m[rspe2.45.a2.m == 10] <- 1

rspe2.85.a2.m <- (rspe2.85.a2*10) - unsuit
rspe2.85.a2.m[rspe2.85.a2.m < 10] <- 0
rspe2.85.a2.m[rspe2.85.a2.m == 10] <- 1

writeRaster(rspe2.cur.m, filename = "R_speciosa_01/R_speciosa_A2_currentEnv_masked.tif", datatype='INT2S')
writeRaster(rspe2.45.a2.m, filename = "R_speciosa_01/R_speciosa_A2_2070_RCP45_ensem_bin_masked.tif", datatype='INT2S')
writeRaster(rspe2.85.a2.m, filename = "R_speciosa_01/R_speciosa_A2_2070_RCP85_ensem_bin_masked.tif", datatype='INT2S')

#A3
rspe3.cur <- raster("R_speciosa_03B/R_speciosa_03B_currentEnv_bin_ave.tif")

rspe3.cur.m <- (rspe3.cur*10) - unsuit
rspe3.cur.m[rspe3.cur.m < 10] <- 0
rspe3.cur.m[rspe3.cur.m == 10] <- 1

writeRaster(rspe3.cur.m, filename = "R_speciosa_03B/R_speciosa_03B_currentEnv_bin_ave_masked.tif", datatype='INT2S')


#Convert to df
# unmasked
currentA1.df <- as.data.frame(rspe1.cur, na.rm = T)
colnames(currentA1.df) <- "x"
currentA1.sum <- count(currentA1.df, x)

rcp45.70.df <- as.data.frame(rspe1.45, na.rm = T)
colnames(rcp45.70.df) <- "x"
rcp45.70.sum <- count(rcp45.70.df, x)

rcp85.70.df <- as.data.frame(rspe1.85, na.rm = T)
colnames(rcp85.70.df) <- "x"
rcp85.70.sum <- count(rcp85.70.df, x)

currentA2.df <- as.data.frame(rspe2.cur, na.rm = T)
colnames(currentA2.df) <- "x"
currentA2.sum <- count(currentA2.df, x)

rcp45.70.a2.df <- as.data.frame(rspe2.45.a2, na.rm = T)
colnames(rcp45.70.a2.df) <- "x"
rcp45.70.a2.sum <- count(rcp45.70.a2.df, x)

rcp85.70.a2.df <- as.data.frame(rspe2.85.a2, na.rm = T)
colnames(rcp85.70.a2.df) <- "x"
rcp85.70.a2.sum <- count(rcp85.70.a2.df, x)


currentA3.df <- as.data.frame(rspe3.cur, na.rm = T)
colnames(currentA3.df) <- "x"
currentA3.sum <- count(currentA3.df, x)


change.summ <- as.data.frame(cbind(currentA1.sum$x, currentA1.sum$n, currentA2.sum$n, currentA3.sum$n, rcp45.70.sum$n, rcp85.70.sum$n, rcp45.70.a2.sum$n, rcp85.70.a2.sum$n))
colnames(change.summ) <- c("State", "Current-A1", "Current-A2", "Current-A3", "RCP4.5-2070-A1", "RCP8.5-2070-A1", "RCP4.5-2070-A2", "RCP8.5-2070-A2")

change.summ[3,1] <- "%Suitable Area"
change.summ[3,2] <- change.summ[2,2]/(change.summ[1,2] + change.summ[2,2])*100
change.summ[3,3] <- change.summ[2,3]/(change.summ[1,3] + change.summ[2,3])*100
change.summ[3,4] <- change.summ[2,4]/(change.summ[1,4] + change.summ[2,4])*100
change.summ[3,5] <- change.summ[2,5]/(change.summ[1,5] + change.summ[2,5])*100
change.summ[3,6] <- change.summ[2,6]/(change.summ[1,6] + change.summ[2,6])*100
change.summ[3,7] <- change.summ[2,7]/(change.summ[1,7] + change.summ[2,7])*100
change.summ[3,8] <- change.summ[2,8]/(change.summ[1,8] + change.summ[2,8])*100
change.summ

write.csv(change.summ, file= paste0(sp.name,"_unmasked_pixel_changes_a1a2.csv"), row.names = F)


# masked
currentA1.df <- as.data.frame(rspe1.cur.m, na.rm = T)
colnames(currentA1.df) <- "x"
currentA1.sum <- count(currentA1.df, x)

rcp45.70.df <- as.data.frame(rspe1.45.m, na.rm = T)
colnames(rcp45.70.df) <- "x"
rcp45.70.sum <- count(rcp45.70.df, x)

rcp85.70.df <- as.data.frame(rspe1.85.m, na.rm = T)
colnames(rcp85.70.df) <- "x"
rcp85.70.sum <- count(rcp85.70.df, x)

currentA2.df <- as.data.frame(rspe2.cur.m, na.rm = T)
colnames(currentA2.df) <- "x"
currentA2.sum <- count(currentA2.df, x)

rcp45.a2.df <- as.data.frame(rspe2.45.a2.m, na.rm = T)
colnames(rcp45.a2.df) <- "x"
rcp45.a2.sum <- count(rcp45.a2.df, x)

rcp85.a2.df <- as.data.frame(rspe2.85.a2.m, na.rm = T)
colnames(rcp85.a2.df) <- "x"
rcp85.a2.sum <- count(rcp85.a2.df, x)

currentA3.df <- as.data.frame(rspe3.cur.m, na.rm = T)
colnames(currentA3.df) <- "x"
currentA3.sum <- count(currentA3.df, x)


change.summ.m <- as.data.frame(cbind(currentA1.sum$x, currentA1.sum$n, currentA2.sum$n, currentA3.sum$n, rcp45.70.sum$n, rcp85.70.sum$n, rcp45.a2.sum$n, rcp85.a2.sum$n))
colnames(change.summ.m) <- c("State", "Current-A1", "Current-A2", "Current-A3", "RCP4.5-2070-A1", "RCP8.5-2070-A1", "RCP4.5-2070-A2", "RCP8.5-2070-A2")

change.summ.m[3,1] <- "%Suitable Area"
change.summ.m[3,2] <- change.summ.m[2,2]/(change.summ.m[1,2] + change.summ.m[2,2])*100
change.summ.m[3,3] <- change.summ.m[2,3]/(change.summ.m[1,3] + change.summ.m[2,3])*100
change.summ.m[3,4] <- change.summ.m[2,4]/(change.summ.m[1,4] + change.summ.m[2,4])*100
change.summ.m[3,5] <- change.summ.m[2,5]/(change.summ.m[1,5] + change.summ.m[2,5])*100
change.summ.m[3,6] <- change.summ.m[2,6]/(change.summ.m[1,6] + change.summ.m[2,6])*100
change.summ.m[3,7] <- change.summ.m[2,7]/(change.summ.m[1,7] + change.summ.m[2,7])*100
change.summ.m[3,8] <- change.summ.m[2,8]/(change.summ.m[1,8] + change.summ.m[2,8])*100
change.summ.m

write.csv(change.summ.m, file= paste0(sp.name,"_masked_pixel_changes_a1a2.csv"), row.names = F)


#compare geographic range overlaps
#masked
# %overlapping areas
geog.range.overlap(enmtools.species(rspe1.cur.m), enmtools.species(rspe2.cur.m))
geog.range.overlap(enmtools.species(rspe1.cur.m), enmtools.species(rspe3.cur.m))
geog.range.overlap(enmtools.species(rspe2.cur.m), enmtools.species(rspe3.cur.m))

# indices
raster.overlap(rspe1.cur.m, rspe2.cur.m)
raster.overlap(rspe1.cur.m, rspe3.cur.m)
raster.overlap(rspe2.cur.m, rspe3.cur.m)

#unmasked
# %overlapping areas
geog.range.overlap(enmtools.species(rspe1.cur), enmtools.species(rspe2.cur))
geog.range.overlap(enmtools.species(rspe1.cur), enmtools.species(rspe3.cur))
geog.range.overlap(enmtools.species(rspe2.cur), enmtools.species(rspe3.cur))

# indices
raster.overlap(rspe1.cur, rspe2.cur)
raster.overlap(rspe1.cur, rspe3.cur)
raster.overlap(rspe2.cur, rspe3.cur)


#plot
#Approaches, current climate

par(mfrow=c(1,3))
plot(rspe1.cur.m, main = paste0(sp.name," A1"))
plot(rspe2.cur.m, main = paste0(sp.name," A2"))
plot(rspe3.cur.m, main = paste0(sp.name," A3"))

#A1 with future
plot(rspe1.cur.m, main = paste0(sp.name," CURRENT-A1"))
plot(rspe1.45.m, main = paste0(sp.name," RCP4.5-A1"))
plot(rspe1.85.m, main = paste0(sp.name," RCP8.5-A1"))

#A2 with future
plot(rspe2.cur.m, main = paste0(sp.name," CURRENT-A2"))
plot(rspe2.45.a2.m, main = paste0(sp.name," RCP4.5-A2"))
plot(rspe2.85.a2.m, main = paste0(sp.name," RCP8.5-A2"))


#### R. lagascae ####
sp.name <- "R_lagascae"
#A1
rlag1.cur <- raster("R_lagascae_01/R_lagascae_01_currentEnv_bin_ave.tif")
rlag1.45 <- raster("R_lagascae_01/R_lagascae_01_2070_RCP45_ensem_bin.tif")
rlag1.85 <- raster("R_lagascae_01/R_lagascae_01_2070_RCP85_ensem_bin.tif")

rlag1.cur.m <- (rlag1.cur*10) - unsuit
rlag1.cur.m[rlag1.cur.m < 10] <- 0
rlag1.cur.m[rlag1.cur.m == 10] <- 1

rlag1.45.m <- (rlag1.45*10) - unsuit
rlag1.45.m[rlag1.45.m < 10] <- 0
rlag1.45.m[rlag1.45.m == 10] <- 1

rlag1.85.m <- (rlag1.85*10) - unsuit
rlag1.85.m[rlag1.85.m < 10] <- 0
rlag1.85.m[rlag1.85.m == 10] <- 1

writeRaster(rlag1.cur.m, filename = "R_lagascae_01/R_lagascae_01_currentEnv_bin_ave_masked.tif", datatype='INT2S')
writeRaster(rlag1.45.m, filename = "R_lagascae_01/R_lagascae_01_2070_RCP45_ensem_bin_masked.tif", datatype='INT2S')
writeRaster(rlag1.85.m, filename = "R_lagascae_01/R_lagascae_01_2070_RCP85_ensem_bin_masked.tif", datatype='INT2S')

#A2
rlag2.cur <- raster("R_lagascae_01/R_lagascae_A2_currentEnv.tif")
rlag2.45.a2 <- raster("R_lagascae_01/R_lagascae_A2_RCP45.tif")
rlag2.85.a2 <- raster("R_lagascae_01/R_lagascae_A2_RCP85.tif")

rlag2.cur.m <- (rlag2.cur*10) - unsuit
rlag2.cur.m[rlag2.cur.m < 10] <- 0
rlag2.cur.m[rlag2.cur.m == 10] <- 1

rlag2.45.a2.m <- (rlag2.45.a2*10) - unsuit
rlag2.45.a2.m[rlag2.45.a2.m < 10] <- 0
rlag2.45.a2.m[rlag2.45.a2.m == 10] <- 1

rlag2.85.a2.m <- (rlag2.85.a2*10) - unsuit
rlag2.85.a2.m[rlag2.85.a2.m < 10] <- 0
rlag2.85.a2.m[rlag2.85.a2.m == 10] <- 1

writeRaster(rlag2.cur.m, filename = "R_lagascae_01/R_lagascae_A2_currentEnv_masked.tif", datatype='INT2S')
writeRaster(rlag2.45.a2.m, filename = "R_lagascae_01/R_lagascae_A2_2070_RCP45_ensem_bin_masked.tif", datatype='INT2S')
writeRaster(rlag2.85.a2.m, filename = "R_lagascae_01/R_lagascae_A2_2070_RCP85_ensem_bin_masked.tif", datatype='INT2S')

#A3
rlag3.cur <- raster("R_lagascae_03A/R_lagascae_03A_currentEnv_bin_ave.tif")

rlag3.cur.m <- (rlag3.cur*10) - unsuit
rlag3.cur.m[rlag3.cur.m < 10] <- 0
rlag3.cur.m[rlag3.cur.m == 10] <- 1

writeRaster(rlag3.cur.m, filename = "R_lagascae_03A/R_lagascae_03A_currentEnv_bin_ave_masked.tif", datatype='INT2S')

#Convert to df
#unmasked

currentA1.df <- as.data.frame(rlag1.cur, na.rm = T)
colnames(currentA1.df) <- "x"
currentA1.sum <- count(currentA1.df, x)

rcp45.70.df <- as.data.frame(rlag1.45, na.rm = T)
colnames(rcp45.70.df) <- "x"
rcp45.70.sum <- count(rcp45.70.df, x)

rcp85.70.df <- as.data.frame(rlag1.85, na.rm = T)
colnames(rcp85.70.df) <- "x"
rcp85.70.sum <- count(rcp85.70.df, x)

currentA2.df <- as.data.frame(rlag2.cur, na.rm = T)
colnames(currentA2.df) <- "x"
currentA2.sum <- count(currentA2.df, x)

rcp45.70.a2.df <- as.data.frame(rlag2.45.a2, na.rm = T)
colnames(rcp45.70.a2.df) <- "x"
rcp45.70.a2.sum <- count(rcp45.70.a2.df, x)

rcp85.70.a2.df <- as.data.frame(rlag2.85.a2, na.rm = T)
colnames(rcp85.70.a2.df) <- "x"
rcp85.70.a2.sum <- count(rcp85.70.a2.df, x)

currentA3.df <- as.data.frame(rlag3.cur, na.rm = T)
colnames(currentA3.df) <- "x"
currentA3.sum <- count(currentA3.df, x)


change.summ <- as.data.frame(cbind(currentA1.sum$x, currentA1.sum$n, currentA2.sum$n, currentA3.sum$n, rcp45.70.sum$n, rcp85.70.sum$n, rcp45.70.a2.sum$n, rcp85.70.a2.sum$n))
colnames(change.summ) <- c("State", "Current-A1", "Current-A2", "Current-A3", "RCP4.5-2070-A1", "RCP8.5-2070-A1", "RCP4.5-2070-A2", "RCP8.5-2070-A2")

change.summ[3,1] <- "%Suitable Area"
change.summ[3,2] <- change.summ[2,2]/(change.summ[1,2] + change.summ[2,2])*100
change.summ[3,3] <- change.summ[2,3]/(change.summ[1,3] + change.summ[2,3])*100
change.summ[3,4] <- change.summ[2,4]/(change.summ[1,4] + change.summ[2,4])*100
change.summ[3,5] <- change.summ[2,5]/(change.summ[1,5] + change.summ[2,5])*100
change.summ[3,6] <- change.summ[2,6]/(change.summ[1,6] + change.summ[2,6])*100
change.summ[3,7] <- change.summ[2,7]/(change.summ[1,7] + change.summ[2,7])*100
change.summ[3,8] <- change.summ[2,8]/(change.summ[1,8] + change.summ[2,8])*100
change.summ

write.csv(change.summ, file= paste0(sp.name,"_unmasked_pixel_changes_A1A2.csv"), row.names = F)


#masked
currentA1.df <- as.data.frame(rlag1.cur.m, na.rm = T)
colnames(currentA1.df) <- "x"
currentA1.sum <- count(currentA1.df, x)

rcp45.70.df <- as.data.frame(rlag1.45.m, na.rm = T)
colnames(rcp45.70.df) <- "x"
rcp45.70.sum <- count(rcp45.70.df, x)

rcp85.70.df <- as.data.frame(rlag1.85.m, na.rm = T)
colnames(rcp85.70.df) <- "x"
rcp85.70.sum <- count(rcp85.70.df, x)

currentA2.df <- as.data.frame(rlag2.cur.m, na.rm = T)
colnames(currentA2.df) <- "x"
currentA2.sum <- count(currentA2.df, x)

rcp45.70.a2.df <- as.data.frame(rlag2.45.a2.m, na.rm = T)
colnames(rcp45.70.a2.df) <- "x"
rcp45.70.a2.sum <- count(rcp45.70.a2.df, x)

rcp85.70.a2.df <- as.data.frame(rlag2.85.a2.m, na.rm = T)
colnames(rcp85.70.a2.df) <- "x"
rcp85.70.a2.sum <- count(rcp85.70.a2.df, x)

currentA3.df <- as.data.frame(rlag3.cur.m, na.rm = T)
colnames(currentA3.df) <- "x"
currentA3.sum <- count(currentA3.df, x)


change.summ.m <- as.data.frame(cbind(currentA1.sum$x, currentA1.sum$n, currentA2.sum$n, currentA3.sum$n, rcp45.70.sum$n, rcp85.70.sum$n, rcp45.70.a2.sum$n, rcp85.70.a2.sum$n))
colnames(change.summ.m) <- c("State", "Current-A1", "Current-A2", "Current-A3", "RCP4.5-2070-A1", "RCP8.5-2070-A1", "RCP4.5-2070-A2", "RCP8.5-2070-A2")

change.summ.m[3,1] <- "%Suitable Area"
change.summ.m[3,2] <- change.summ.m[2,2]/(change.summ.m[1,2] + change.summ.m[2,2])*100
change.summ.m[3,3] <- change.summ.m[2,3]/(change.summ.m[1,3] + change.summ.m[2,3])*100
change.summ.m[3,4] <- change.summ.m[2,4]/(change.summ.m[1,4] + change.summ.m[2,4])*100
change.summ.m[3,5] <- change.summ.m[2,5]/(change.summ.m[1,5] + change.summ.m[2,5])*100
change.summ.m[3,6] <- change.summ.m[2,6]/(change.summ.m[1,6] + change.summ.m[2,6])*100
change.summ.m[3,7] <- change.summ.m[2,7]/(change.summ.m[1,7] + change.summ.m[2,7])*100
change.summ.m[3,8] <- change.summ.m[2,8]/(change.summ.m[1,8] + change.summ.m[2,8])*100
change.summ.m

write.csv(change.summ.m, file= paste0(sp.name,"_masked_pixel_changes_A1A2.csv"), row.names = F)

#compare geographic range overlaps
#unmasked
# %overlapping areas
geog.range.overlap(enmtools.species(rlag1.cur), enmtools.species(rlag2.cur))
geog.range.overlap(enmtools.species(rlag1.cur), enmtools.species(rlag3.cur))
geog.range.overlap(enmtools.species(rlag2.cur), enmtools.species(rlag3.cur))

# indices
raster.overlap(rlag1.cur, rlag2.cur)
raster.overlap(rlag1.cur, rlag3.cur)
raster.overlap(rlag2.cur, rlag3.cur)


#masked
# %overlapping areas
geog.range.overlap(enmtools.species(rlag1.cur.m), enmtools.species(rlag2.cur.m))
geog.range.overlap(enmtools.species(rlag1.cur.m), enmtools.species(rlag3.cur.m))
geog.range.overlap(enmtools.species(rlag2.cur.m), enmtools.species(rlag3.cur.m))

# indices
raster.overlap(rlag1.cur.m, rlag2.cur.m)
raster.overlap(rlag1.cur.m, rlag3.cur.m)
raster.overlap(rlag2.cur.m, rlag3.cur.m)


#plot
#Approaches, current climate

par(mfrow=c(1,3))
plot(rlag1.cur.m, main = paste0(sp.name," A1"))
plot(rlag2.cur.m, main = paste0(sp.name," A2"))
plot(rlag3.cur.m, main = paste0(sp.name," A3"))

#A1 with future
plot(rlag1.cur.m, main = paste0(sp.name," CURRENT-A1"))
plot(rlag1.45.m, main = paste0(sp.name," RCP4.5-A1"))
plot(rlag1.85.m, main = paste0(sp.name," RCP8.5-A1"))


#A2 with future
plot(rlag2.cur.m, main = paste0(sp.name," CURRENT-A2"))
plot(rlag2.45.a2.m, main = paste0(sp.name," RCP4.5-A2"))
plot(rlag2.85.a2.m, main = paste0(sp.name," RCP8.5-A2"))

#### R. lobata ####
sp.name <- "R_lobata"
#A1
rlob1.cur <- raster("R_lobata_01/R_lobata_01_currentEnv_bin_ave.tif")
rlob1.45 <- raster("R_lobata_01/R_lobata_01_2070_RCP45_ensem_bin.tif")
rlob1.85 <- raster("R_lobata_01/R_lobata_01_2070_RCP85_ensem_bin.tif")

rlob1.cur.m <- (rlob1.cur*10) - unsuit
rlob1.cur.m[rlob1.cur.m < 10] <- 0
rlob1.cur.m[rlob1.cur.m == 10] <- 1

rlob1.45.m <- (rlob1.45*10) - unsuit
rlob1.45.m[rlob1.45.m < 10] <- 0
rlob1.45.m[rlob1.45.m == 10] <- 1

rlob1.85.m <- (rlob1.85*10) - unsuit
rlob1.85.m[rlob1.85.m < 10] <- 0
rlob1.85.m[rlob1.85.m == 10] <- 1

writeRaster(rlob1.cur.m, filename = "R_lobata_01/R_lobata_01_currentEnv_bin_ave_masked.tif", datatype='INT2S')
writeRaster(rlob1.45.m, filename = "R_lobata_01/R_lobata_01_2070_RCP45_ensem_bin_masked.tif", datatype='INT2S')
writeRaster(rlob1.85.m, filename = "R_lobata_01/R_lobata_01_2070_RCP85_ensem_bin_masked.tif", datatype='INT2S')

#A2
rlob2.cur <- raster("R_lobata_01/R_lobata_A2_currentEnv.tif")
rlob2.45 <- raster("R_lobata_01/R_lobata_A2_RCP45.tif")
rlob2.85 <- raster("R_lobata_01/R_lobata_A2_RCP85.tif")


rlob2.cur.m <- (rlob2.cur*10) - unsuit
rlob2.cur.m[rlob2.cur.m < 10] <- 0
rlob2.cur.m[rlob2.cur.m == 10] <- 1

rlob2.45.m <- (rlob2.45*10) - unsuit
rlob2.45.m[rlob2.45.m < 10] <- 0
rlob2.45.m[rlob2.45.m == 10] <- 1

rlob2.85.m <- (rlob2.85*10) - unsuit
rlob2.85.m[rlob2.85.m < 10] <- 0
rlob2.85.m[rlob2.85.m == 10] <- 1

writeRaster(rlob2.cur.m, filename = "R_lobata_01/R_lobata_A2_currentEnv_masked.tif", datatype='INT2S')
writeRaster(rlob2.45.m, filename = "R_lobata_01/R_lobata_A2_RCP45_masked.tif", datatype='INT2S')
writeRaster(rlob2.85.m, filename = "R_lobata_01/R_lobata_A2_RCP85_masked.tif", datatype='INT2S')

#A3
rlob3.cur <- raster("R_lobata_03A/R_lobata_03A_currentEnv_bin_ave.tif")

rlob3.cur.m <- (rlob3.cur*10) - unsuit
rlob3.cur.m[rlob3.cur.m < 10] <- 0
rlob3.cur.m[rlob3.cur.m == 10] <- 1

writeRaster(rlob3.cur.m, filename = "R_lobata_03A/R_lobata_03A_currentEnv_bin_ave_masked.tif", datatype='INT2S')

#Convert to df
#unmasked

currentA1.df <- as.data.frame(rlob1.cur, na.rm = T)
colnames(currentA1.df) <- "x"
currentA1.sum <- count(currentA1.df, x)

rcp45.70.df <- as.data.frame(rlob1.45, na.rm = T)
colnames(rcp45.70.df) <- "x"
rcp45.70.sum <- count(rcp45.70.df, x)

rcp85.70.df <- as.data.frame(rlob1.85, na.rm = T)
colnames(rcp85.70.df) <- "x"
rcp85.70.sum <- count(rcp85.70.df, x)

currentA2.df <- as.data.frame(rlob2.cur, na.rm = T)
colnames(currentA2.df) <- "x"
currentA2.sum <- count(currentA2.df, x)

rcp45.70.a2.df <- as.data.frame(rlob2.45, na.rm = T)
colnames(rcp45.70.a2.df) <- "x"
rcp45.70.a2.sum <- count(rcp45.70.a2.df, x)

rcp85.70.a2.df <- as.data.frame(rlob2.85, na.rm = T)
colnames(rcp85.70.a2.df) <- "x"
rcp85.70.a2.sum <- count(rcp85.70.a2.df, x)

currentA3.df <- as.data.frame(rlob3.cur, na.rm = T)
colnames(currentA3.df) <- "x"
currentA3.sum <- count(currentA3.df, x)


change.summ <- as.data.frame(cbind(currentA1.sum$x, currentA1.sum$n, currentA2.sum$n, currentA3.sum$n, rcp45.70.sum$n, rcp85.70.sum$n, rcp45.70.a2.sum$n, rcp85.70.a2.sum$n))
colnames(change.summ) <- c("State", "Current-A1", "Current-A2", "Current-A3", "RCP4.5-2070-A1", "RCP8.5-2070-A1",  "RCP4.5-2070-A2", "RCP8.5-2070-A2")

change.summ[3,1] <- "%Suitable Area"
change.summ[3,2] <- change.summ[2,2]/(change.summ[1,2] + change.summ[2,2])*100
change.summ[3,3] <- change.summ[2,3]/(change.summ[1,3] + change.summ[2,3])*100
change.summ[3,4] <- change.summ[2,4]/(change.summ[1,4] + change.summ[2,4])*100
change.summ[3,5] <- change.summ[2,5]/(change.summ[1,5] + change.summ[2,5])*100
change.summ[3,6] <- change.summ[2,6]/(change.summ[1,6] + change.summ[2,6])*100
change.summ[3,7] <- change.summ[2,7]/(change.summ[1,7] + change.summ[2,7])*100
change.summ[3,8] <- change.summ[2,8]/(change.summ[1,8] + change.summ[2,8])*100
change.summ

write.csv(change.summ, file= paste0(sp.name,"_unmasked_pixel_changes_A1A2.csv"), row.names = F)


#masked
currentA1.df <- as.data.frame(rlob1.cur.m, na.rm = T)
colnames(currentA1.df) <- "x"
currentA1.sum <- count(currentA1.df, x)

rcp45.70.df <- as.data.frame(rlob1.45.m, na.rm = T)
colnames(rcp45.70.df) <- "x"
rcp45.70.sum <- count(rcp45.70.df, x)

rcp85.70.df <- as.data.frame(rlob1.85.m, na.rm = T)
colnames(rcp85.70.df) <- "x"
rcp85.70.sum <- count(rcp85.70.df, x)

currentA2.df <- as.data.frame(rlob2.cur.m, na.rm = T)
colnames(currentA2.df) <- "x"
currentA2.sum <- count(currentA2.df, x)

rcp45.70.a2.df <- as.data.frame(rlob2.45.m, na.rm = T)
colnames(rcp45.70.a2.df) <- "x"
rcp45.70.a2.sum <- count(rcp45.70.a2.df, x)

rcp85.70.a2.df <- as.data.frame(rlob2.85.m, na.rm = T)
colnames(rcp85.70.a2.df) <- "x"
rcp85.70.a2.sum <- count(rcp85.70.a2.df, x)

currentA3.df <- as.data.frame(rlob3.cur.m, na.rm = T)
colnames(currentA3.df) <- "x"
currentA3.sum <- count(currentA3.df, x)


change.summ.m <- as.data.frame(cbind(currentA1.sum$x, currentA1.sum$n, currentA2.sum$n, currentA3.sum$n, rcp45.70.sum$n, rcp85.70.sum$n, rcp45.70.a2.sum$n, rcp85.70.a2.sum$n))
colnames(change.summ.m) <- c("State", "Current-A1", "Current-A2", "Current-A3", "RCP4.5-2070-A1", "RCP8.5-2070-A1", "RCP4.5-2070-A2", "RCP8.5-2070-A2")

change.summ.m[3,1] <- "%Suitable Area"
change.summ.m[3,2] <- change.summ.m[2,2]/(change.summ.m[1,2] + change.summ.m[2,2])*100
change.summ.m[3,3] <- change.summ.m[2,3]/(change.summ.m[1,3] + change.summ.m[2,3])*100
change.summ.m[3,4] <- change.summ.m[2,4]/(change.summ.m[1,4] + change.summ.m[2,4])*100
change.summ.m[3,5] <- change.summ.m[2,5]/(change.summ.m[1,5] + change.summ.m[2,5])*100
change.summ.m[3,6] <- change.summ.m[2,6]/(change.summ.m[1,6] + change.summ.m[2,6])*100
change.summ.m[3,7] <- change.summ.m[2,7]/(change.summ.m[1,7] + change.summ.m[2,7])*100
change.summ.m[3,8] <- change.summ.m[2,8]/(change.summ.m[1,8] + change.summ.m[2,8])*100
change.summ.m

write.csv(change.summ.m, file= paste0(sp.name,"_masked_pixel_changes_A1A2.csv"), row.names = F)

#compare geographic range overlaps
#unmasked
# %overlapping areas
geog.range.overlap(enmtools.species(rlob1.cur), enmtools.species(rlob2.cur))
geog.range.overlap(enmtools.species(rlob1.cur), enmtools.species(rlob3.cur))
geog.range.overlap(enmtools.species(rlob2.cur), enmtools.species(rlob3.cur))

# indices
raster.overlap(rlob1.cur, rlob2.cur)
raster.overlap(rlob1.cur, rlob3.cur)
raster.overlap(rlob2.cur, rlob3.cur)

#masked
# %overlapping areas
geog.range.overlap(enmtools.species(rlob1.cur.m), enmtools.species(rlob2.cur.m))
geog.range.overlap(enmtools.species(rlob1.cur.m), enmtools.species(rlob3.cur.m))
geog.range.overlap(enmtools.species(rlob2.cur.m), enmtools.species(rlob3.cur.m))

# indices
raster.overlap(rlob1.cur.m, rlob2.cur.m)
raster.overlap(rlob1.cur.m, rlob3.cur.m)
raster.overlap(rlob2.cur.m, rlob3.cur.m)

#plot
#Approaches, current climate

par(mfrow=c(1,3))
plot(rlob1.cur.m, main = paste0(sp.name," A1"))
plot(rlob2.cur.m, main = paste0(sp.name," A2"))
plot(rlob3.cur.m, main = paste0(sp.name," A3"))

#A1 with future
plot(rlob1.cur.m, main = paste0(sp.name," CURRENT-A1"))
plot(rlob1.45.m, main = paste0(sp.name," RCP4.5-A1"))
plot(rlob1.85.m, main = paste0(sp.name," RCP8.5-A1"))

#A2 with future
plot(rlob2.cur.m, main = paste0(sp.name," CURRENT-A2"))
plot(rlob2.45.m, main = paste0(sp.name," RCP4.5-A2"))
plot(rlob2.85.m, main = paste0(sp.name," RCP8.5-A2"))

### EXTRACTING CLIMATE DATA ####
setwd("Envi")

#get PH shapefile
phl.shp <- getData('GADM', country = 'PHL', level=0)

#read and crop altitude data
altitude <- raster("D:/DATA/GIS_Layers/wc2.1_30s_elev/wc2.1_30s_elev.tif")
alt.PH <- crop(altitude, phl.shp)
alt.PH <- mask(alt.PH, phl.shp) 
plot(alt.PH)

#extract data for each data point 
sppoints <- datapoints[,-1]

envdata <- stack(modelclimEnv, alt.PH)
names(envdata)[12] <- "altitude"
climate.df <- raster::extract(envdata, sppoints, df = T)
climate.df <- cbind(datapoints, climate.df)
climate.df <- climate.df[,-4]
write.csv(climate.df, file = "maxent_clean2_Tetrastigma_Rafflesia_30s_envdata.csv", row.names = F)

#extract data for each binary result
#masked
## R. speciosa
#current
spe <- as.data.frame(rspe1.cur.m, na.rm = T, xy=T)
colnames(spe) <- c("lon", "lat", "bin")
spe <- filter(spe, bin==1)
spe <- cbind(spe$lon, spe$lat)

spe.df <- raster::extract(envdata, spe, method = "simple", na.rm = T)
spe.df <- as.data.frame(spe.df)
spe.df$Scenario <- rep("Current", dim(spe.df)[1])

#RCP4.5
spe45 <- as.data.frame(rspe1.45.m, na.rm = T, xy=T)
colnames(spe45) <- c("lon", "lat", "bin")
spe45 <- filter(spe45, bin==1)
spe45 <- cbind(spe45$lon, spe45$lat)

spe45.df <- raster::extract(envdata, spe45, method = "simple", na.rm = T)
spe45.df <- as.data.frame(spe45.df)
spe45.df$Scenario <- rep("RCP4.5", dim(spe45.df)[1])


#RCP 8.5
spe85 <- as.data.frame(rspe1.85.m, na.rm = T, xy=T)
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
lag <- as.data.frame(rlag1.cur.m, na.rm = T, xy=T)
colnames(lag) <- c("lon", "lat", "bin")
lag <- filter(lag, bin==1)
lag <- cbind(lag$lon, lag$lat)

lag.df <- raster::extract(envdata, lag, method = "simple", na.rm = T)
lag.df <- as.data.frame(lag.df)
lag.df$Scenario <- rep("Current", dim(lag.df)[1])

#RCP4.5
lag45 <- as.data.frame(rlag1.45.m, na.rm = T, xy=T)
colnames(lag45) <- c("lon", "lat", "bin")
lag45 <- filter(lag45, bin==1)
lag45 <- cbind(lag45$lon, lag45$lat)

lag45.df <- raster::extract(envdata, lag45, method = "simple", na.rm = T)
lag45.df <- as.data.frame(lag45.df)
lag45.df$Scenario <- rep("RCP4.5", dim(lag45.df)[1])


#RCP 8.5
lag85 <- as.data.frame(rlag1.85.m, na.rm = T, xy=T)
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
lob <- as.data.frame(rlob1.cur.m, na.rm = T, xy=T)
colnames(lob) <- c("lon", "lat", "bin")
lob <- filter(lob, bin==1)
lob <- cbind(lob$lon, lob$lat)

lob.df <- raster::extract(envdata, lob, method = "simple", na.rm = T)
lob.df <- as.data.frame(lob.df)
lob.df$Scenario <- rep("Current", dim(lob.df)[1])

#RCP4.5
lob45 <- as.data.frame(rlob1.45.m, na.rm = T, xy=T)
colnames(lob45) <- c("lon", "lat", "bin")
lob45 <- filter(lob45, bin==1)
lob45 <- cbind(lob45$lon, lob45$lat)

lob45.df <- raster::extract(envdata, lob45, method = "simple", na.rm = T)
lob45.df <- as.data.frame(lob45.df)
lob45.df$Scenario <- rep("RCP4.5", dim(lob45.df)[1])


#RCP 8.5
lob85 <- as.data.frame(rlob1.85.m, na.rm = T, xy=T)
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
all.alt <- rbind(allspe, alllag, alllob)
all.alt$Species <- as.factor(all.alt$Species)
all.alt$Scenario <- as.factor(all.alt$Scenario)
summary(all.alt)

write.csv(all.alt, file = "All_Rafflesia_spp_suitable_masked_envdata.csv", row.names = F)

#boxplots
#altitude
altboxgg <- ggplot() + 
  geom_boxplot(aes(y = altitude, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Altitude (masl)")
ggsave(altboxgg, file="Plots/Boxplot-alt-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)

tempboxgg <- ggplot() + 
  geom_boxplot(aes(y = bio1/10, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Mean Annual Temp (°C)")
ggsave(tempboxgg, file="Plots/Boxplot-matemp-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)

precboxgg <- ggplot() + 
  geom_boxplot(aes(y = bio12, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Mean Annual Precipitation (mm/yr)")
ggsave(precboxgg, file="Plots/Boxplot-prec-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)


tempseaboxgg <- ggplot() + 
  geom_boxplot(aes(y = bio4, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Temperature Seasonality (sd)")
ggsave(tempseaboxgg, file="Plots/Boxplot-tempsea-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)

precseaboxgg <- ggplot() + 
  geom_boxplot(aes(y = bio15, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Precipitation Seasonality (cv)")
ggsave(precseaboxgg, file="Plots/Boxplot-precsea-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)

precoldboxgg <- ggplot() + 
  geom_boxplot(aes(y = bio19, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Precipitation of Coldest Quarter (mm/quarter)")
ggsave(precoldboxgg, file="Plots/Boxplot-precoldq-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)

#soils
bdboxgg <- ggplot() + 
  geom_boxplot(aes(y = bulk_density, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Soil Bulk Density (kg/m^3)")
ggsave(bdboxgg, file="Plots/Boxplot-bulkdens-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)

cecboxgg <- ggplot() + 
  geom_boxplot(aes(y = cec, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Soil Cation Exchange Capacity (cmolc/kg)")
ggsave(cecboxgg, file="Plots/Boxplot-cec-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)

soilphgg <- ggplot() + 
  geom_boxplot(aes(y = soilph/10, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Soil pH")
ggsave(soilphgg, file="Plots/Boxplot-soilph-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)

claygg <- ggplot() + 
  geom_boxplot(aes(y = clay, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Soil Clay Content (%)")
ggsave(claygg, file="Plots/Boxplot-clay-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)

sandgg <- ggplot() + 
  geom_boxplot(aes(y = sand, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Soil Sand Content (%)")
ggsave(sandgg, file="Plots/Boxplot-sand-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)

siltgg <- ggplot() + 
  geom_boxplot(aes(y = silt, x = Species, fill=Scenario), 
               data=all.alt , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Soil Silt Content (%)")
ggsave(siltgg, file="Plots/Boxplot-silt-Rafflesia_spp_masked.png", width=19.89, height=15, units="cm", dpi=300)

#LINEAR MODELS ####
# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

#Speciosa

#check and transform distribution of altitude data
hist(allspe$altitude)
hist(sqrt(allspe$altitude))
allspe$alt.sqrt <- sqrt(allspe$altitude)
allspe <- na.omit(allspe)

a1 <- lm(alt.sqrt ~ Scenario, data = allspe)
summ(a1)
a2 <- aov(a1)

a3 <- TukeyHSD(a2, 'Scenario', conf.level = 0.95)
a3
plot(a3, las=1 , col="brown")


# Generate groups from Tukey test
a3lab <- generate_label_df(a3 , 'Scenario')
a3lab


#Lagascae

#check and transform distribution of altitude data
hist(alllag$altitude)
hist(sqrt(alllag$altitude))
alllag$alt.sqrt <- sqrt(alllag$altitude)
alllag <- na.omit(alllag)

a4 <- lm(alt.sqrt ~ Scenario, data = alllag)
summ(a4)
a5 <- aov(a4)

a6 <- TukeyHSD(a5, 'Scenario', conf.level = 0.95)
a6
plot(a6, las=1 , col="brown")


# Generate groups from Tukey test
a6lab <- generate_label_df(a6 , 'Scenario')
a6lab

#Lobata

#check and transform distribution of altitude data

alllob2 <- alllob %>% slice_sample(n=10000)

hist(alllob2$altitude)
hist(sqrt(alllob2$altitude))
alllob2$alt.sqrt <- sqrt(alllob2$altitude)
alllob2 <- na.omit(alllob2)

a7 <- lm(alt.sqrt ~ Scenario, data = alllob2)
summ(a7)
a8 <- aov(a7)

a9 <- TukeyHSD(a8, 'Scenario', conf.level = 0.95)
a9
plot(a9, las=1 , col="brown")


# Generate groups from Tukey test
a9lab <- generate_label_df(a9 , 'Scenario')
a9lab



#### PLOTS ####
#pal <- colorRampPalette(c("grey", "darkgreen", "red"))
#cuts <- c(0,1,2,3)

#R. speciosa

windows()
par(mfrow = c(1,3))
plot(rs1)
plot(rs2)
plot(rs3)

#R. lagascae
par(mfrow = c(1,3))
plot(rlag1)
plot(rlag2)
plot(rlag3)

#R. lobata
par(mfrow = c(1,3))
plot(rlob1)
plot(rlob2)
plot(rlob3)


### END OF CODE ####