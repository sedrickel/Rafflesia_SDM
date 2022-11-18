require(tidyverse)
library(forcats)
library(hrbrthemes)
library(viridis)
library(multcompView)
library(visreg)
library(jtools)
library(readxl)
library(rsq)
library(performance)

# Function to group the treatments that are not sig. different together.
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

### Input data ####

alldata <- read_excel("All_Rafflesia_spp_suitable_masked_envdata.xlsx")
summary(alldata)

alldata$Species <- as.factor(alldata$Species)
alldata$Scenario <- as.factor(alldata$Scenario)

#separate per species 
allspe <- filter(alldata, Species == "Rspe") #speciosa
alllag <- filter(alldata, Species == "Rlag") #lagascae
alllob <- filter(alldata, Species == "Rlob") #lobata

summary(allspe)
summary(alllag)
summary(alllob)

#### LINEAR MODELS ####

### Single models ####
# Using all data points under one linear model only #

#### Speciosa ####

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

visreg(a1)

# Generate groups from Tukey test
a3lab <- generate_label_df(a3, 'Scenario')
a3lab

#partial R-squared

r2.all <- rsq.partial(objF = a1, objR = NULL, adj = F, type = "v")
r2.all


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

alllob2 <- alllob
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

visreg(a7)




### Multiple regression models ###
# This is needed to verify the significance of the scenario changes, especially for R. lobata
# Need to sample randomly from the data set and test the significance each time

### Looping 1000 glm's ###

options(width = 80)
n <- 1000

for (i in 1:n) {
  
  #sample 10,000 points from the df
  sam1 <- slice_sample(.data = filter(ras8616.ext.1,  Class == 1), n=5000, replace=T)
  sam0 <- slice_sample(.data = filter(ras8616.ext.1,  Class == 0), n=5000, replace=T)
  
  sam1 <- rbind(sam1, sam0)
  
  #GLM
  glm1 <- glm(Class ~ MATemp + Aridity.sqrt + Slope.sqrt + FirePres.all, sam1, family='binomial')
  s <- summ(glm1)
  
  #partial R-squared
  a1 <- rsq.partial(objF = glm1, objR = NULL, adj = F, type = "v")
  
  #extract and combine coefficients
  c <- coefficients(glm1) 
  
  if (i==1){
    coef <- c
  } else {
    coef <- as.data.frame(rbind(coef, c))  
  }
  
  #extract and combine partial R2
  r2 <- a1$partial.rsq
  
  if (i==1){
    R2.df <- r2
  } else {
    R2.df <- as.data.frame(rbind(R2.df, r2))  
  }
  
  #extract and combine p-values
  pv <- s$coeftable[16:20]
  
  if (i==1){
    pv.df <- pv
  } else {
    pv.df <- as.data.frame(rbind(pv.df, pv))  
  }
  
  #Sys.sleep(0.01)
  #print(i)
  
  extra <- nchar('||100%')
  width <- options()$width
  step <- round(i / n * (width - extra))
  text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
                  strrep(' ', width - step - extra), round(i / n * 100))
  cat(text)
  Sys.sleep(0.5)
  cat(if (i == n) '\n' else '\014')
  
  
}
# end loop #

#rename column name of R2 data frame
colnames(R2.df) <- c(a1$variable)

#rename column name of P-values data frame
colnames(pv.df) <- c("Intercept", a1$variable)

#write it out
write.xlsx(coef, file = "GLM_Prk_Wdl_bin_10k.xlsx", sheetName = "Coefficients", col.names = T, row.names = F, append = T)
write.xlsx(R2.df, file = "GLM_Prk_Wdl_bin_10k.xlsx", sheetName = "Partial R2", col.names = T, row.names = F, append = T)
write.xlsx(pv.df, file = "GLM_Prk_Wdl_bin_10k.xlsx", sheetName = "P-values", col.names = T, row.names = F, append = T)


#getting quantiles
#coefficients
qc.int <- quantile(coef$`(Intercept)`, c(0.025, 0.975))
qc.temp <- quantile(coef$MATemp, c(0.025, 0.975))
qc.arid <- quantile(coef$Aridity.sqrt, c(0.025, 0.975))
qc.slope <- quantile(coef$Slope.sqrt, c(0.025, 0.975))
qc.fire <- quantile(coef$FirePres.all, c(0.025, 0.975))


qc.df <- as.data.frame(rbind(qc.int, qc.temp, qc.arid, qc.slope, qc.fire))
qc.df$Var <- c("Intercept", a1$variable)

#partial R-squared
qr.temp <- quantile(R2.df $MATemp, c(0.025, 0.5, 0.975))
qr.arid <- quantile(R2.df$Aridity.sqrt, c(0.025, 0.5, 0.975))
qr.slope <- quantile(R2.df$Slope.sqrt, c(0.025, 0.5, 0.975))
qr.fire <- quantile(R2.df$FirePres.all, c(0.025, 0.5, 0.975))


qr.df <- as.data.frame(rbind(qr.temp, qr.arid, qr.slope, qr.fire))
qr.df$Var <- c(a1$variable)

#write it out
write.xlsx(qc.df, file = "GLM_Prk_Wdl_bin_10k.xlsx", sheetName = "Quant-Coeffs", col.names = T, row.names = F, append = T)
write.xlsx(qr.df, file = "GLM_Prk_Wdl_bin_10k.xlsx", sheetName = "Quant-R2", col.names = T, row.names = F, append = T)

table(ras8616.ext.1$Class)
summ(glm.all)
r2.all


