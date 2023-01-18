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
library(xlsx)

setwd("D:/PROJECTS/Rafflesia_SDM/Altitude_analysis")

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

alldata <- read_excel("All_Rafflesia_spp_suitable_masked_envdata.xlsx") #entire PH
alldata <- read.csv("D:/PROJECTS/Rafflesia_SDM/Island_analysis/All_Rafflesia_spp_a1_suitable_masked_native_ranges_envdata.csv") #native habitats only
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

visreg(a4)

#Lobata

#check and transform distribution of altitude data
hist(alllob$altitude)
hist(sqrt(alllob$altitude))
alllob$alt.sqrt <- sqrt(alllob$altitude)
alllob <- na.omit(alllob)

alllob2 <- alllob
alllob2 <- alllob %>% slice_sample(n=10000)

hist(alllob2$altitude)
hist(sqrt(alllob2$altitude))
alllob2$alt.sqrt <- sqrt(alllob2$altitude)
alllob2 <- na.omit(alllob2)


a7 <- lm(alt.sqrt ~ Scenario, data = alllob)
summ(a7)
a8 <- aov(a7)

a9 <- TukeyHSD(a8, 'Scenario', conf.level = 0.95)
a9
plot(a9, las=1 , col="brown")


# Generate groups from Tukey test
a9lab <- generate_label_df(a9 , 'Scenario')
a9lab

visreg(a7)


#### Multiple regression models ####
# This is needed to verify the significance of the scenario changes, especially for R. lobata
# Need to sample randomly from the data set and test the significance each time

### Looping 1000 glm's ###

options(width = 80)
n <- 1000

df <- alllob #alllob, allspe

for (i in 1:n) {
  
  #sample 5,000 points from the df of each species
  sam1 <- slice_sample(.data = df, n=5000, replace=T)
  
  #linear model
  lm1 <- lm(alt.sqrt ~ Scenario, data = sam1)
  s <- summ(lm1)
  
  #Tukey
  lm2 <- aov(lm1)
  lm3 <- TukeyHSD(lm2, 'Scenario', conf.level = 0.95)
  lm3
  
  
  #extract and combine coefficients
  c <- coefficients(lm1) 
  
  if (i==1){
    coef <- c
  } else {
    coef <- as.data.frame(rbind(coef, c))  
  }
  

  #extract and combine p-values
  pv <- lm3$Scenario[10:12]
  
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

lm3lab <- generate_label_df(lm3 , 'Scenario')
lm3lab

#rename column name of R2 data frame
#colnames(R2.df) <- c(a1$variable)

#rename column name of P-values data frame
colnames(pv.df) <- c("RCP4.5-Current", "RCP8.5-Current", "RCP8.5-RCP4.5")

#write it out
write.xlsx(coef, file = "Rspeciosa_bin_5k_tukey.xlsx", sheetName = "Coefficients", col.names = T, row.names = F, append = T)
#write.xlsx(R2.df, file = "Rlagascae_bin_5k_tukey.xlsx", sheetName = "Partial R2", col.names = T, row.names = F, append = T)
write.xlsx(pv.df, file = "Rspeciosa_bin_5k_tukey.xlsx", sheetName = "P-values", col.names = T, row.names = F, append = T)


#getting quantiles
#coefficients
qc.int <- quantile(coef$`(Intercept)`, c(0.025, 0.975))
qc.rcp45 <- quantile(coef$ScenarioRCP4.5, c(0.025, 0.975))
qc.rcp85 <- quantile(coef$ScenarioRCP8.5, c(0.025, 0.975))

qc.df <- as.data.frame(rbind(qc.int, qc.rcp45, qc.rcp85))
qc.df$Var <- c("Intercept", "RCP45", "RCP85")

#P-values
qp.int <- quantile(pv.df$Intercept, c(0.025, 0.5, 0.975))
qp.r45 <- quantile(pv.df$RCP45, c(0.025, 0.5, 0.975))
qp.r85 <- quantile(pv.df$RCP85, c(0.025, 0.5, 0.975))

qp.df <- as.data.frame(rbind(qp.int, qp.r45, qp.r85))
qp.df$Var <- c("RCP4.5-Current", "RCP8.5-Current", "RCP8.5-RCP4.5")


#partial R-squared
#qr.sce <- quantile(R2.df$Scenario, c(0.025, 0.5, 0.975))


#write it out
write.xlsx(qc.df, file = "Rspeciosa_bin_5k_tukey.xlsx", sheetName = "Quant-Coeffs", col.names = T, row.names = F, append = T)
write.xlsx(qp.df, file = "Rspeciosa_bin_5k_tukey.xlsx", sheetName = "Quant-PVals", col.names = T, row.names = F, append = T)


#### BOXPLOTS ####

#boxplot
altboxgg <- ggplot() + 
  geom_boxplot(aes(y = altitude, x = Species, fill=Scenario), 
               data=alldata , outlier.shape = 1, outlier.size = 1) +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Altitude (masl)")
ggsave(altboxgg, file="Boxplot-alt-spp.png", width=19.89, height=15, units="cm", dpi=300)

#violin plot 
altgg <- ggplot(aes(y = altitude, x = Species, fill=Scenario), data=alldata) +
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  scale_fill_viridis(discrete=T, option = "E", name="", direction = -1) +
  theme_ipsum() +
  ylab("Altitude (masl)")

ggsave(altgg, file="Vioplot-alt-spp.png", width=19.89, height=15, units="cm", dpi=300)

