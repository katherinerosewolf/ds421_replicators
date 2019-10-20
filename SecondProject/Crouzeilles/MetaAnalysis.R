#### Script provided by Renato Crouzeilles (International Institute for Sustainability) - renatocrouzeilles@gmail.com

##########################################################################################################################

# working directory
library(here)

setwd(here("SecondProject", "Crouzeilles"))

getwd()

## Read data for RESTORED SYSTEMS
data<-read.table("Meta_Analysis.txt", header=T, sep="\t")
data<-subset(data, data$Disturbance_Conversion=="secondary" | data$Disturbance_Conversion=="selectively_logged")

## Plants
Plants_all<-subset(data, data$Plants==1)
outlyers<-boxplot(Plants_all$RR)
sort(outlyers$out) # -1.7918  1.2040
Plants<-subset(Plants_all, Plants_all$RR<1.2040 & Plants_all$RR>-1.7918)
hist(Plants$RR, breaks=30)
qqnorm(Plants$RR)
qqline(Plants$RR)

## Invertebrates
Invertebrates_all<-subset(data, data$Invertebrates==1)
outlyers<-boxplot(Invertebrates_all$RR)
sort(outlyers$out) # -1.7065  1.2040
Invertebrates<-subset(Invertebrates_all, Invertebrates_all$RR<1.2040 & Invertebrates_all$RR>-1.7065)
hist(Invertebrates$RR, breaks=30)
qqnorm(Invertebrates$RR)
qqline(Invertebrates$RR)

## Birds
Birds_all<-subset(data, data$Birds==1)
outlyers<-boxplot(Birds_all$RR)
sort(outlyers$out) # -0.8979  0.6739
Birds<-subset(Birds_all, Birds_all$RR<0.6739 & Birds_all$RR>-0.8979)
hist(Birds$RR, breaks=30)
qqnorm(Birds$RR)
qqline(Birds$RR)

## Herpetofaua
Herpeto_all<-subset(data, data$Herpetofauna==1)
outlyers<-boxplot(Herpeto_all$RR)
sort(outlyers$out) # -1.6094  1.1337
Herpeto<-subset(Herpeto_all, Herpeto_all$RR<1.1337 & Herpeto_all$RR>-1.6094)
hist(Herpeto$RR, breaks=30)
qqnorm(Herpeto$RR)
qqline(Herpeto$RR)

## Mammals
Mammals_all<-subset(data, data$Mammals==1)
outlyers<-boxplot(Mammals_all$RR)
sort(outlyers$out) # -1.6094  1.3138
Mammals<-subset(Mammals_all, Mammals_all$RR<1.3138 & Mammals_all$RR>-1.6094)
hist(Mammals$RR, breaks=30)
qqnorm(Mammals$RR)
qqline(Mammals$RR)

## Cover
Cover_all<-subset(data, data$Cover==1)
outlyers<-boxplot(Cover_all$RR)
sort(outlyers$out) # -1.8003  1.0647
Cover<-subset(Cover_all, Cover_all$RR<1.0647 & Cover_all$RR>-1.8003)
hist(Cover$RR, breaks=30)
qqnorm(Cover$RR)
qqline(Cover$RR)

## Density
Density_all<-subset(data, data$Density==1)
outlyers<-boxplot(Density_all$RR)
sort(outlyers$out) # -1.5782  1.3740
Density<-subset(Density_all, Density_all$RR<1.3740 & Density_all$RR>-1.5782)
hist(Density$RR, breaks=30)
qqnorm(Density$RR)
qqline(Density$RR)

## Biomass
Biomass_all<-subset(data, data$Biomass==1)
outlyers<-boxplot(Biomass_all$RR)
sort(outlyers$out) # -2.1972  1.2585
Biomass<-subset(Biomass_all, Biomass_all$RR<1.2585 & Biomass_all$RR>-2.1972)
hist(Biomass$RR, breaks=30)
qqnorm(Biomass$RR)
qqline(Biomass$RR)

## Height
Height_all<-subset(data, data$Height==1)
outlyers<-boxplot(Height_all$RR)
sort(outlyers$out) # -1.5275
Height<-subset(Height_all, Height_all$RR>-1.5275)
hist(Height$RR, breaks=30)
qqnorm(Height$RR)
qqline(Height$RR)

## Litter
Litter_all<-subset(data, data$Litter==1)
outlyers<-boxplot(Litter_all$RR)
sort(outlyers$out)
Litter<-na.omit(Litter_all)
hist(Litter$RR, breaks=30)
qqnorm(Litter$RR)
qqline(Litter$RR)

## Sample 1 comparisons per study landscape
randomRows = function(df,n){
  return(df[sample(nrow(df),n),])
}
library("plyr")

## Bootstrap for meta-analysis ALL THESE ORIGINALLY 10000 RUNS
## Plants
b.Plants <- c()
for (i in 1:1000){
  # First take 1 comparisons per study landscape
  bsample <- ddply(Plants,.(Site),randomRows,1)
  #now calculate the bootstrap estimate for mean RR
  bestimate <- mean(bsample$RR)
  b.Plants <- c(b.Plants,bestimate)}

## Invertebrates
b.Invertebrates <- c()
for (i in 1:1000){
  bsample <- ddply(Invertebrates,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Invertebrates <- c(b.Invertebrates,bestimate)}

## Birds
b.Birds <- c()
for (i in 1:1000){
  bsample <- ddply(Birds,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Birds <- c(b.Birds,bestimate)}

## Herpeto
b.Herpeto <- c()
for (i in 1:1000){
  bsample <- ddply(Herpeto,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Herpeto <- c(b.Herpeto,bestimate)}

## Mammals
b.Mammals <- c()
for (i in 1:1000){
  bsample <- ddply(Mammals,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Mammals <- c(b.Mammals,bestimate)}

## Cover
b.Cover <- c()
for (i in 1:1000){
  bsample <- ddply(Cover,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Cover <- c(b.Cover,bestimate)}

## Density
b.Density <- c()
for (i in 1:1000){
  bsample <- ddply(Density,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Density <- c(b.Density,bestimate)}

## Biomass
b.Biomass <- c()
for (i in 1:1000){
  bsample <- ddply(Biomass,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Biomass <- c(b.Biomass,bestimate)}

## Height
b.Height <- c()
for (i in 1:1000){
  bsample <- ddply(Height,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Height <- c(b.Height,bestimate)}

## Litter
b.Litter <- c()
for (i in 1:1000){
  bsample <- ddply(Litter,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Litter <- c(b.Litter,bestimate)}

## Join results
meanRR_Restored<-matrix(c(b.Plants,b.Invertebrates,b.Birds,b.Herpeto,b.Mammals,b.Cover,b.Density,b.Biomass,b.Height,b.Litter),,10)

write.table(meanRR_Restored, "Meta_Analysis_Restored_systems_results.txt", quote=F, sep="\t")

##########################################################################################################################
## Read data for DEGRADED SYSTEMS
data<-read.table("Meta_Analysis.txt", header=T, sep="\t")
data<-subset(data, data$Disturbance_Conversion=="shaded_plantation" | data$Disturbance_Conversion=="agriculture" | data$Disturbance_Conversion=="pasture" | data$Disturbance_Conversion=="plantation" | data$Disturbance_Conversion=="open")

## Plants
Plants_all<-subset(data, data$Plants==1)
outlyers<-boxplot(Plants_all$RR)
sort(outlyers$out) # -2.7314  1.4069
Plants<-subset(Plants_all, Plants_all$RR<1.4069 & Plants_all$RR>-2.7314)
hist(Plants$RR, breaks=30)
qqnorm(Plants$RR)
qqline(Plants$RR)

## Invertebrates
Invertebrates_all<-subset(data, data$Invertebrates==1)
outlyers<-boxplot(Invertebrates_all$RR)
sort(outlyers$out) # -2.5530  2.0794
Invertebrates<-subset(Invertebrates_all, Invertebrates_all$RR<2.0794 & Invertebrates_all$RR>-2.5530)
hist(Invertebrates$RR, breaks=30)
qqnorm(Invertebrates$RR)
qqline(Invertebrates$RR)

## Birds
Birds_all<-subset(data, data$Birds==1)
outlyers<-boxplot(Birds_all$RR)
sort(outlyers$out) # -2.7555  1.6529
Birds<-subset(Birds_all, Birds_all$RR<1.6529 & Birds_all$RR>-2.7555)
hist(Birds$RR, breaks=30)
qqnorm(Birds$RR)
qqline(Birds$RR)

## Herpetofaua
Herpeto_all<-subset(data, data$Herpetofauna==1)
outlyers<-boxplot(Herpeto_all$RR)
sort(outlyers$out) # -2.2073
Herpeto<-subset(Herpeto_all, Herpeto_all$RR>-2.2073)
hist(Herpeto$RR, breaks=30)
qqnorm(Herpeto$RR)
qqline(Herpeto$RR)

## Mammals
Mammals_all<-subset(data, data$Mammals==1)
outlyers<-boxplot(Mammals_all$RR)
sort(outlyers$out) # -2.7300
Mammals<-subset(Mammals_all, Mammals_all$RR>-2.7300)
hist(Mammals$RR, breaks=30)
qqnorm(Mammals$RR)
qqline(Mammals$RR)

## Cover
Cover_all<-subset(data, data$Cover==1)
outlyers<-boxplot(Cover_all$RR)
sort(outlyers$out) # -4.2627
Cover<-subset(Cover_all, Cover_all$RR>-4.2627)
hist(Cover$RR, breaks=30)
qqnorm(Cover$RR)
qqline(Cover$RR)

## Density
Density_all<-subset(data, data$Density==1)
outlyers<-boxplot(Density_all$RR)
sort(outlyers$out) # -3.6630
Density<-subset(Density_all, Density_all$RR>-3.6630)
hist(Density$RR, breaks=30)
qqnorm(Density$RR)
qqline(Density$RR)

## Biomass
Biomass_all<-subset(data, data$Biomass==1)
outlyers<-boxplot(Biomass_all$RR)
sort(outlyers$out) # -6.6724
Biomass<-subset(Biomass_all, Biomass_all$RR>-6.6724)
hist(Biomass$RR, breaks=30)
qqnorm(Biomass$RR)
qqline(Biomass$RR)

## Height
Height_all<-subset(data, data$Height==1)
outlyers<-boxplot(Height_all$RR)
sort(outlyers$out) # -3.1499
Height<-subset(Height_all, Height_all$RR>-3.1499)
hist(Height$RR, breaks=30)
qqnorm(Height$RR)
qqline(Height$RR)

## Litter
Litter_all<-subset(data, data$Litter==1)
outlyers<-boxplot(Litter_all$RR)
sort(outlyers$out)
Litter<-Litter_all
hist(Litter$RR, breaks=30)
qqnorm(Litter$RR)
qqline(Litter$RR)

## Sample 1 comparisons per study landscape
randomRows = function(df,n){
  return(df[sample(nrow(df),n),])
}
library("plyr")

## Bootstrap for meta-analysis AGAIN 10000 TIMES
## Plants
b.Plants <- c()
for (i in 1:1000){
  # First take 1 comparisons per study landscape
  bsample <- ddply(Plants,.(Site),randomRows,1)
  #now calculate the bootstrap estimate for mean RR
  bestimate <- mean(bsample$RR)
  b.Plants <- c(b.Plants,bestimate)}

## Invertebrates
b.Invertebrates <- c()
for (i in 1:1000){
  bsample <- ddply(Invertebrates,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Invertebrates <- c(b.Invertebrates,bestimate)}

## Birds
b.Birds <- c()
for (i in 1:1000){
  bsample <- ddply(Birds,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Birds <- c(b.Birds,bestimate)}

## Herpeto
b.Herpeto <- c()
for (i in 1:1000){
  bsample <- ddply(Herpeto,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Herpeto <- c(b.Herpeto,bestimate)}

## Mammals
b.Mammals <- c()
for (i in 1:1000){
  bsample <- ddply(Mammals,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Mammals <- c(b.Mammals,bestimate)}

## Cover
b.Cover <- c()
for (i in 1:1000){
  bsample <- ddply(Cover,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Cover <- c(b.Cover,bestimate)}

## Density
b.Density <- c()
for (i in 1:1000){
  bsample <- ddply(Density,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Density <- c(b.Density,bestimate)}

## Biomass
b.Biomass <- c()
for (i in 1:1000){
  bsample <- ddply(Biomass,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Biomass <- c(b.Biomass,bestimate)}

## Height
b.Height <- c()
for (i in 1:1000){
  bsample <- ddply(Height,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Height <- c(b.Height,bestimate)}

## Litter
b.Litter <- c()
for (i in 1:1000){
  bsample <- ddply(Litter,.(Site),randomRows,1)
  bestimate <- mean(bsample$RR)
  b.Litter <- c(b.Litter,bestimate)}

## Join results
meanRR_Degraded<-matrix(c(b.Plants,b.Invertebrates,b.Birds,b.Herpeto,b.Mammals,b.Cover,b.Density,b.Biomass,b.Height,b.Litter),,10)

write.table(meanRR_Degraded, "Meta_Analysis_Degraded_systems_results.txt", quote=F, sep="\t")

#### New code

##Packages
library(tidyverse)

##Setting up data tables
#Degraded
RR_deg <- meanRR_Degraded %>%
  as_tibble() %>%
  select("Deg_Plants" = V1, "Deg_Invertebrates" = V2, "Deg_Birds" = V3, "Deg_Herpeto" = V4, "Deg_Mammals" = V5, "Deg_Cover" = V6, "Deg_Density" = V7, "Deg_Biomass" = V8, "Deg_Height" = V9, "Deg_Litter" = V10)

deg_bio <- RR_deg %>%
  select(Deg_Plants, Deg_Invertebrates, Deg_Birds, Deg_Herpeto, Deg_Mammals)

deg_veg <- RR_deg %>%
  select(Deg_Cover, Deg_Density, Deg_Biomass, Deg_Height, Deg_Litter)

#Restored
RR_res <- meanRR_Restored %>%
  as_tibble() %>%
  select("Res_Plants" = V1, "Res_Invertebrates" = V2, "Res_Birds" = V3, "Res_Herpeto" = V4, "Res_Mammals" = V5, "Res_Cover" = V6, "Res_Density" = V7, "Res_Biomass" = V8, "Res_Height" = V9, "Res_Litter" = V10)

res_bio <- RR_res %>%
  select(Res_Plants, Res_Invertebrates, Res_Birds, Res_Herpeto, Res_Mammals)

res_veg <- RR_res %>%
  select(Res_Cover, Res_Density, Res_Biomass, Res_Height, Res_Litter)

# # Combine by metric
# biodiversity <- bind_cols(deg_bio, res_bio)
# 
# vegetation <- bind_cols(deg_veg, res_veg)


# Plots

pdf("biodiversity.pdf")

biodiverse_deg_plot <- ggplot(stack(deg_bio), aes(x=ind, y=values)) +
  geom_boxplot(notch = TRUE) +
  coord_flip() + 
  ylab("Response Ratio") + 
  xlab("Biodiversity Metric") + 
  theme_minimal() +
  scale_y_continuous(limits = c(-1.7, 0.1))

biodiverse_res_plot <- ggplot(stack(res_bio), aes(x=ind, y=values)) +
  geom_boxplot(notch = TRUE) +
  coord_flip() + 
  ylab("Response Ratio") + 
  xlab("Biodiversity Metric") + 
  theme_minimal() +
  scale_y_continuous(limits = c(-1.7, 0.1))

print(biodiverse_deg_plot)
print(biodiverse_res_plot)
dev.off()

pdf("vegetation.pdf")

veg_deg_plot <- ggplot(stack(deg_veg), aes(x=ind, y=values)) + 
  geom_boxplot(notch = TRUE) + 
  coord_flip() + 
  ylab("Response Ratio") + 
  xlab("Vegetation Metric") + 
  theme_minimal() +
  scale_y_continuous(limits = c(-1.7, 0.1))

veg_res_plot <- ggplot(stack(res_veg), aes(x=ind, y=values)) + 
  geom_boxplot(notch = TRUE) + 
  coord_flip() + 
  ylab("Response Ratio") + 
  xlab("Vegetation Metric") + 
  theme_minimal() +
  scale_y_continuous(limits = c(-1.7, 0.1))
 
print(veg_deg_plot)
print(veg_res_plot)
dev.off()
