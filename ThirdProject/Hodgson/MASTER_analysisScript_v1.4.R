## this script analyses the counts from the #EpicDuckChallenge. It takes as input only
 # the master data-file, "MASTER_AllCountData.csv".

setwd("~/ds421_replicators/ThirdProject/Hodgson")

data <- read.csv("MASTER_AllCountData.csv")
data$Date <- as.Date(data$Date, format = "%d/%m/%Y")
data$Height_m[data$Height_m == "n/a"] <- NA
data$Start_time[data$Start_time== "n/a"] <- NA
data$End_time[data$End_time== "n/a"] <- NA
data$Elapsed_time[data$Elapsed_time== "n/a"] <- NA

###  The following are the analyses detailed in OSF Pre-registration at "osf.io/7unqw"

library(MASS)

## Question One: What is the absolute difference between the true number of animals in a
 # colony and the counted number of animals in the colony?
autoCount10s <- subset(data, data$Count_type == "Auto", drop = TRUE)
autoCount10s <- subset(autoCount10s, autoCount10s$percent_Input == 10, drop = TRUE)
## 10% training data was chosen as the level of training to present; disregard other
 # estimates for different amounts of semi-automated detection counts.

counts <- rbind(subset(data, data$Count_type == "Ground", drop = TRUE),
                subset(data, data$Count_type == "UAV_manual", drop = TRUE),
                autoCount10s)
counts <- droplevels(counts)
colonies <- droplevels(subset(data, data$Count_type == "Retreived_skewers"))

counts$trueCounts <- c(rep(NaN, nrow(counts)))
for (i in 1:10) {
    counts$trueCounts[counts$Colony == i] <- colonies$Count[colonies$Colony == i][1]
}
counts$absdiff <- abs(counts$Count - counts$trueCounts)
counts$tech <- paste(counts$Count_type, counts$Height_m, sep = "_")

m1 <- with(counts, glm(absdiff ~ tech + as.factor(Colony), family = "quasipoisson"))
summary(m1)
summary(aov(m1))

## Results so far indicate that we have strong significant effects of tech and colony on
 # absolute error of counts.

# Figures for Question One

# with(counts, boxplot(absdiff ~ tech + as.factor(Colony)))  ## Too crowded. Split by colony.

pdf(file = "Q1_abs_error_v2.pdf", width = 8, height = 8)
par(mfrow = c(3,2), oma = c(0,0,0,0), mar = c(2.5,4.5,1,1))
#for(i in c(2,4,3,10,1,7,8,6,9,5)) { ## all colonies
for(i in c(2,4,3,1,7,6)) { ## in-focus colonies
    activeColony <- droplevels(subset(counts, counts$Colony == i))
    activeColony$tech <- factor(activeColony$tech,
                                levels = c("Ground_NA", "UAV_manual_30", "UAV_manual_60",
                                           "UAV_manual_90", "UAV_manual_120",
                                           "Auto_30", "Auto_60", "Auto_90", "Auto_120"))
    with(activeColony, boxplot(absdiff~ tech, ylim = c(0, 320), axes=FALSE,
                               ylab = "Absolute error", at = c(1,3,4,5,6,8,9,10,11)))
    axis(side = 2)
    axis(side = 1, at = c(1,4.5, 9.5),
         labels = c("Ground", "UAV manual", "UAV automatic"))
    box()
    colonyN <- activeColony$trueCounts[1]
    text(0.3, 310, paste("Colony ", i, ", N = ", colonyN, sep = ""), cex = 1.5, pos = 4)
}
dev.off()

### Post-hoc analyses for Question One: there appears to be an
 ## interaction effect of tech and colony on absolute error?
m1a <- with(counts, glm(absdiff ~ tech * as.factor(Colony), family = "quasipoisson"))
summary(m1a)
summary(aov(m1a))

counts$tech <- factor(counts$tech, levels = c("Ground_NA", "UAV_manual_30", "UAV_manual_60",
                                           "UAV_manual_90", "UAV_manual_120",
                                           "Auto_30", "Auto_60", "Auto_90", "Auto_120"))
m1r <- with(counts, glmmPQL(absdiff ~ tech, random = ~1 |as.factor(Colony), family = "quasipoisson"))
summary(m1r)
summary(aov(m1r))
TukeyHSD(aov(m1r))

#### plots of the fixed effects in m1r:
## table

##NOTE: I got to here before an unignorable error appeared - Sam

fixedEffects <- as.data.frame(summary(m1r)$tTable)
plotFrame <- data.frame(tech =
                        c("Ground", "GSD 0.82 cm", "GSD 1.64 cm", "GSD 2.47 cm", "GSD 3.29 cm"),
                        est = c(fixedEffects$Value[1],
                                fixedEffects$Value[1]+fixedEffects$Value[2:5]),
                        SE = fixedEffects$Std.Error,
                        P = fixedEffects[5])

## Question Two: What is the difference in variances between ground counts and drone
 #               counts at different heights?

# library(Rcmdr)
test <- with(counts, glmmPQL(Count ~ tech, random = ~1 | as.factor(Colony),
                family = "quasipoisson"))
leveneTest(residuals(test), group = counts$tech)  ## this is the test as pre-registered. It
 # returns a significant effect of tech on count variance. It will need a post-hoc test
 # added to it, though.
counts$splitter <- counts$Colony %in% c(1,2,3,4,6,7)

test2 <- with(counts, lm(Count ~ as.factor(Colony)*tech+splitter))
counts$deviation <- residuals(test2)
m2 <- glm(abs(counts$deviation) ~ counts$tech)
TukeyHSD(aov(m2))

## Upshot: absolute deviation is higher for ground counts than for any drone altitude
 # (max P < 0.001). Also, drone deviation at 30m is less than at 120m (P = 0.005), and
 # drone deviation at 60m is less than 120m (P = 0.029), but deviation at 90m is not
 # significantly less than 120m (P = 0.329).

## Figures for Question Two:

pdf(file = "Q2_variance_v2.pdf", width = 8, height = 8)
par(mfrow = c(3,2), oma = c(0,0,0,0), mar = c(2.5,4.5,1,1))
# for(i in c(2,4,3,10,1,7,8,6,9,5)) { # for all colonies
for(i in c(2,4,3,1,7,6)) { ## in-focus colonies
    activeColony <- droplevels(subset(counts, counts$Colony == i))
    activeColony$tech <- factor(activeColony$tech,
                                levels = c("Ground_NA", "UAV_manual_30", "UAV_manual_60",
                                           "UAV_manual_90", "UAV_manual_120",
                                           "Auto_30", "Auto_60", "Auto_90", "Auto_120"))
    with(activeColony, boxplot(deviation~ tech,
                               ylim = c(-200, 470),
                               axes=FALSE,
                               ylab = "Count deviation from mean",
                               at = c(1,3,4,5,6,8,9,10,11)))
    axis(side = 2)
    axis(side = 1, at = c(1, 4.5, 9.5),
         labels = c("Ground", "UAV manual", "UAV automatic"))
    box()
    colonyN <- activeColony$trueCounts[1]
    text(5.5, 440, paste("Colony ", i, ", N = ", colonyN, sep = ""), cex = 1.5, pos = 2)
}
dev.off()


## Question Three: How are the different counting techniques biased?

counts$diff <- counts$Count - counts$trueCounts

m3 <- with(counts, glmmPQL(diff ~ tech, random = ~1|as.factor(Colony), family = "gaussian"))
summary(m3)

## upshot of the pre-declared analysis: All significant biases are negative biases.
 # Ground-counts are
 # significantly negatively biased (95% CI: -116 - -44), ditto UAV counts at 120m
 # (95% CI: -113 - -60), and UAV counts at 90m (95% CI: -60 - -8). UAV counts at 60m
 # were not significantly negatively biased (95% CI: -35 - 18), and nor were UAV counts
 # at 30m (95% CI: -30 - 23)

## undeclared exploratory analysis: what techniques have significantly different degrees
 # of bias?

TukeyHSD(aov((m3)))

## Figures for Question 3


pdf(file = "Q3_bias.pdf", width = 8, height = 11)
par(mfrow = c(5,2), oma = c(0,0,0,0), mar = c(2.5,4.5,1,1))
for(i in c(2,4,3,10,1,7,8,6,9,5)) {
    activeColony <- droplevels(subset(counts, counts$Colony == i))
    activeColony$tech <- factor(activeColony$tech,
                                levels = c("Ground_NA", "UAV_manual_30", "UAV_manual_60",
                                           "UAV_manual_90", "UAV_manual_120"))
    with(activeColony, boxplot(diff~ tech,
                               ylim = c(-500, 300),
                               axes=FALSE
#                               ylab = "Count deviation from mean")
                               ,type = "n"))
    abline(h = 0, lty = 2)   
    with(activeColony, boxplot(diff~ tech,
                               ylim = c(-500, 300),
                               axes=FALSE,
                               ylab = "Difference from true count", add = TRUE,
                               col = "white"))
    axis(side = 2)
    axis(side = 1, at = c(1,2,3,4,5),
         labels = c("Ground", "UAV 30m", "UAV 60m", "UAV 90m", "UAV 120m"))
    box()
    colonyN <- activeColony$trueCounts[1]
    text(5.5, 250, paste("Colony ", i, ", N = ", colonyN, sep = ""), cex = 1.5, pos = 2)
}
dev.off()


##############  Further plots, tightening to publication versions.

## Simpler plots (not directly based on the model):
boxplot(counts$absdiff~counts$tech, axes = FALSE, ylab = "Absolute error", xlab = "",
     col = "gray80")
axis(side = 2)
axis(side = 1, at = c(1:5),
     labels = c("Ground", "GSD 0.82 cm", "GSD 1.64 cm", "GSD 2.47 cm", "GSD 3.29 cm"))
box()

## Coursed box plots:
inFocus <- subset(counts, counts$Colony %in% c(1,2,3,4,6,7))
inFocus$focus <- c(rep("inFocus", nrow(inFocus)))
allColonies <- counts
allColonies$focus <- c(rep("allColonies", nrow(allColonies)))
focusCounts <- rbind(inFocus, allColonies)

with(focusCounts, boxplot(absdiff ~ focus + tech, axes = FALSE, col = c("grey70", "white"),
                          at = c(1,2,4,5,7,8,10,11,13,14),
                          ylab = "Absolute error"))
axis(side = 2)
axis(side = 1, at = c(seq(1.5, 13.5, 3)),
          labels = c("Ground", "GSD 0.82 cm", "GSD 1.64 cm", "GSD 2.47 cm", "GSD 3.29 cm"))
box()

## Below here are the versions
 # of plots and models that are included in the MS, incorporating semi-automated
 # count data.

## New coursed box plots, incorporating data from semi-automated counts

counts$tech <- factor(counts$tech,
                      levels = c("Ground_NA", "UAV_manual_30", "Auto_30",
                                 "UAV_manual_60", "Auto_60",
                                 "UAV_manual_90", "Auto_90",
                                 "UAV_manual_120", "Auto_120"))
focusCounts$tech <- factor(focusCounts$tech,
                      levels = c("Ground_NA", "UAV_manual_30", "Auto_30",
                                 "UAV_manual_60", "Auto_60",
                                 "UAV_manual_90", "Auto_90",
                                 "UAV_manual_120", "Auto_120"))
inFocus <- subset(focusCounts, focusCounts$focus == "inFocus")
# allColonies <- subset(focusCounts, focusCounts$focus == "allColonies")

## Question One plot (absolute error)
pdf("Q1AbsoluteErrorByFocusAndTech_v2.pdf", width = 5, height = 11, fonts = "Helvetica",
    family = "Helvetica")
par(mar = c(4.5,5.5,1,1))
with(focusCounts, boxplot(absdiff ~ tech,
                          axes = FALSE,
                          col = c("grey70"),
                          at = c(1,4,5,9,10,14,15,19,20),
                          xlab = "Absolute error", horizontal = T, xlim = c(0, 22)))
with(inFocus, boxplot(absdiff ~ tech,
                          axes = FALSE,
                          col = c("white"),
                          at = c(2,6,7,11,12,16,17,21,22),
                          xlab = "Absolute error", horizontal = T, add = TRUE))

axis(side = 1)
axis(side = 2, at = c(1.5, 5.5, 10.5, 15.5, 20.5),
          labels = c("Ground",
                     paste("0.82 cm\n", "(30 m)"),
                     paste("1.64 cm\n", "(60 m)"),
                     paste("2.47 cm\n", "(90 m)"),
                     paste("3.29 cm\n", "(120 m)")),
     las = 2)
box()
dev.off()

## Question Two plot (variability)

pdf("Q2VarianceByFocusAndTech_v2.pdf", width = 5, height = 11, fonts = "Helvetica",
    family = "Helvetica")
par(mar = c(4.5,5.5,1,1))
with(focusCounts, boxplot(deviation ~ tech,
                          axes = FALSE,
                          col = c("grey70"),
                          at = c(1,4,5,9,10,14,15,19,20),
                          xlab = "Count deviation from mean", horizontal = T, xlim = c(0, 22)))
with(inFocus, boxplot(deviation ~ tech,
                          axes = FALSE,
                          col = c("white"),
                          at = c(2,6,7,11,12,16,17,21,22),
                          xlab = "Count deviation from mean", horizontal = T, add = TRUE))

axis(side = 1)
axis(side = 2, at = c(1.5, 5.5, 10.5, 15.5, 20.5),
          labels = c("Ground",
                     paste("0.82 cm\n", "(30 m)"),
                     paste("1.64 cm\n", "(60 m)"),
                     paste("2.47 cm\n", "(90 m)"),
                     paste("3.29 cm\n", "(120 m)")),
     las = 2)
box()
dev.off()

## Question Three plot (bias)

pdf("Q3BiasByFocusAndTech_zeroline_v2.pdf", width = 5, height = 11, fonts = "Helvetica",
    family = "Helvetica")
par(mar = c(4.5,5.5,1,1))
with(focusCounts, boxplot(diff ~ tech,
                          axes = FALSE,
                          col = c("grey70"),
                          at = c(1,4,5,9,10,14,15,19,20),
                          xlab = "Difference from true count", horizontal = T, xlim = c(0, 22),
                          type = "n"))
abline(v = 0, lty =1, col = "grey70")
with(focusCounts, boxplot(diff ~ tech,
                          axes = FALSE,
                          col = c("grey70"),
                          at = c(1,4,5,9,10,14,15,19,20),
                          horizontal = T, xlim = c(0, 22),
                          add = TRUE))
with(inFocus, boxplot(deviation ~ tech,
                          axes = FALSE,
                          col = c("white"),
                          at = c(2,6,7,11,12,16,17,21,22),
                          horizontal = T, add = TRUE))

axis(side = 1)
axis(side = 2, at = c(1.5, 5.5, 10.5, 15.5, 20.5),
          labels = c("Ground",
                     paste("0.82 cm\n", "(30 m)"),
                     paste("1.64 cm\n", "(60 m)"),
                     paste("2.47 cm\n", "(90 m)"),
                     paste("3.29 cm\n", "(120 m)")),
     las = 2)
box()
dev.off()

## New hypothesis tests  ## note that focusCounts is the dataframe containing all
 # colonies.

# question 1a: absolute error, all colonies
m1a <- with(focusCounts,
           glmmPQL(absdiff ~ tech, random = ~1|as.factor(Colony), family = "quasipoisson")
           )
summary(m1a)
summary(aov(m1a))
TukeyHSD(aov(m1a))
         
## question 1b: absolute error, in-focus colonies

m1b <- with(inFocus,
           glmmPQL(absdiff ~ tech, random = ~1|as.factor(Colony), family = "quasipoisson")
           )
summary(m1b)
summary(aov(m1b))
TukeyHSD(aov(m1b))

# question 2a: variance, all colonies - no semi-auto counts

focusCountsM <- focusCounts[c(1:214, 239:589),]
library(Rcmdr)
m2 <- with(focusCountsM,
           glmmPQL(deviation ~ tech, random = ~1|as.factor(Colony), family = "gaussian")
           )

m2a <- glm(abs(residuals(m2)) ~ focusCountsM$tech)
summary(m2a)
TukeyHSD(aov(m2a))

# question 2b: variance, in-focus colonies - no semi-auto counts

inFocusM <- inFocus[1:214,]
library(Rcmdr)
m2 <- with(inFocusM,
           glmmPQL(deviation ~ tech, random = ~1|as.factor(Colony), family = "gaussian")
           )

m2b <- glm(abs(residuals(m2)) ~ inFocusM$tech)
summary(m2b)
TukeyHSD(aov(m2b))

# question 3a: bias, all colonies

m3a <- with(focusCounts,
           glmmPQL(diff ~ tech, random = ~1|as.factor(Colony), family = "gaussian")
           )
summary(m3a)
summary(aov(m3a))
TukeyHSD(aov(m3a))

# question 3b: bias, in-focus colonies

m3b <- with(inFocus,
           glmmPQL(diff ~ tech, random = ~1|as.factor(Colony), family = "gaussian")
           )
summary(m3b)
summary(aov(m3b))
TukeyHSD(aov(m3b))


## calculating the 'now with 95% more accuracy!' figure for the headline.

counts$colonyCT <- as.factor(paste(counts$Colony, counts$tech))

ctRMSE <- data.frame(colonyCT = levels(droplevels(inFocus$colonyCT)),
                     RMSE = c(rep(1, length(levels(droplevels(inFocus$colonyCT))))),
                     sRMSE = c(rep(1, length(levels(droplevels(inFocus$colonyCT)))))
                     ) ## change to counts$colonyCT each time, for all colonies

for (i in levels(ctRMSE$colonyCT)) {
    activeColonyCT <- subset(inFocus, inFocus$colonyCT == i) ## counts here too!
    activeRMSE <- sqrt(sum((activeColonyCT$Count - activeColonyCT$trueCounts)^2)/
                       nrow(activeColonyCT))
    activesRMSE <- activeRMSE / activeColonyCT$trueCounts[1]
    ctRMSE$RMSE[ctRMSE$colonyCT == i] <- activeRMSE
    ctRMSE$sRMSE[ctRMSE$colonyCT == i] <- activesRMSE
}

groundsRMSEs <- ctRMSE[grep("Ground", ctRMSE$colonyCT),]

autosRMSEs <- ctRMSE[grep("Auto", ctRMSE$colonyCT),]
autosRMSEs30 <- ctRMSE[grep("Auto_30", ctRMSE$colonyCT),]
autosRMSEs60 <- ctRMSE[grep("Auto_60", ctRMSE$colonyCT),]
autosRMSEs90 <- ctRMSE[grep("Auto_90", ctRMSE$colonyCT),]
autosRMSEs120 <- ctRMSE[grep("Auto_120", ctRMSE$colonyCT),]

manualsRMSEs <- ctRMSE[grep("manual", ctRMSE$colonyCT),]
manualsRMSEs30 <- ctRMSE[grep("manual_30", ctRMSE$colonyCT),]
manualsRMSEs60 <- ctRMSE[grep("manual_60", ctRMSE$colonyCT),]
manualsRMSEs90 <- ctRMSE[grep("manual_90", ctRMSE$colonyCT),]
manualsRMSEs120 <- ctRMSE[grep("manual_120", ctRMSE$colonyCT),]

comparator <- function(drone, ground) {
    holder <- c(rep(1, nrow(drone)))
    for (i in 1: nrow(drone)) {
        holder[i] <- (1-(drone$sRMSE[i] / ground$sRMSE[i]))*100
    }
    return(mean(holder))
    cat("Drone-counts are", mean(holder), "% more accurate than ground-counts") 
}







## Comparison of drone counts - manual vs semi-auto
## read in data

x <- read.csv("UAVfinal2.csv")
head(x)
as.factor(x$technique)

## test for differences between UAV and auto

x.image <- subset(x, image == "Y")
image.test <- glm(count ~ counter, data=x.image, offset = log(actual), family=poisson())
summary(image.test)

exp(-0.06)#= 0.94 --> 94% as good as humans
exp(0.013)#= 1.01 --> humans typically overcounting by 1%

## test for differences between UAV and auto just in-focus

x.image.focus <- subset(x.image, focus == "Y")
image.focus <- glm(count ~ counter, data=x.image.focus, offset = log(actual), family=poisson())
summary(image.focus)

exp(-0.019)#=  --> 98% as good as humans
exp(0.015)#= 1.01 --> humans typically overcounting by 1% not sure why this is the same as all colonies
