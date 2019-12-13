# Final Project MATH 8850

rm(list = ls())
set.seed(1126)

library(MASS)
library(survival)
library(GlobalDeviance)
library(flexsurv)
library(xtable)

# Retreive data
library("Epi")
data(thoro)

n <- dim(thoro)[1]

# Create age and dosage columns
age <- as.numeric((thoro$exitdat - thoro$birthdat)/365.25) # in year
dosage <- c(rep(0,n-1),1)

for (i in 1:n){
  if (thoro$volume[i] == 0){
    dosage[i] = 0
  } else if (thoro$volume[i] <= 20){
    dosage[i] = 1
  } else {
    dosage[i] = 2
  }
}

# Create Cancer Cause of Death column
cancer <- thoro$cause == 2

# Create new thoro2
thoro2 <- data.frame(cbind(thoro$sex, thoro$contrast, thoro$exitstat, thoro$cause,
                cancer, dosage, age))
colnames(thoro2) <- c("Gender", "Treatment", "Exit_Stat", "Cause_Death", "Cancer", "Dosage", "Age")

# Remove cause of death of 15 and 16
thoro2<- thoro2[ -which(thoro2$Cause_Death >14), ]

# Combine Exit_Stat 2 and 3 -> only 2
for (i in 1:nrow(thoro2)){
  if (thoro2$Exit_Stat[i] == 3){
    thoro2$Exit_Stat[i] = 2
  }
}


# Change Treatment to 0:untreated, 1:treated
thoro2$Treatment <- (thoro2$Treatment - 2)*(-1)

# Make thoro3, which includes every subject except those who havn't died
thoro3 <- thoro2[-which(is.na(thoro2$Cancer)),]

####################
# Dats is now clean
####################

# thoro2.1 <- thoro[which(thoro$sex == 1),]
# mean((thoro2.1$injecdat[which(thoro2.1$volume>20)]- thoro2.1$exitdat[which(thoro2.1$volume>20)])/365.25)

####################
# K-M Plot
####################

# K-M Plot for each dosage level x gender WITHOUT SUBSETTING
fit <- survfit(Surv(Age,Exit_Stat == 1) ~ strata(Gender, Dosage), data = thoro2)
#print(fit)
plot(fit, main="K-M Plot for (Dosage Level x Gender) in Thoro Data", xlab="Years", ylab="Proportion Surviving",
     lty=c(1,1,1,2,2,2), col=rep(c(2,3,4),2))
legend('bottomleft', c("Male", "Female",
                       "Untreated", "<= 20 Dosage", "> 20 Dosage"), 
       lty=c(1,2,1,1,1), col=c(1,1,2,3,4))

# Cut after age 40 for thoro2 and thoro3
thoro2<- thoro2[ -which(thoro2$Age <= 40), ]
thoro2$Age <- thoro2$Age - 40

thoro3<- thoro3[ -which(thoro3$Age <= 40), ]

## Check that there are no zeros
# min(thoro2$Age)
# thoro2$Age[which(thoro2$Age ==0)]

# K-M Plot for each dosage level x gender WITH SUBSETTING
fit <- survfit(Surv(Age,Exit_Stat == 1) ~ strata(Gender, Dosage), data = thoro2)
#summary(fit)
plot(fit, main="K-M Plot for (Dosage Level x Gender) in Thoro Data \n Conditioned for Age > 40", 
     xlab="Years After 40", ylab="Proportion Surviving",
     lty=c(1,1,1,2,2,2), col=rep(c(2,3,4),2))
legend('topright', c("Male", "Female",
                       "Untreated", "<= 20 Dosage", "> 20 Dosage"), 
       lty=c(1,2,1,1,1), col=c(1,1,2,3,4))

#####################
# Logistic Regression
#####################

# Using thoro3
# Treatement: 0 untreated, 1 treated
# Cancer: 0 not diagnosed, 1 cancer related cause of death

# table(thoro3$Treatment,thoro3$Cancer)

logit <- glm(Cancer ~ Treatment, data = thoro3, family = "binomial")
summ <- summary(logit)
xtable(summ)
intercept <- summ$coefficients[1,1]
Trt_Effect <- summ$coefficients[2,1]

plogis(intercept)
plogis(intercept + Trt_Effect)


######################
# AFT
######################

expAFT <- flexsurvreg(Surv(Age,Exit_Stat == 1) ~ as.factor(Gender) * as.factor(Dosage), data = thoro2, dist = "exponential")
#coef(expAFT)

weiAFT <- flexsurvreg(Surv(Age,Exit_Stat == 1) ~ as.factor(Gender) * as.factor(Dosage), data = thoro2, dist = "weibull")
coef(weiAFT)

lnoAFT <- flexsurvreg(Surv(Age,Exit_Stat == 1) ~ as.factor(Gender) * as.factor(Dosage), data = thoro2, dist = "lnorm")
#coef(lnoAFT)

gamAFT <- flexsurvreg(Surv(Age,Exit_Stat == 1) ~ as.factor(Gender) * as.factor(Dosage), data = thoro2, dist = "gamma")
#coef(gamAFT)

# Weibull has lowest AIC:
AIC <- matrix(c(expAFT$AIC, weiAFT$AIC, lnoAFT$AIC, gamAFT$AIC), nrow = 1)
colnames(AIC) <- c("exponential", "weibull", "log-normal","gamma")
xtable(AIC)

#################
# wei curves
#################

thoro_men <- thoro2[-which(thoro2$Gender==2),]
fit2 <- survfit(Surv(Age,Exit_Stat == 1) ~ strata(Dosage), data = thoro_men)
plot(fit2, main="K-M Plot for Dosage Level of Men in Thoro Data \n With Weibull AFT Curves Overlapped",
     xlab="Years", ylab="Proportion Surviving", lty=c(1,1,1), col=c(2,3,4))
legend('bottomleft', c("Untreated", "<= 20 Dosage", "> 20 Dosage"), 
       lty=c(1,1,1), col=c(2,3,4))

# Male Untreated
curve(1-pweibull(x, scale = exp(coef(weiAFT)[2]),shape = exp(coef(weiAFT)[1])),0,70, col=1, add = T) 
# Male Treated with dosage < 20
curve(1-pweibull(x, scale = exp(coef(weiAFT)[2] + coef(weiAFT)[4]),shape = exp(coef(weiAFT)[1])),0,70, add = T) 
# Male Treated with dosage > 20
curve(1-pweibull(x, scale = exp(coef(weiAFT)[2] + coef(weiAFT)[5]),shape = exp(coef(weiAFT)[1])),0,70, add = T)

thoro_women <- thoro2[-which(thoro2$Gender==1),]
fit2 <- survfit(Surv(Age,Exit_Stat == 1) ~ strata(Dosage), data = thoro_women)
plot(fit2, main="K-M Plot for Dosage Level of Women in Thoro Data \n With Weibull AFT Curves Overlapped",
     xlab="Years", ylab="Proportion Surviving", lty=c(1,1,1), col=c(2,3,4))
legend('bottomleft', c("Untreated", "<= 20 Dosage", "> 20 Dosage"), 
       lty=c(1,1,1), col=c(2,3,4))

# Female Untreated
curve(1-pweibull(x, scale = exp(coef(weiAFT)[2] + coef(weiAFT)[3]),shape = exp(coef(weiAFT)[1])),0,70, add = T) 
# Female Treated with dosage < 20
curve(1-pweibull(x, scale = exp(coef(weiAFT)[2] + coef(weiAFT)[3] + coef(weiAFT)[4] + coef(weiAFT)[6]),shape = exp(coef(weiAFT)[1])),0,70, add = T) 
# Female Treated with dosage > 20
curve(1-pweibull(x, scale = exp(coef(weiAFT)[2] + coef(weiAFT)[3] + coef(weiAFT)[5] + coef(weiAFT)[7]),shape = exp(coef(weiAFT)[1])),0,70, add = T)

#############################################
# Obtain Means and CI's for Expected Survival
#############################################

coef <- coef(weiAFT)

# Means for Weibull
male.dosage0 <- integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]),
                                                     scale = exp(coef[2]))}, 0, Inf)$value
male.dosage1 <- integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]),
                                                     scale = exp(coef[2]+coef[4]))}, 0, Inf)$value
male.dosage2 <- integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]),
                                                     scale = exp(coef[2]+coef[5]))}, 0, Inf)$value
female.dosage0 <- integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]),
                                                       scale = exp(coef[2]+coef[3]))}, 0, Inf)$value
female.dosage1 <- integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]),
                                                       scale = exp(coef[2]+coef[3]+coef[4]+coef[6]))}, 0, Inf)$value
female.dosage2 <- integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]),
                                                       scale = exp(coef[2]+coef[3]+coef[5]+coef[7]))}, 0, Inf)$value

# mean matrix
mean.weibull <- matrix(rbind(c(male.dosage2, male.dosage1, male.dosage0),
                             c(female.dosage2, female.dosage1, female.dosage0)),
                       ncol = 3,nrow = 2)
colnames(mean.weibull) <- c("dosage2","dosage1","no trt")
rownames(mean.weibull) <- c("Male", "Female")
mean.weibull <- mean.weibull + 40
mean.weibull
xtable(mean.weibull)

# Set up for getting CI's
g0 <- c(0,1,0,0,0,0,0)
g1 <- c(0,1,0,1,0,0,0)
g2 <- c(0,1,0,0,1,0,0)
g3 <- c(0,1,1,0,0,0,0)
g4 <- c(0,1,1,1,0,1,0)
g5 <- c(0,1,1,0,1,0,1)

# Obtaining SE for scale
cov <- weiAFT$cov

sd0 <- sqrt(g0 %*% cov %*% g0)
sd1 <- sqrt(g1 %*% cov %*% g1)
sd2 <- sqrt(g2 %*% cov %*% g2)
sd3 <- sqrt(g3 %*% cov %*% g3)
sd4 <- sqrt(g4 %*% cov %*% g4)
sd5 <- sqrt(g5 %*% cov %*% g5)

# CI's
CI.m.dosage0 <- c(male.dosage0 /exp(1.96*sd0), male.dosage0 *exp(1.96*sd0)) + 40
CI.m.dosage1 <- c(male.dosage1/exp(1.96*sd1), male.dosage1*exp(1.96*sd1)) + 40
CI.m.dosage2 <- c(male.dosage2/exp(1.96*sd2), male.dosage2*exp(1.96*sd2)) + 40
CI.f.dosage0 <- c(female.dosage0/exp(1.96*sd3), female.dosage0*exp(1.96*sd3)) + 40
CI.f.dosage1 <- c(female.dosage1/exp(1.96*sd3), female.dosage1*exp(1.96*sd3)) + 40
CI.f.dosage2 <- c(female.dosage2/exp(1.96*sd3), female.dosage2*exp(1.96*sd3)) + 40

CI.weibull <- matrix(rbind(c(CI.m.dosage2, CI.m.dosage1, CI.m.dosage0),
                             c(CI.f.dosage2, CI.f.dosage1, CI.f.dosage0)),
                       ncol = 6,nrow = 2)
colnames(CI.weibull) <- c("dosage >= 20 LCL", "dosage >= 20 UCL", "dosage < 20 LCL", "dosage < 20 UCL", "no trt LCL", "no trt UCL")
rownames(CI.weibull) <- c("Male", "Female")
CI.weibull
xtable(CI.weibull)

######################
# Cox (not used)
######################

## doing a cox model fit for our stratified gender and dosage
# coxFit <- coxph(survObj ~ as.factor(thoro2$Gender) * as.factor(thoro2$Dosage), method = "efron")
# summary(coxFit)

## Use cox.zph to analyze if assumptions tenable
# cox.zph(coxFit)

## High p-value says our proportional hazards assumption from the cox model is not violated 
# plot(fit, log="xy", fun="cumhaz", col=1:3)
