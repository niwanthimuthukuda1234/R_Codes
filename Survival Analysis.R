library(survival)
library(My.stepwise)
library(MASS)
library(survivalMPL)
library(survminer)
library(flexsurv)
library(SurvRegCensCov)
library(reliaR)
getwd()
setwd("F:/Analysis")

#Importing dataset
data=read.csv("Data.csv",header = TRUE)
attach(data)
print(data)
Surv(data$SURVIVAL,data$STATUS==1)

#Define Variables
x=cbind(GENDER,BLOOD_GP,RH_FACTOR,KTB,ferritin,hb,AD,ANOT)

#Descriptive Statistics
summary(SURVIVAL)
summary(hb)
summary(ferritin)
summary(Nblod)

#Non-paraetric Analysis
#Kaplan-Meier method
#Overall Kaplan-Meier Survival Curve
kmsurvival1=survfit(Surv(SURVIVAL,STATUS)~1)
summary(kmsurvival1)
plot(kmsurvival1,xlab = " Survival Time ( Year )",ylab ="Survival Probability")
title("Overall Kaplan Meier Survival Curve")
print(kmsurvival1)

#Kaplan-Meier Survival Curve by Gender
kmsurvival2=survfit(Surv(SURVIVAL,STATUS)~GENDER)
summary(kmsurvival2)
print(kmsurvival2)
plot(kmsurvival2,col = c("blue","red"))
plot(kmsurvival2,conf.int = "none",col = c("blue","red"),xlab = "Survival Time ( Year )",ylab = "Survival Probability")
legend("topright",c("Male","Female"),col = c("blue","red"),lty = 1)
title("Kaplan Meier Survival curves by Gender")

#Kaplan-Meier Survival Curve by Blood Group
kmsurvival3=survfit(Surv(SURVIVAL,STATUS)~BLOOD_GP)
summary(kmsurvival3)
print(kmsurvival3)
plot(kmsurvival3,col = c("red","blue","green","black"))
plot(kmsurvival3,conf.int = "none",col = c("red","blue","green","black"),xlab = "Survival Time ( Year )",ylab = "Survival Probability")
legend("bottomleft",c("A","B","AB","O"),col = c("red","blue","green","black"),lty = 1)
title("Kaplan Meier Survival curves by Blood Group")

#Kaplan-Meier Survival Curve by Rh Factor
kmsurvival4=survfit(Surv(SURVIVAL,STATUS)~RH_FACTOR)
summary(kmsurvival4)
print(kmsurvival4)
plot(kmsurvival4,col = c("red","blue"))
plot(kmsurvival4,conf.int = "none",col = c("red","blue"),xlab = "Survival Time ( Year )",ylab = "Survival Probability")
legend("bottomleft",c("Positive","Negative"),col = c("red","blue"),lty = 1)
title("Kaplan Meier Survival curves by Rh Factor")

#Kaplan-Meier Survival Curve by Kind of Blood Transfused
kmsurvival5=survfit(Surv(SURVIVAL,STATUS)~KTB)
summary(kmsurvival5)
print(kmsurvival5)
plot(kmsurvival5,col = c("red","blue"))
plot(kmsurvival5,conf.int = "none",col = c("red","blue"),xlab = "Survival Time ( Year )",ylab = "Survival Probability")
legend("bottomleft",c("Filtrated","Washed"),col = c("red","blue"),lty = 1)
title("Kaplan Meier Survival curves by kind of blood transfused")

#Kaplan-Meier Survival Curve by annual number of transfusions
kmsurvival8=survfit(Surv(SURVIVAL,STATUS)~ANOT)
summary(kmsurvival8)
print(kmsurvival8)
plot(kmsurvival8,col = c("red","blue"))
plot(kmsurvival8,conf.int = "none",col = c("red","blue"),xlab = "Survival Time ( Year )",ylab = "Survival Probability")
legend("bottomleft",c("<= 12","> 12"),col = c("red","blue"),lty = 1)
title("Kaplan Meier Survival curves by Annual number of transfusions")

#Kaplan-Meier Survival Curve by Hb level
kmsurvival6=survfit(Surv(SURVIVAL,STATUS)~HB)
summary(kmsurvival6)
print(kmsurvival6)
plot(kmsurvival6,col = c("red","blue"))
plot(kmsurvival6,conf.int = "none",col = c("red","blue"),xlab = "Survival Time ( Year )",ylab = "Survival Probability")
legend("bottomleft",c("<= 9","> 9"),col = c("red","blue"),lty = 1)
title("Kaplan Meier Survival curves by Hb level")

#Kaplan-Meier Survival Curve by Accompanied Diseases
kmsurvival7=survfit(Surv(SURVIVAL,STATUS)~AD)
summary(kmsurvival7)
print(kmsurvival7)
plot(kmsurvival7,col = c("red","blue"))
plot(kmsurvival7,conf.int = "none",col = c("red","blue","green","black","purple","orange"),xlab = "Survival Time ( Year )",ylab = "Survival Probability")
legend("bottomleft",c("None","Heart disease","Kidney disease","Hepatitis","Diabetes","Liver disorder"),col = c("red","blue","green","black","purple","orange"),lty = 1)
title("Kaplan Meier Survival curves by Accompanied diseases")


#Log-rank test
#for Gender
survdiff(Surv(SURVIVAL,STATUS)~GENDER)

#for Blood group
survdiff(Surv(SURVIVAL,STATUS)~BLOOD_GP)

#for Rh factor
survdiff(Surv(SURVIVAL,STATUS)~RH_FACTOR)

#for Kind of transfused blood
survdiff(Surv(SURVIVAL,STATUS)~KTB)

#for Annual number of transfusions
survdiff(Surv(SURVIVAL,STATUS)~ANOT)

#for Hb levels
survdiff(Surv(SURVIVAL,STATUS)~HB)

#for Accompanied diseases
survdiff(Surv(SURVIVAL,STATUS)~AD)


#Semi-parametric analysis
#univariate cox regression

covariates=c("AGE","GENDER","BLOOD_GP","RH_FACTOR","KTB","hb","Nblod","ferritin","AD")
univ_formulas=sapply(covariates,function(x)as.formula(paste('Surv(SURVIVAL,STATUS)~',x)))
univ_models=lapply(univ_formulas,function(x){coxph(x,data = data)})
#Extract data
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

#multivariate Cox regression analysis

mydata=na.omit(data)
dim(data)
dim(mydata)
head(mydata)
Surv(mydata$SURVIVAL,mydata$STATUS==1)
myvariablelist1=c("AGE","GENDER","BLOOD_GP","RH_FACTOR","KTB","Nblod","hb","ferritin","AD")
FinalCoxPHmodel=My.stepwise.coxph(Time ="SURVIVAL",Status = "STATUS",variable.list = myvariablelist1,in.variable = c("KTB","hb","Nblod","ferritin","AD"),data = data,sle = 0.15,sls = 0.15)

#baseline hazard function
CoxpHmodel=coxph(Surv(SURVIVAL,STATUS)~BLOOD_GP+hb+AD+Nblod+KTB+ferritin,data = mydata)
CoxpHmodel
basehaz(CoxpHmodel)
AIC(CoxpHmodel)

#model diagnostics
#Graphical method ( for categorical covariates )
#for kind of blood transfused
plot(kmsurvival5,fun = function(s)-log(-log(s)),xlab = "Survival Time ( Year )",ylab = "-log(-log(Survival Probability))",col = c("red","blue"))
legend("topright",c("Filtrated","Washed"),col = c("red","blue"),lty = 1)

#for accompanied diseases
plot(kmsurvival7,fun = function(s)-log(-log(s)),xlab = "Survival Time ( Year )",ylab = "-log(-log(Survival Probability))",col = c("red","blue","green","black","purple","orange"))
legend("topright",c("None","Heart disease","Kidney disease","Hepatitis","Diabetes","Liver disorder"),col = c("red","blue","green","black","purple","orange"),lty = 1)

#for blood group
plot(kmsurvival3,fun = function(s)-log(-log(s)),xlab = "Survival Time ( Year )",ylab = "-log(-log(Survival Probability))",col = c("red","blue","green","black"))
legend("topright",c("A","B","AB","O"),col = c("red","blue","green","black"),lty = 1)

#if the lines are not parallel, the assumption is ivalid
#for continuous covariates
Schoenfeld=cox.zph(CoxpHmodel,transform="km")
Schoenfeld
plot(Schoenfeld)
##null hypothesis is the PH assumption is valid

#Parametric analysis
#Should have to eliminate the observations with zero survival times
#AFT models
data1=read.csv("Data1.csv",header = TRUE)
attach(data1)
print(data1)
Surv(data1$SURVIVAL,data1$STATUS==1)
x=c("GENDER","BLOOD_GP","RH_FACTOR","KTB","ferritin","hb","AD","Nblod")

#full model
Weibullmdl1=survreg(Surv(SURVIVAL,STATUS)~KTB+hb+Nblod+ferritin+AD+BLOOD_GP+GENDER+RH_FACTOR,data = data1,dist = "weibull")
summary(Weibullmdl1)
AIC(Weibullmdl1)

#removed gender
Weibullmdl2=survreg(Surv(SURVIVAL,STATUS)~KTB+hb+Nblod+ferritin+AD+BLOOD_GP+RH_FACTOR,data = data1,dist = "weibull")
summary(Weibullmdl2)
AIC(Weibullmdl2)

#Removed Rh factor
Weibullmdl3=survreg(Surv(SURVIVAL,STATUS)~KTB+hb+Nblod+ferritin+AD+BLOOD_GP,data = data1,dist = "weibull")
summary(Weibullmdl3)
AIC(Weibullmdl3)

#Removed blood group
Weibullmdl4=survreg(Surv(SURVIVAL,STATUS)~KTB+hb+Nblod+ferritin+AD,data = data1,dist = "weibull")
summary(Weibullmdl4)
AIC(Weibullmdl4)

#Removed ferritin level
Weibullmdl5=survreg(Surv(SURVIVAL,STATUS)~KTB+hb+Nblod+AD,data = data1,dist = "weibull")
summary(Weibullmdl5)
AIC(Weibullmdl5)

#other AFT models for model 5
#exponential
expmdl5=survreg(Surv(SURVIVAL,STATUS)~KTB+hb+Nblod+AD,data = data1,dist = "exponential")
summary(expmdl5)
AIC(expmdl5)

#loglogistic 
llmdl5=survreg(Surv(SURVIVAL,STATUS)~KTB+hb+Nblod+AD,data = data1,dist = "loglogistic")
summary(llmdl5)
AIC(llmdl5)

#lognormal
lnmdl5=survreg(Surv(SURVIVAL,STATUS)~KTB+hb+Nblod+AD,data = data1,dist = "lognormal")
summary(lnmdl5)
AIC(lnmdl5)

# model diagnostics of the fitted weibull model ( by deviance residuals)
resid.deviance = residuals(Weibullmdl5, type="deviance")
par(mfrow=c(2,2))
boxplot(resid.deviance ~ HB,xlab="Hemoglobin level")
title("Deviance residuals vs hemoglobin level")
boxplot(resid.deviance ~ KTB,xlab = "Kind of transfused blood")
title("Deviance residuals vs Kind of transfused blood")
boxplot(resid.deviance ~ AD,xlab = "accompanied diseases")
title("Deviance residuals vs accompanied diseases")
boxplot(resid.deviance ~ ANOT,xlab = "annual number of transfusions")
title("Deviance residuals vs annual number of transfusions")







