###########################################################################
## OJT 2014 in KOTI
## R script for Sec 4
###########################################################################
loadedOnly <- scan(what="list") 
"extrafontdb" "grid"        "lattice"     "MASS"        "nnet"        "Rttf2pt1"    "tools"   



otherPkgs <- scan(what="list") 
"snow"         "ff"           "bit"          "glmnet"       "Matrix"       "ISLR"         "leaps"        "car"         
"extrafont"    "RColorBrewer"


destdir <- "C:/Users/jang/Documents/OJT/R/packages"
download.packages(otherPkgs, destdir)
download.packages(loadedOnly, destdir)


###########################################################################
###########################################################################
###########################################################################
###########################################################################
reqpkgs <- c("RColorBrewer", "car", "ISLR", "leaps", "glmnet", "parallel", "snow", "Rmpi", "extrafont", "ff")
insedpkgs <- installed.packages()[,1]
inspkg <- reqpkgs[!reqpkgs %in% insedpkgs]
if(length(inspkg) > 0){
	install.packages(inspkg)
}


## ---- page 007 ----------------------------------------------------------------------------------------------
ojtpath <- "~/Dropbox/OJT/"
if(.Platform$OS.type=="windows"){
  user <- shell("echo %username%", intern=TRUE)
  ojtpath <- paste("C:/Users/", user, "/Documents/OJT/", sep="")
}
adv <- read.table(paste(ojtpath, "R/Data/Advertising.data.txt", sep=""), header=TRUE)
attach(adv)
par(mfrow=c(1,3), cex.lab=1.5, mar=c(4, 4.5, 1, 1))
plot(TV, Sales, col=2)
abline( lm(Sales ~ TV), col=4, lwd=2)
plot(Radio, Sales, col=2)
abline( lm(Sales ~ Radio), col=4, lwd=2)
plot(Newspaper, Sales, col=2)
abline( lm(Sales ~ Newspaper), col=4, lwd=2)
detach(adv)


## ---- page 010 ----------------------------------------------------------------------------------------------
attach(adv)
par(mfrow=c(1,1))
plot(TV, Sales, col=2, pch=19)
fit <- lm(Sales ~ TV)
abline(fit, col=4, lwd=2)
segments(TV, Sales, TV, fit$fitted.values, col=4, lty=2)
detach(adv)


## ---- page 011 ----------------------------------------------------------------------------------------------
fit
names(fit)


## ---- page 012 ----------------------------------------------------------------------------------------------
summary(fit)


## ---- page 013 ---------------------------------------------------------------------------------------------- 
names(summary(fit))


## ---- page 014 ---------------------------------------------------------------------------------------------- 
anova(fit)


## ---- page 015 ---------------------------------------------------------------------------------------------- 
attach(adv); par(mar=c(4,4,0,2))
plot(TV, Sales, col=2, pch=19)
fit <- lm(Sales ~ TV)
abline(fit, col=4, lwd=2)
newx <- seq(min(TV)-20, max(TV)+20, , 100)
ci <- predict(fit, data.frame(TV=newx), interval = "confidence")
lines(newx, ci[,2], col=3, lwd=2, lty=2)
lines(newx, ci[,3], col=3, lwd=2, lty=2)
legend("topleft", legend="95% CI", lwd=2, lty=2, col=3)
detach(adv)


## ---- page 016 ---------------------------------------------------------------------------------------------- 
par(mfrow=c(2,2))
plot(fit)


## ---- page 022 ---------------------------------------------------------------------------------------------- 
attach(adv)
fit <- lm(Sales ~ TV + Radio + Newspaper)
fit


## ---- page 023 ---------------------------------------------------------------------------------------------- 
anova(fit)


## ---- page 024 ---------------------------------------------------------------------------------------------- 
summary(fit)


## ---- page 025 ---------------------------------------------------------------------------------------------- 
par(mfrow=c(2,2))
plot(fit)
detach(adv)

###############################################################################################################
## Qualitative predictors
###############################################################################################################
## ---- page 033 ---------------------------------------------------------------------------------------------- 
cred <- read.csv(paste(ojtpath, "R/Data/Credit.csv", sep=""), header=TRUE)
head(cred)


## ---- page 034 ---------------------------------------------------------------------------------------------- 
cred <- read.csv(paste(ojtpath, "R/Data/Credit.csv", sep=""), header=TRUE)
par(mar=c(4, 4.5, 1, 1))
plot(cred[ , c("Balance", "Age", "Cards", "Education", "Income", "Limit", "Rating")], col=4, cex=.6)


## ---- page 036 ---------------------------------------------------------------------------------------------- 
fit <- lm(Balance ~ Gender, data=cred)
summary(fit)


## ---- page 037 ---------------------------------------------------------------------------------------------- 
fit <- lm(Balance ~ Ethnicity, data=cred)
summary(fit)

###############################################################################################################
## Extensions of the linear model
###############################################################################################################
## ---- page 043 ---------------------------------------------------------------------------------------------- 
fit <- lm(Sales ~ TV + Radio + TV:Radio, data=adv) #lm(Sales ~ TV*Radio, data=adv)
summary(fit)


## ---- page 050 ---------------------------------------------------------------------------------------------- 
fit1 <- lm(Balance ~ Income + Student, data=cred) 
summary(fit1)


## ---- page 051 ---------------------------------------------------------------------------------------------- 
fit2 <- lm(Balance ~ Income + Student + Income:Student, data=cred)
summary(fit2)


## ---- page 052 ---------------------------------------------------------------------------------------------- 
newx <- seq(min(cred$Income), max(cred$Income), , 100)
pred1_st <- predict(fit1, newdata=data.frame(Income=newx, Student="Yes"))
pred1_nst <- predict(fit1, newdata=data.frame(Income=newx, Student="No"))
pred2_st <- predict(fit2, newdata=data.frame(Income=newx, Student="Yes"))
pred2_nst <- predict(fit2, newdata=data.frame(Income=newx, Student="No"))
par(mfrow=c(1,2))
plot(cred$Income, cred$Balance, col=cred$Student, pch=19, xlab="Income", ylab="Balance")
lines(newx, pred1_st, col=2, lwd=4); lines(newx, pred1_nst, col=1, lwd=4)
title("no interaction \n between income and student")
legend("topleft", legend=c("student", "non-student"), col=2:1, lwd=2)
plot(cred$Income, cred$Balance, col=cred$Student, pch=19, xlab="Income", ylab="Balance")
lines(newx, pred2_st, col=2, lwd=4); lines(newx, pred2_nst, col=1, lwd=4)
title("with an interaction term \n between income and student")
legend("topleft", legend=c("student", "non-student"), col=2:1, lwd=2)


## ---- page 052 ---------------------------------------------------------------------------------------------- 
vehuse <- read.csv(paste(ojtpath, "R/Data/Vehicleuse_subset.csv", sep=""), 
                   fileEncoding="UTF-8", encoding="UTF-8")
#vehuse <- read.csv(paste(ojtpath, "R/Data/Vehicleuse_subset.csv", sep=""))
colnames(vehuse) <- c("weekday", "time", "sudo", "household", "pur", "trip", "gender",
                      "ln_dis", "ln_time")
vehuse <- vehuse[vehuse$pur !=0 & vehuse$trip !=0, ]
head(vehuse, 4)


#vehuse$weekday <- factor(vehuse$weekday, labels=c("주말", "주중"))
#vehuse$time <- factor(vehuse$time, labels=c("기타시간", "출퇴근시간"))
#vehuse$sudo <- factor(vehuse$sudo, labels=c("비수도권", "수도권"))
#vehuse$gender <- factor(vehuse$gender, labels=c("여성", "남성"))
#vehuse$pur <- factor(ifelse(vehuse$pur <7, 1, 0), labels=c("비일상", "일상"))
#vehuse$household <- cut(as.numeric(vehuse$household), c(0, 2, 4 ,7), 
#                        labels=c("2인이하", "3-4인","5인이상")) 


#vehuse$weekday <- factor(vehuse$weekday, labels=iconv(c("주말", "주중"), "CP949", "UTF-8"))
#vehuse$time <- factor(vehuse$time, labels=iconv(c("기타시간", "출퇴근시간"), "CP949", "UTF-8"))
#vehuse$sudo <- factor(vehuse$sudo, labels=iconv(c("비수도권", "수도권"), "CP949", "UTF-8"))
#vehuse$gender <- factor(vehuse$gender, labels=iconv(c("여성", "남성"), "CP949", "UTF-8"))
#vehuse$pur <- factor(ifelse(vehuse$pur <7, 1, 0), labels=iconv(c("비일상", "일상"), "CP949", "UTF-8"))
#vehuse$household <- cut(as.numeric(vehuse$household), c(0, 2, 4 ,7), 
#                        labels=iconv(c("2인이하", "3~4인","5인이상"), "CP949", "UTF-8")) 


vehuse$weekday <- factor(vehuse$weekday, labels=c("Weekend", "Weekday"))
vehuse$time <- factor(vehuse$time, labels=c("OtherHour", "RushHour"))
vehuse$sudo <- factor(vehuse$sudo, labels=c("NonMetroRegion", "MetroRegion"))
vehuse$gender <- factor(vehuse$gender, labels=c("Female", "Male"))
vehuse$pur <- factor(ifelse(vehuse$pur <7, 1, 0), labels=c("Unusually", "Usually"))
vehuse$household <- cut(as.numeric(vehuse$household), c(0, 2, 4 ,7), 
                        labels=c("LE2", "3-4","GE5")) 



## ---- page 056 ---------------------------------------------------------------------------------------------- 
if(.Platform$OS.type == "windows"){
  library(extrafont)
  font_import(pattern="NanumGothic")

  fonts()
  loadfonts(device="win")
}

attach(vehuse) 
library(RColorBrewer)
par(mar=c(4.5,4.6,4,2), cex.lab=1.5, cex.axis=1.5, family="NanumGothic")
cols <- brewer.pal(8, "Set1")[1:2]
cells <- tapply(ln_time, list(household, weekday), mean)
matplot(exp(cells), type="b", lwd=4, pch=19, col=cols, axes=FALSE, xlab=expression(bold("가구인수")), 
        ylab=expression(bold("통행당 시간")))
axis(1, at=1:3, labels=c("2인이하", "3~4인","5인이상")); axis(2); box()
legend("topright", legend=c("주말", "주중"), col=cols, lwd=4, pch=19, lty=1:2, cex=2)




## ---- page 057 ---------------------------------------------------------------------------------------------- 
bartlett.test(ln_time ~ paste(as.character(household), as.character(weekday), sep="_"), data=vehuse)

# Analysis of Variance 
fit <- lm(ln_time ~ household + weekday + household:weekday, data=vehuse) 

# Anova Tables with  heteroscedasticity-corrected coefficient covariance matrix
library(car)
Anova(fit, type="III", white.adjust="hc3")


## ---- page 058 ---------------------------------------------------------------------------------------------- 
fit <- lm(ln_time ~ (household + sudo + weekday + pur + time + gender)^2, data=vehuse)
# Anova(fit, type="III", white.adjust="hc3")
fit <- lm(ln_time ~ household + sudo + weekday + pur + time + gender + 
                    household:sudo +  household:pur  + household:time + household:gender +
                    sudo:weekday + sudo:pur + 
                    weekday:pur +  weekday:time + 
                    pur:time + pur:gender, data=vehuse)

predmat <- expand.grid(weekday=levels(weekday), time=levels(time), sudo=levels(sudo), 
                       household = levels(household), pur=levels(pur), gender =levels(gender))

predval <- exp(predict(fit, predmat))
out <- cbind(predmat, predval)
head(out)


## ---- page 059 ---------------------------------------------------------------------------------------------- 
par(mfrow=c(2,2), cex.main=2)
for(sudind in c("MetroRegion", "NonMetroRegion")){
  for(weeind in c("Weekday", "Weekend")){
    out1 <- subset(out, sudo==sudind & weekday == weeind)
    out3 <- NULL
    for(genind in c("Male", "Female")){
      for(timind in c("OtherHour", "RushHour")){
        for(purind in c("Unusually", "Usually")){
          out2 <- subset(out1, gender==genind & time==timind & pur==purind)
          out3 <- cbind(out3, out2[,"predval"])
        }
      }
    }
    ksudind <- c("수도권", "비수도권")[c("MetroRegion", "NonMetroRegion") == sudind]
    kweeind <- c("주중", "주말")[c("Weekday", "Weekend") == weeind]
    matplot(out3, type="o", col=rep(cols, each=4), lty=rep(c(1,1,2,2), 2), pch=rep(c(19, 4), 4), lwd=4, 
            axes=FALSE, xlab=expression(bold("가구인수")), ylab=expression(bold("통행당 시간")),
            xlim=c(1,4), ylim=range(predval), main=paste(ksudind, "&", kweeind))
    axis(1, at=1:3, labels=c("2인이하", "3~4인","5인이상")); axis(2); box()
    legend("topright", legend=c("남성", "여성", "기타시간", "출퇴근시간", "비일상", "일상"),
           col=c(cols[1], cols[2], 1, 1, 1, 1), lty=c(1, 1, 1, 2, NA, NA), 
           pch=c(NA, NA, NA, NA, 19, 4), lwd=3)
  }
}
detach(vehuse)


###############################################################################################################
## Non-linear effects of predictors
###############################################################################################################

## ---- page 061 ---------------------------------------------------------------------------------------------- 
library(ISLR)
data(Auto)
par(mfrow=c(1,1)
plot(Auto$horsepower, Auto$mpg, col="grey", xlab="Horsepower", ylab="Miles per gallon")


## ---- page 062 ---------------------------------------------------------------------------------------------- 
fit1 <- lm(mpg ~ horsepower, data=Auto)
fit2 <- lm(mpg ~ horsepower + I(horsepower^2), data=Auto)
summary(fit2)


## ---- page 063 ---------------------------------------------------------------------------------------------- 
plot(Auto$horsepower, Auto$mpg, col="grey", xlab="Horsepower", ylab="Miles per gallon")
newx <- seq(40, 250, , 100)
pred1 <- predict(fit1, data.frame(horsepower=newx))
pred2 <- predict(fit2, data.frame(horsepower=newx))
lines(newx, pred1, col=2, lwd=2)
lines(newx, pred2, col=4, lwd=2)
legend("topright", legend=c("Linear", "Degree 2"), lwd=2, col=c(2,4))


###############################################################################################################
## Qualitative variables
###############################################################################################################

## ---- page 066 ---------------------------------------------------------------------------------------------- 
library(ISLR)
data(Default)
cols <- brewer.pal(8, "Set1")[2:1]
layout(matrix(1:3, ncol=3), width=c(2,1,1))
plot(Default$balance, Default$income, col=cols[as.numeric(Default$default)], 
     pch=c(1,3)[as.numeric(Default$default)], xlab="Balance", ylab="Income")
boxplot(Default$balance ~ Default$default, col=cols, xlab="Default", ylab="Balance")
boxplot(Default$income ~ Default$default, col=cols, xlab="Default", ylab="Income")


###############################################################################################################
## Logistic regression
###############################################################################################################

## ---- page 072 ---------------------------------------------------------------------------------------------- 
fit <- glm (default ~ balance, family = binomial(link = "logit"), data=Default)
summary(fit)



## ---- page 073 ---------------------------------------------------------------------------------------------- 
par(mfrow=c(1,1))
newx <- seq(0, 3000)
pred <- predict(fit, newdata=data.frame(balance=newx), type="response")
cols <- brewer.pal(8, "Set1")[1:2]
plot(Default$balance, as.numeric(Default$default) -1 + runif(length(Default$default), -.05, .05), 
     col=cols[2], pch=4, xlab="Balance", ylab="Probability of Default")
lines(newx, pred, lwd=3, col=cols[1])



## ---- page 075 ----------------------------------------------------------------------------------------------
library(ISLR)
data(Default)
fit <- glm (default ~ balance + income + student, family = binomial(link = "logit"), data=Default)


###############################################################################################################
## Model selection
###############################################################################################################

## ---- page 083 ---------------------------------------------------------------------------------------------- 
cred <- read.csv(paste(ojtpath, "R/Data/Credit.csv", sep=""), header=TRUE)
cred <- cred[,-1]
fit <- lm(Balance ~ ., data=cred)
library(leaps)
fullfit <- regsubsets(Balance ~ ., nvmax=11, data=cred)
summary(fullfit)$rss
summary(fullfit)$rsq


## ---- page 084 ---------------------------------------------------------------------------------------------- 
summary(fullfit)$which[1:8,]


## ---- page 085 ---------------------------------------------------------------------------------------------- 
modset <- NULL
for(i in 1:nrow(summary(fullfit)$which)){
  modset[i] <- paste(colnames(summary(fullfit)$which)[summary(fullfit)$which[i,]][-1], collapse=", ")
}
data.frame(subset=modset, rss=summary(fullfit)$rss, rsq=summary(fullfit)$rsq)[1:6, ]


## ---- page 089 ---------------------------------------------------------------------------------------------- 
forwfit <- regsubsets(Balance ~ ., nvmax=11, method="forward", data=cred)
modset <- NULL
for(i in 1:nrow(summary(forwfit)$which)){
  modset[i] <- paste(colnames(summary(forwfit)$which)[summary(forwfit)$which[i,]][-1], collapse=", ")
}
data.frame(subset=modset, rss=summary(forwfit)$rss, rsq=summary(forwfit)$rsq)[1:4, ]


## ---- page 095 ---------------------------------------------------------------------------------------------- 
adjfit <- regsubsets(Balance ~ ., nvmax=11, method="exhaustive", data=cred)
cols <- brewer.pal(8, "Set1")[c(2,4)]
cp <- summary(adjfit)$cp
bic <- summary(adjfit)$bic
adr2 <- summary(adjfit)$adjr2
par(mfrow=c(1,3), mar=c(4,4.5,4,2))
# Mallow's CP
plot(summary(adjfit)$cp, type="b", col=cols[2], lwd=2, xlab="Number of Predictors",
     ylab=expression(bold(C[p])))
points(cp, col=cols[1], pch=19, cex=2)
points(which.min(cp), min(cp), col=cols[1], pch=4, cex=3)
# BIC
plot(bic, type="b", col=cols[2], lwd=2, xlab="Number of Predictors", ylab="BIC")
points(bic, col=cols[1], pch=19, cex=2)
points(which.min(bic), min(bic), col=cols[1], pch=4, cex=3)
# Adjusted R-square
plot(adr2, type="b", col=cols[2], lwd=2, xlab="Number of Predictors",
     ylab=expression(bold(Adjusted~R^2)))
points(adr2, col=cols[1], pch=19, cex=2)
points(which.max(adr2), max(adr2), col=cols[1], pch=4, cex=3)


## ---- page 098 ---------------------------------------------------------------------------------------------- 
plot(bic, type="b", col=cols[2], lwd=2, xlab="Number of Predictors", ylab="BIC")
points(bic, col=cols[1], pch=19, cex=2)
points(which.min(bic), min(bic), col=cols[1], pch=4, cex=3)

# Validation Set Error
cols <- brewer.pal(8, "Set1")[c(2,4)]
n <- nrow(cred)
set.seed(324)
train <- sample(seq(n), round(3*n/4), replace=FALSE)
regfit <- regsubsets(Balance ~., data=cred[train,], nvmax=11)
 
val.errors <- rep(NA, 11)
x.test <- model.matrix(Balance ~., data=cred[-train, ])# notice the -index!
for(i in 1:11){
  coefi <- coef(regfit, id=i)
  pred <- x.test[ , names(coefi)] %*% coefi
  val.errors[i] = mean((cred$Balance[-train] - pred)^2)
}
plot(sqrt(val.errors), pch=19, type="b", xlab="Number of Predictors",
     ylab="Validation Set Error", col=cols[2])
points(sqrt(val.errors), col=cols[1], pch=19, cex=2)
points(which.min(val.errors), min(sqrt(val.errors)), pch=4, col=cols[1], cex=5)
#points(sqrt(regfit.fwd$rss[-1]/180), col="blue", pch=19, type="b")
#legend("topright", legend=c("Training", "Validation"), col=c("blue","black"), pch=19)


## ---- page 099 ---------------------------------------------------------------------------------------------- 
# Cross-Validation Error
predict.regsubsets <- function(object,newdata,id,...){
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form, newdata)
  coefi <- coef(object, id=id)
  mat[ , names(coefi)]%*%coefi
}
set.seed(11)
folds <- sample(rep(1:10, length=nrow(cred)))
table(folds)
cv.errors <- matrix(NA, 10, 11)
for(k in 1:10){
  best.fit <- regsubsets(Balance ~ ., data=cred[folds!=k,], nvmax=11,method="forward")
  for(i in 1:11){
    pred <- predict(best.fit, cred[folds==k,], id=i)
    cv.errors[k, i] <- mean( (cred$Balance[folds==k] - pred)^2)
  }
}
rmse.cv <- sqrt(apply(cv.errors, 2, mean))
plot(rmse.cv, pch=19, type="b", xlab="Number of Predictors",
     ylab="Cross-Validation Error", col=cols[2])
points(rmse.cv, col=cols[1], pch=19, cex=2)
points(which.min(rmse.cv), min(rmse.cv), pch=4, col=cols[1], cex=5)


## ---- page 106 ---------------------------------------------------------------------------------------------- 
## This is achieved by calling `glmnet` with `alpha=0` (see the helpfile).
## There is also a `cv.glmnet` function which will do the cross-validation for us.
library(glmnet)
x <- model.matrix(Balance ~ . -1, data=cred); y <- cred$Balance
n <- nrow(x)
one <- rep(1, n)
meanx <- drop(one %*% x)/n; x <- scale(x, meanx, FALSE)
normx <- sqrt(drop(one %*% (x^2))); x <- scale(x, FALSE, normx)
fit.ridge <- glmnet(x, y, alpha=0, lambda = exp(seq(log(1e-02), log(1e+05), , 100)) )
ind <- pmatch(c("Income", "Limit", "Rating","StudentYes"), rownames(coefficients(fit.ridge))[-1])
cols <- rep("grey", nrow(coefficients(fit.ridge))-1)
cols[ind] <- c("black", brewer.pal(8, "Set1")[1:3])

par(mfrow=c(1,2))
plot(fit.ridge, xvar="lambda", label=TRUE, col=cols, lwd=2, xlab=expression(lambda),
     ylab="Standardized Coefficients")
legend("topright", legend=c("Income", "Limit", "Rating","Student"), col=cols[ind], lwd=2, bty="n")
plot(fit.ridge, xvar="norm", label=TRUE, col=cols, lwd=2, ylab="Standardized Coefficients",
     xlab=expression(group("||", hat(beta)[lambda]^R, "||")[2] /group("||", hat(beta), "||")[2]))
legend("topleft", legend=c("Income", "Limit", "Rating","Student"), col=cols[ind], lwd=2, bty="n")
# CV plot
cv.ridge <- cv.glmnet(x, y, alpha=0)
plot(cv.ridge)


## ---- page 111 ---------------------------------------------------------------------------------------------- 
## Now we fit a lasso model; for this we use the default `alpha=1`
fit.lasso <- glmnet(x, y)
ind <- pmatch(c("Income", "Limit", "Rating","StudentYes"), rownames(coefficients(fit.lasso))[-1])
cols <- rep("grey", nrow(coefficients(fit.lasso))-1)
cols[ind] <- c("black", brewer.pal(8, "Set1")[1:3])

par(mfrow=c(1,2))
plot(fit.lasso, xvar="lambda", label=TRUE, col=cols, lwd=2, xlab=expression(lambda),
     ylab="Standardized Coefficients")
legend("topright", legend=c("Income", "Limit", "Rating","Student"), col=cols[ind], lwd=2, bty="n")
plot(fit.lasso, xvar="norm", label=TRUE, col=cols, lwd=2, ylab="Standardized Coefficients",
     xlab=expression(group("||", hat(beta)[lambda]^L, "||")[2] /group("||", hat(beta), "||")[2]))
legend("topleft", legend=c("Income", "Limit", "Rating","Student"), col=cols[ind], lwd=2, bty="n")
# CV
cv.lasso <- cv.glmnet(x, y)
plot(cv.lasso)
coef(cv.lasso)


###############################################################################################################
## Low level language with R
###############################################################################################################

## ---- page 123 ---------------------------------------------------------------------------------------------- 
path <-paste(ojtpath, "R/Code/", sep="")
if(.Platform$OS.type=="windows"){
  if(.Platform$r_arch=="i386"){
    dyn.load(paste(path, "fibc32.dll", sep=""))
  } else{
    dyn.load(paste(path, "fibc64.dll", sep=""))
  }
} else{
  dyn.load(paste(path, "fibc.so", sep=""))
}

fibC <- function(n){
    outval <- .C("fibc", as.integer(n), outval=double(1))
    return(outval$outval)  
}
system.time(fibC(10^7))


## ---- page 124 ---------------------------------------------------------------------------------------------- 
fibR <- function(n){
    tmpval2 <- rep(0,2)
    tmpval2[1] <- 0.0
    tmpval2[2] <- 1.0

    if(n < 3){
        outval <- tmpval2[n]
    }else{
        for(i in 3:n){
            outval <- tmpval2[1] + tmpval2[2]
            tmpval2[1] <- tmpval2[2]
            tmpval2[2] <- outval
        }
    }
    return(outval)
}
system.time(fibR(10^7))


###############################################################################################################
## Big data with R
###############################################################################################################

## ---- page 128 ---------------------------------------------------------------------------------------------- 
library(ff)
ffOb <- ff(vmode="double", 1:80)



###############################################################################################################
## Parallels computing
###############################################################################################################

## ---- page 00 ---------------------------------------------------------------------------------------------- 
library(parallel)
workerFunc <- function(n) {
  return(n^2)
}
values <- 1:10000000

## Serial version:
system.time(res <- lapply(values, workerFunc))


## Parallel calculation (mclapply):
system.time(res <- mclapply(values, workerFunc, mc.cores = 2))



## ---- page 134 ---------------------------------------------------------------------------------------------- 
workerFunc <- function(n) { return(n^2) }
values <- 1:100
## Number of workers (R processes) to use:
numWorkers <- 8
## Set up the "cluster"
library(parallel)
cl <- parallel:::makeCluster(numWorkers, type = "PSOCK")
## Parallel calculation (parLapply):
res <- parLapply(cl, values, workerFunc)
## Shut down cluster
stopCluster(cl)
print(unlist(res))


## ---- page 135 ---------------------------------------------------------------------------------------------- 
workerFunc <- function(n) { return(n^2) }
values <- 1:100
library(parallel)
library(snow)
library(Rmpi)
numWorkers <- 2
cl <- makeCluster(numWorkers, type = "MPI")
res <- parLapply(cl, values, workerFunc)
stopCluster(cl)
mpi.exit() # or mpi.quit(), which quits R as well
print(unlist(res))




