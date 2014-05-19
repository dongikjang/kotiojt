###########################################################################
## OJT 2014 in KOTI
## R script for Sec 2
###########################################################################

otherPkgs <- scan(what="list")
"XML"             "colorspace"      "RGtk2"           "rattle"          "topicmodels"     "wordcloud"      
"RColorBrewer"    "Rcpp"            "SnowballC"       "tm"              "twitteR"         "rjson"          
"ROAuth"          "digest"          "RCurl"           "bitops"          "arulesViz"       "arules"         
"Matrix"          "amap"            "fpc"             "flexmix"         "mclust"          "cluster"        
"animation"       "sqldf"           "RSQLite.extfuns" "RSQLite"         "DBI"             "gsubfn"         
"proto"           "doBy"            "survival"        "car"             "gdata"           "RODBC"          
"xlsx"            "xlsxjars"        "rJava"           "MASS"            "DAAG"            "lattice"        
"party"           "modeltools"      "strucchange"     "sandwich"        "zoo"             "plyr"           
"data.table" 


loadedOnly <- scan(what="list")
"chron"         "coin"          "gclus"         "gtools"        "igraph"        "latticeExtra"  "lme4"         
"minqa"         "mvtnorm"       "nlme"          "nnet"          "parallel"      "RcppEigen"     "reshape2"     
"scatterplot3d" "seriation"     "slam"          "stringr"       "tools"         "TSP"           "vcd"   


user <- shell("echo %username%", intern=TRUE)
ojtpath <- paste("C:/Users/", user, "/Documents/OJT/", sep="")

destdir <- paste(ojtpath, "R/packages", sep="")
download.packages(otherPkgs, destdir)
download.packages(loadedOnly, destdir)

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################


## Install required packages
reqpkgs <- c("xlsx", "RODBC", "gdata", "DAAG", "car", "doBy", "sqldf", "tcltk", "plyr",
             "doMC", "fpc", "arules", "arulesViz", "twitteR", "rattle", "animation", 
             "amap", "tm", "SnowballC", "wordcloud", "topicmodels",
             "XML", "RColorBrewer", "party", "data.table", "MASS")
insedpkgs <- installed.packages()[,1]
inspkg <- reqpkgs[!reqpkgs %in% insedpkgs]
if(length(inspkg) > 0){
  install.packages(inspkg)
}


## ---- page 010 ----------------------------------------------------------------------------------------------
install.packages("fields")
update.packages() #to update all your package
update.packages("fields") #to update `fields' package


## ---- page 016 ----------------------------------------------------------------------------------------------
2 + 2
2^2
(1 - 2) * 3
1 - 2 * 3


## ---- page 017 ----------------------------------------------------------------------------------------------
sqrt(2)
sin(pi)
exp(1)
log(10)
log(10,10)


## ---- page 018 ----------------------------------------------------------------------------------------------
squareroot(2)
sqrt(-2)

sqrt 2


## ------------------------------------------------------------------------
x <- 2
x + 3
pi
e^2


## ------------------------------------------------------------------------
whales <- c(74, 122, 235, 111, 292, 111, 211, 133, 156, 79)
whales


x <- c(74, 122, 235, 111, 292)
y <- c(111, 211, 133, 156, 79)
c(x,y)


## ---- page 021 ----------------------------------------------------------------------------------------------
simpsons <- c("Homer",'Marge',"Bart","Lisa","Maggie")


## ---- page 022 ----------------------------------------------------------------------------------------------
names(simpsons) <- c("dad","mom","son","daughter 1","daughter 2")
names(simpsons)
simpsons


## ---- page 023 ----------------------------------------------------------------------------------------------
x <- c(5,12,13) 
x
length(x)
mode(x)


## ---- page 024 ----------------------------------------------------------------------------------------------
y <- "abc"
y
length(y)
mode(y)
z <- c("abc","29 88") 
length(z)
mode(z)


## ---- page 025 ----------------------------------------------------------------------------------------------
m <- rbind(c(1,4), c(2,2))  ## m <- matrix(c(1,4,2,2), 2, 2)
m
m %*% c(1,1)


## ---- page 026 ----------------------------------------------------------------------------------------------
x <- list(u=2, v="abc") 
x
names(x)
x$u 


## ---- page 027 ----------------------------------------------------------------------------------------------
d <- data.frame(list(kids=c("Jack","Jill"),ages=c(12,10)))
d
d$ages
d[,2]
d[,"ages"]


## ---- page 028 ----------------------------------------------------------------------------------------------
# Nile
hn <- hist(Nile, plot=FALSE)
str(hn)


## ---- page 029 ----------------------------------------------------------------------------------------------
plot(hn)
getS3method("plot", "histogram")
# graphics:::plot.histogram


## ---- page 030 ----------------------------------------------------------------------------------------------
whales
sum(whales)
length(whales)
sum(whales)/length(whales)
mean(whales)


## ---- page 031 ----------------------------------------------------------------------------------------------
whales
sort(whales)
min(whales)
max(whales)


## ---- page 032 ----------------------------------------------------------------------------------------------
range(whales)
diff(whales)
cumsum(whales)


## ---- page 033 ----------------------------------------------------------------------------------------------
a <- c(1,2,3,4)
b <- c(5,6,7,8)
a + b
a - b
a - mean(a)


## ---- page 037 ----------------------------------------------------------------------------------------------
1:10
rev(1:10)
10:1
a=1; h=4; n=5
a + h*( 0:(n-1) )


## ---- page 038 ----------------------------------------------------------------------------------------------
seq(1, 9, by=2)
seq(1, 9, length=5)


## ---- page 039 ----------------------------------------------------------------------------------------------
rep(1, 10)
rep(1:3, 3)
rep(1:3, each=3)
rep(c("long", "short"), c(1, 2))


## ---- page 040 ----------------------------------------------------------------------------------------------
ebay <- c(88.8, 88.3, 90.2, 93.5, 95.2, 94.7, 99.2, 99.4, 101.6)
length(ebay)
ebay[1]
ebay[9]
ebay[length(ebay)]    # in case length isn't known


## ---- page 041 ----------------------------------------------------------------------------------------------
ebay
ebay[1:4]
ebay[c(1,5,9)]


ebay[-1]		# all but the first
ebay[-(1:4)]		# all but the 1st - 4th


## ---- page 042 ----------------------------------------------------------------------------------------------
x <- 1:3
names(x) <- c("one","two","three")		# set the names
x["one"]


## ---- page 043 ----------------------------------------------------------------------------------------------
ebay
ebay[1] <- 88.0
ebay
ebay[10:13] <- c(97.0, 99.3, 102.0, 101.8)
ebay

ebay[1:5] <- c(10.1,20.2, 30.1)


## ---- page 044 ----------------------------------------------------------------------------------------------
ebay
ebay>100
ebay[ebay>100]   # values bigger than 100


## ---- page 045 ----------------------------------------------------------------------------------------------
which(ebay>100)   #which indices
ebay[c(9,12,13)]   # directly
sum(ebay>100)   # number bigger than 100
sum(ebay>100)/length(ebay)   # proportion bigger


## ---- page 046 ----------------------------------------------------------------------------------------------
x <- 1:5   
x<5   # is x less than 5   
x>1   # is x more than 1   
x>1 & x<5   # is x bigger than 1 and less than 5   
x>1 && x<5    # First one is false   



## ---- page 047 ----------------------------------------------------------------------------------------------
x>1 | x<5   # is x bigger than 1 or less than 5  
x>1||x<5   # First one true   
x==3    # is x equal to 3
x!=3    # is x not equal to 3
!x==3    # NOT (X EQUAL TO 3)


## ---- page 048 ----------------------------------------------------------------------------------------------
x==c(2,4)   
x %in% c(2,4)   


## ---- page 049 ----------------------------------------------------------------------------------------------
shuttle <- c(0, 1, 0, NA, 0, 0)   
shuttle>0    # note NA is answer   
shuttle==NA    # doesn't work!   
is.na(shuttle)   
mean(shuttle)    # can't add to get the mean   
mean(shuttle, na.rm=TRUE)    # na.rm means remove NA   




## ---- page 054 ----------------------------------------------------------------------------------------------
library(MASS)
names(geyser)
head(geyser$waiting)


## ---- page 056 ----------------------------------------------------------------------------------------------
#Insurance <- read.csv("Insurance.csv", header=TRUE)
library(xlsx)
#read.xlsx("Insurance.xlsx", sheetIndex=1)


# library(RODBC)
# con <- dbConnect(driver, user, password, host, dbname)
# Insurance <- dbSendQuery(con, "SELECT * FROM claims")


# con <- url("http://labs.datacenter.com/test.txt")
# Insurace <- read.csv(con, header=TRUE)


# load("Insurace.RData")


## ---- page 057 ----------------------------------------------------------------------------------------------
tmp <- scan()
74 122 235 11 292 11 211 133 156 79



## ---- page 058 ----------------------------------------------------------------------------------------------
ojtpath <- "~/Dropbox/OJT/"
if(.Platform$OS.type=="windows"){
  user <- shell("echo %username%", intern=TRUE)
  ojtpath <- paste("C:/Users/", user, "/Documents/OJT/", sep="")
}
#tmp <- scan(file=paste(ojtpath, "R/Data/temp.txt", sep=""))



## ---- page 059 ----------------------------------------------------------------------------------------------
#read.table(file=file.choose())


## ---- page 060 ----------------------------------------------------------------------------------------------
read.table("clipboard", sep="\t", header=FALSE) # only for Windows
# read.table(pipe('pbpaste'), header=FASET) # for Mac OS X


library(gdata)
#read.xls(file.choose(), sheet=1, header=FALSE)  # Perl-based solution


library(xlsx)
#read.xlsx(file.choose(), sheetIndex=1, header=FALSE)  # Java-based solution


## ---- page 061 ----------------------------------------------------------------------------------------------
demo(graphics)


## ---- page 063 ----------------------------------------------------------------------------------------------
plot((0:20)*pi/10, sin((0:20)*pi/10))
plot((1:30)*0.92, sin((1:30)*0.92))



## ---- page 065 ----------------------------------------------------------------------------------------------
library(DAAG)
attach(elasticband) # R now knows where to find distance & stretch
par(mfrow=c(1,3), mar=c(4.7,4,.3,2), cex.lab=1.5, cex.axis=1.5)
plot(distance ~ stretch, pch=19)
plot(ACT ~ year, data=austpop, type="l", lwd=2)
plot(ACT ~ year, data=austpop, type="b", lwd=2)
detach(elasticband) 


## ---- page 066 ----------------------------------------------------------------------------------------------
data(austpop)
attach(austpop)
plot(year, ACT)
fit <- spline(year, ACT)	# Fit smooth curve through points
lines(fit$x, fit$y, type="l") 	
text(locator(1), label="Qubic spline fit", cex=2, col=4)
detach(austpop) 	


## ---- page 068 ----------------------------------------------------------------------------------------------
plot(hills)  # Has the same effect as pairs(hills)


## ---- page 070 ----------------------------------------------------------------------------------------------
mat <- matrix(1:4, 2, 2)
mat
layout(mat)
layout.show(4)  #par(mfcol=c(2,2))


## ---- page 071 ----------------------------------------------------------------------------------------------
layout(matrix(1:6, 3, 2))
layout.show(6)	#par(mfcol=c(3,2))


## ---- page 072 ----------------------------------------------------------------------------------------------
layout(matrix(1:6, 2, 3, byrow=T))
layout.show(6)	#par(mfrow=c(2,3))

m <- matrix(c(1:3,3), 2, 2)
layout(m)
layout.show(3)

m <- matrix(c(1,1,2,1), 2, 2)
layout(m, widths=c(2,1), heights=c(1,2))
layout.show(2)


## ---- page 073 ----------------------------------------------------------------------------------------------
m <- matrix(c(0:3), 2, 2)
layout(m, c(1,3), c(1,3))
layout.show(3)
m <- matrix(scan(), 5, 5)
0 0 3 3 3 1 1 3 3 3 0 0
3 3 3 0 2 2 0 5 4 2 2 0 5


layout(m)
layout.show(5)



## ---- page 074 ----------------------------------------------------------------------------------------------
par(cex=1.25, mex=1.25) # character (cex) & margin (mex) expansion


## ---- page 075 ----------------------------------------------------------------------------------------------
oldpar <- par(cex=1.25, mex=1.25)


## ---- page 076 ----------------------------------------------------------------------------------------------
attach(elasticband)
oldpar <- par(cex=1.5, mex=1.5)
plot(distance ~ stretch)
par(oldpar) # Restores the earlier settings
detach(elasticband)


oldpar <- par(cex=1.25, mex=1.25)
on.exit(par(oldpar))


## ---- page 077 ----------------------------------------------------------------------------------------------
par(mfrow=c(2,2), pch=16, cex.lab=1, cex=1, mar=c(4,4,4,2))
library(MASS)
attach(Animals)
plot(body, brain)
plot(sqrt(body), sqrt(brain))
plot((body)^.1, (brain)^0.1)
plot(body, brain, log=c("xy"))
detach(Animals)
par(mfrow=c(1,1), pch=1)    # Restore to 1 figure per page


## ---- page 080 ----------------------------------------------------------------------------------------------
library(DAAG)
attach(primates) 	
plot(Bodywt, Brainwt)
text(x=Bodywt, y=Brainwt, labels=1:nrow(primates), adj=0)
# adj=0 implies left adjusted text


## ---- page 082 ----------------------------------------------------------------------------------------------
plot(x=Bodywt, y=Brainwt, pch=16, xlab="Body weight (kg)",
    ylab="Brain weight (g)", xlim=c(0,310), ylim=c(0,1300))
 
# Specify xlim so that there is room for the labels
text(x=Bodywt, y=Brainwt, labels=1:nrow(primates), pos=4)

text(x=Bodywt, y=Brainwt, labels=row.names(primates), pos=2,
    col="red")
detach(primates)


## ---- page 083 ----------------------------------------------------------------------------------------------
par(mar=c(1,1,1,1) + 0.1)
plot(1, 1, xlim=c(1, 7.5), ylim=c(0,5),
    type="n", xlab="", ylab="", axes=F)  	# Do not plot points
points(1:7, rep(4.5, 7), cex=c(1,1:6), col=c(1,1:6), pch=0:6)
text(1:7, rep(3.5, 7), labels=paste(0:6), cex=1:7, col=c(1,1:6))


points(1:7, rep(2,7), pch=(0:6)+7, cex=2) # Plot symbols 7 to 13
text((1:7)+0.25, rep(2,7), paste((0:6)+7), cex=2)
# Label with symbol number
points(1:7, rep(1,7), pch=(0:6)+14, cex=2) # Plot symbols 14 to 20
text((1:7)+0.25, rep(1,7), paste((0:6)+14), cex=2)


## ---- page 085 ----------------------------------------------------------------------------------------------
view.colours <- function(){
 plot(1, 1, xlim=c(0,14), ylim=c(0,3), type="n", axes=F,
      xlab="",ylab="")
 text(1:6, rep(2.5,6), paste(1:6), col=palette()[1:6], cex=2.5)
 text(10, 2.5, "Default palette", adj=0)
 rainchars <- c("R","O","Y","G","B","I","V")
 text(1:7, rep(1.5,7), rainchars, col=rainbow(7), cex=2.5)
 text(10, 1.5, "rainbow(7)", adj=0)
 cmtxt <- substring("cm.colors", 1:9, 1:9)
 text(1:9, rep(0.5,9), cmtxt, col=cm.colors(9), cex=3)
 text(10, 0.5, "cm.colors(9)", adj=0)
}


## ---- page 086 ----------------------------------------------------------------------------------------------
view.colours()


## ---- page 089 ----------------------------------------------------------------------------------------------
library(car)
attach(Florida)
plot(BUSH, BUCHANAN, xlab="Bush", ylab="Buchanan")
County <- rownames(Florida)
identify(BUSH, BUCHANAN, County, col=4)
detach(Florida)


## ---- page 091 ----------------------------------------------------------------------------------------------
attach(Florida) # if not already attached
plot(BUSH, BUCHANAN, xlab="Bush", ylab="Buchanan")
locator()
detach(Florida)


## ---- page 099 ----------------------------------------------------------------------------------------------
x <- seq(0,100,by=.5)
y <- pi*x^2
plot(x, y, xlab="Radius", ylab=expression(Area == pi*r^2), type='l')



## ---- page 102 ----------------------------------------------------------------------------------------------
staz <- function(x){
  m <- mean(x, na.rm=TRUE)
  v <- var(x, na.rm=TRUE)
  val <- (x-m)/sqrt(v)
  return(val)
}
 
a <- 1:10
staz(a)



## ---- page 106 ----------------------------------------------------------------------------------------------
fact.3 <- function(x){  
  if (x <= 1) {
    1 # termination condition
  } else{
    x * fact.3(x - 1) # recursive call
  }
} 
fact.3(5)  


## ---- page 107 ----------------------------------------------------------------------------------------------
fact.4 <- function(x){  
  if (x <= 1) {
    1 # termination condition
  } else{
    x * Recall(x - 1) # recursive call
  }
} 
fact.4(5) 




## ---- page 108 ----------------------------------------------------------------------------------------------
head(iris, 2)
library(doBy)
summaryBy( Sepal.Width + Sepal.Length ~ Species , iris )


## ---- page 109 ----------------------------------------------------------------------------------------------
orderBy(~ Sepal.Width, iris )


## ---- page 110 ----------------------------------------------------------------------------------------------
orderBy(~ Species + Sepal.Width, iris )


## ---- page 111 ----------------------------------------------------------------------------------------------
sampleBy(~  Species, frac =0.1, data = iris )


## ---- page 112 ----------------------------------------------------------------------------------------------
split(iris, iris$Species)


## ---- page 113 ----------------------------------------------------------------------------------------------
lapply(split( iris$Sepal.Length , iris$Species ), mean)


## ---- page 114 ----------------------------------------------------------------------------------------------
subset(iris, Species == "setosa")
subset(iris, Species == "setosa" & Sepal.Length > 5.0)


## ---- page 115 ----------------------------------------------------------------------------------------------
# include Sepal.Length, Species column
subset(iris, select = c(Sepal.Length, Species))
# exclude Sepal.Length, Species column
subset(iris, select = -c(Sepal.Length, Species))


## ---- page 116 ----------------------------------------------------------------------------------------------
x <- data.frame(name = c("a", "b", "c"), math = c(1, 2, 3))
y <- data.frame(name = c("c", "b", "a"), english = c(4, 5, 6))
merge(x, y)
x <- data.frame(name = c("a", "b", "c"), math = c(1, 2, 3))
y <- data.frame(name = c("c", "b", "a"), english = c(4, 5, 6))
merge(x, y, all = TRUE)


## ---- page 117 ----------------------------------------------------------------------------------------------
x <- c(20, 11, 33, 50, 47)
sort(x)
sort(x, decreasing = TRUE)
order(x)
rank(x)


## ---- page 118 ----------------------------------------------------------------------------------------------
iris[order(iris$Sepal.Length), ]
iris[order(iris$Sepal.Length, iris$Petal.Length), ]


## ---- page 119 ----------------------------------------------------------------------------------------------
with( iris, {
  print(mean(Sepal.Length))
  print(mean(Sepal.Width))
})


x <- data.frame(val = c(NA, 1, NA, 2, 3, 4))
x <- within(x, {
  val <- ifelse(is.na(val), median(val, na.rm = TRUE), val)
})
head(x, 4)


## ---- page 120 ----------------------------------------------------------------------------------------------
aggregate(Sepal.Width ~ Species, iris, mean)


## ---- page 121 ----------------------------------------------------------------------------------------------
x <- data.frame(medicine = c("a", "b", "c"), ctl = c(5, 3, 2),
                exp = c(4, 5, 7))
x
stacked_x <- stack(x)
head(stacked_x, 4)


## ---- page 122 ----------------------------------------------------------------------------------------------
library(doBy)
summaryBy(values ~ ind, stacked_x)

unstack(stacked_x, values ~ ind)



## ---- page 123 ----------------------------------------------------------------------------------------------
library(sqldf)
sqldf("select distinct Species from iris")
sqldf("select avg(Sepal_Length) from iris where Species = 'setosa'" )
sqldf("select species, avg(sepal_length) from iris group by species" )


## ---- page 125 ----------------------------------------------------------------------------------------------
head(iris, 2)
apply(iris[ , 1:4], 2, sum)
(x <- list(a=1:5, b=6:9))
lapply(x, mean)



## ---- page 130 ----------------------------------------------------------------------------------------------
set.seed(1)
d <- data.frame(year = rep(2000:2002, each = 3), count = round(runif(9, 0, 20)))
print(d)




## ---- page 131 ----------------------------------------------------------------------------------------------
library(plyr)
ddply(d, "year", function(x) {
  mean.count <- mean(x$count)
  sd.count <- sd(x$count)
  cv <- sd.count/mean.count
  data.frame(cv.count = cv)
})


## ---- page 132 ----------------------------------------------------------------------------------------------
ddply(d, "year", summarise, mean.count = mean(count))


## ---- page 133 ----------------------------------------------------------------------------------------------
ddply(d, "year", transform, total.count = sum(count))


## ---- page 134 ----------------------------------------------------------------------------------------------
ddply(d, "year", mutate, mu = mean(count), sigma = sd(count), 
      cv = sigma/mu)


## ---- page 135 ----------------------------------------------------------------------------------------------
par(mfrow = c(1, 3), mar = c(2, 2, 1, 1), oma = c(3, 3, 0, 0))
d_ply(d, "year", transform, plot(count, main = unique(year), type = "o"))
mtext("count", side = 1, outer = TRUE, line = 1)
mtext("frequency", side = 2, outer = TRUE, line = 1)


## ---- page 136 ----------------------------------------------------------------------------------------------
head(baseball, 2)
baseball.dat <- subset(baseball, year > 2000) # data from the plyr package
x <- ddply(baseball.dat, c("year", "team"), summarize, homeruns = sum(hr))
head(x, 3)


## ---- page 137 ----------------------------------------------------------------------------------------------
f <- function(x) if (x == 1) stop("Error!") else 1
safe.f <- failwith(NA, f, quiet = TRUE)
#llply(1:2, f)
llply(1:2, safe.f)




## ---- page 138 ----------------------------------------------------------------------------------------------
# not available under Windows
library(doMC)
x <- c(1:10)
wait <- function(i) Sys.sleep(0.1)
system.time(llply(x, wait))
system.time(sapply(x, wait))
registerDoMC(2)
system.time(llply(x, wait, .parallel = TRUE))


## ---- page 139 ----------------------------------------------------------------------------------------------
system.time(ddply(baseball, "id", summarize, length(year)))
system.time(tapply(baseball$year, baseball$id, function(x) length(x)))
system.time(ddply(idata.frame(baseball), "id", summarize, length(year)))
library(data.table)
dt <- data.table(baseball, key="id")
system.time(dt[, length(year), by=list(id)])


## ---- page 142 ----------------------------------------------------------------------------------------------
# iris data
str(iris)

# split into training and test datasets
set.seed(1234)
ind <- sample(2, nrow(iris), replace=T, prob=c(0.7, 0.3))
iris.train <- iris[ind==1, ]
iris.test <- iris[ind==2, ]




## ---- page 143 ----------------------------------------------------------------------------------------------
# build a decision tree
library(party)
iris.formula <- Species ~ Sepal.Length + Sepal.Width + 
                          Petal.Length + Petal.Width
iris.ctree <- ctree(iris.formula, data=iris.train)


## ---- page 144 ----------------------------------------------------------------------------------------------
plot(iris.ctree)


## ---- page 145 ----------------------------------------------------------------------------------------------
# predict on test data
pred <- predict(iris.ctree, newdata = iris.test)

# check prediction result
table(pred, iris.test$Species)


## ---- page 147 ----------------------------------------------------------------------------------------------
library(animation)
par(mar=c(3, 3, 1,3), mgp=c(1.5, 0.5, 0), bg="white")
layout(matrix(1:8, ncol=4))
cent <- 1.5 * c(1, 1, -1, -1, 1, -1, 1, -1)
x <- NULL
set.seed(131)
for (i in 1:8) {
  x <- c(x, rnorm(25, mean=cent[i]))
}
x <- matrix(x, ncol=2)
colnames(x) <- c("X1", "X2")
kmeans.ani(x, centers=4, pch=1:4, col=1:4)


## ---- page 149 ----------------------------------------------------------------------------------------------
set.seed(8953)
iris2 <- iris
# remove class IDs
iris2$Species <- NULL
# k-means clustering
iris.kmeans <- kmeans(iris2, 3)
# check result
table(iris$Species, iris.kmeans$cluster)




## ---- page 150 ----------------------------------------------------------------------------------------------
# plot clusters and their centers
plot(iris2[c("Sepal.Length", "Sepal.Width")], col=iris.kmeans$cluster)
points(iris.kmeans$centers[, c("Sepal.Length", "Sepal.Width")],
       col=1:3, pch="*", cex=5)




## ---- page 151 ----------------------------------------------------------------------------------------------
library(fpc)
iris2 <- iris[-5] # remove class IDs
# DBSCAN clustering
ds <- dbscan(iris2, eps = 0.42, MinPts = 5)
# compare clusters with original class IDs
table(ds$cluster, iris$Species)


## ---- page 152 ----------------------------------------------------------------------------------------------
# 1-3: clusters; 0: outliers or noise
plotcluster(iris2, ds$cluster)


## ---- page 153 ----------------------------------------------------------------------------------------------
library(amap)
model <- hclusterpar(iris[,-5], method="euclidean", link="ward", nbproc=1)
plot(model, main="Cluster Dendrogram", xlab="",hang=-1, labels=iris$Species)
#Add in rectangles to show the clusters.
rect.hclust(model, k=3) #cutree(model, k=3)


## ---- page 155 ----------------------------------------------------------------------------------------------
#ojtpath <- "~/Dropbox/OJT/"
#if(.Platform$OS.type=="windows")  ojtpath<-"C:/Users/jang/Documents/OJT/"
load(paste(ojtpath, "R/Data/titanic.raw.rdata", sep=""))
dim(titanic.raw)
idx <- sample(1:nrow(titanic.raw), 8)
titanic.raw[idx, ]




## ---- page 156 ----------------------------------------------------------------------------------------------
# find association rules with the APRIORI algorithm
library(arules)
rules <- apriori(titanic.raw, control=list(verbose=F), 
                parameter=list(minlen=2, supp=0.005, conf=0.8),
appearance=list(rhs=c("Survived=No", "Survived=Yes"), default="lhs"))
# sort rules
quality(rules) <- round(quality(rules), digits=3)
rules.sorted <- sort(rules, by="lift")
# have a look at rules inspect(rules.sorted)


## ---- page 157 ----------------------------------------------------------------------------------------------
# have a look at rules 
inspect(rules.sorted)




## ---- page 158 ----------------------------------------------------------------------------------------------
library(arulesViz)
plot(rules, method = "graph")




## ---- page 163 ----------------------------------------------------------------------------------------------
library(twitteR)
#consumerKey <- "YOUR API key"
#consumerSecret <- "API secret"

credential <- OAuthFactory$new(
  consumerKey=consumerKey,
  consumerSecret=consumerSecret,
  requestURL="https://api.twitter.com/oauth/request_token",
  accessURL="https://api.twitter.com/oauth/access_token",
   authURL="https://api.twitter.com/oauth/authorize")

download.file(url="http://curl.haxx.se/ca/cacert.pem",
              destfile="cacert.pem")


## ---- page 164 ----------------------------------------------------------------------------------------------
credential$handshake(cainfo="cacert.pem")



## ---- page 165 ----------------------------------------------------------------------------------------------
registerTwitterOAuth(credential)



## ---- page 166 ----------------------------------------------------------------------------------------------
val <- iconv("교통연구원", "CP949", "UTF-8")
result <- searchTwitter(val, since='2014-01-01', until='2014-05-20', n=1000, cainfo="cacert.pem")
result[1:3]



## ---- page 167 ----------------------------------------------------------------------------------------------
## retrieve tweets from the user timeline of @rdatammining
library(twitteR)

tweets <- userTimeline('rdatamining', n=400, cainfo="cacert.pem")
# save(tweets, file=paste(ojtpath, "R/Data/rdmTweets-201405.RData", sep="" ))
# load(paste(ojtpath, "R/Data/rdmTweets-201405.RData", sep="" ))

(nDocs <- length(tweets))
strwrap(tweets[[320]]$text, width = 50)
# convert tweets to a data frame
df <- do.call("rbind", lapply(tweets, as.data.frame))




## ---- page 168 ----------------------------------------------------------------------------------------------
library(tm)
# build a corpus
myCorpus <- Corpus(VectorSource(df$text))
# convert to lower case
myCorpus <- tm_map(myCorpus, tolower)
# remove punctuation & numbers
myCorpus <- tm_map(myCorpus, removePunctuation)
myCorpus <- tm_map(myCorpus, removeNumbers)
# remove URLs
removeURL <- function(x) gsub("http[[:alnum:]]*", "", x)
myCorpus <- tm_map(myCorpus, removeURL)
# remove 'r' and 'big' from stopwords
myStopwords <- setdiff(stopwords("english"), c("r", "big"))
# remove stopwords
myCorpus <- tm_map(myCorpus, removeWords, myStopwords)


## ---- page 169 ----------------------------------------------------------------------------------------------
# keep a copy of corpus
library(SnowballC)
myCorpusCopy <- myCorpus
# stem words
myCorpus <- tm_map(myCorpus, stemDocument)
# stem completion
myCorpus <- tm_map(myCorpus, stemCompletion,
dictionary = myCorpusCopy)
# replace "miners" with "mining", because "mining" was
# first stemmed to "mine" and then completed to "miners"
myCorpus <- tm_map(myCorpus, gsub, pattern="miners",
replacement="mining")
strwrap(myCorpus[320], width=50)


## ---- page 170 ----------------------------------------------------------------------------------------------
myTdm <- TermDocumentMatrix(myCorpus, control=list(wordLengths=c(1,Inf)))

# inspect frequent words
(freq.terms <- findFreqTerms(myTdm, lowfreq=20))


## ---- page 171 ----------------------------------------------------------------------------------------------
# which words are associated with 'r'?
findAssocs(myTdm, "r", 0.2)
# which words are associated with 'mining'?
findAssocs(myTdm, "mining", 0.25)




## ---- page 172 ----------------------------------------------------------------------------------------------
library(wordcloud); library(RColorBrewer); par(mar=c(0,0,0,0))
m <- as.matrix(myTdm)
freq <- sort(rowSums(m), decreasing=TRUE)
wordcloud(words=names(freq), freq=freq, min.freq=4, random.order=FALSE, 
          scale=c(3, .99), colors=c("grey", rev(brewer.pal(8, "Set1"))))


## ---- page 173 ----------------------------------------------------------------------------------------------
library(topicmodels)
set.seed(123)
myLda <- LDA(as.DocumentTermMatrix(myTdm), k=8)
terms(myLda, 3)


## ---- page 174 ----------------------------------------------------------------------------------------------
library(rattle)
rattle()




## ---- page 187 ----------------------------------------------------------------------------------------------
library(XML)
# OAuthid <- "Your Oauthid"
doc3 <- xmlTreeParse(paste("http://openapi.seoul.go.kr:8088/", OAuthid, 
            "/xml/CardBusStatisticsService/1/200/201402/5516/", sep=""))
docroot <- xmlRoot(doc3)
doclist <- xmlToList(docroot)
doclist <- doclist[3:length(doclist)]
doclist <- lapply(doclist, unlist)
doclist <- Reduce(rbind, doclist)
doclist <- as.data.frame(doclist, row.names = NA, stringsAsFactors = FALSE)
doclist[, c(1:3, 7:9)] <- apply(doclist[, c(1:3, 7:9)], 2, as.numeric)
Encoding(doclist[, "BUS_ROUTE_NM"]) <- "UTF-8"
Encoding(doclist[, "BUS_STA_NM"]) <- "UTF-8"

head(doclist, 2)



## ---- page 188 ----------------------------------------------------------------------------------------------
library(RColorBrewer)
cols <- brewer.pal(8, "Set1")[1:2]
tmp <- doclist[, c("RIDE_PASGR_NUM",  "ALIGHT_PASGR_NUM")]
rownames(tmp) <- paste(1:nrow(doclist), doclist[,"BUS_STA_NM"], sep="_")
par(cex.axis=.5, mar=c(9,4,4,2))
barplot(t(as.matrix(tmp)), col=cols, las=3, axes=FALSE)
par(cex.axis=2)
axis(2)
legend("topleft", legend=c("RIDE", "ALIGHT"), col=cols, pch=15, cex=2)




