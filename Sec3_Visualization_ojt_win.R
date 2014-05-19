###########################################################################
## OJT 2014 in KOTI
## R script for Sec 3
###########################################################################



loadedOnly <- scan(what="list") 
"bitops"     "boot"       "caTools"    "coda"       "colorspace" "deldir"     "digest"     "gtable"     "httpuv"     "labeling"  
"LearnBayes" "MASS"       "munsell"    "nlme"       "plyr"       "png"        "proto"      "Rcpp"       "splines"    "stringr"   
"tools"      "whisker"    "xtable"     "yaml"    



otherPkgs <- scan(what="list") 
"leafletR"       "shapefiles"     "foreign"        "plotGoogleMaps" "lattice"        "scales"         "shiny"         
"WDI"            "RJSONIO"        "rjson"          "googleVis"      "reshape2"       "rCharts"        "RColorBrewer"  
"spdep"          "Matrix"         "rgeos"          "raster"         "mapproj"        "maps"           "ggmap"         
"ggplot2"        "PBSmapping"     "rgdal"          "maptools"       "sp"             "RgoogleMaps"  


destdir <- "C:/Users/jang/Documents/OJT/R/packages"
download.packages(otherPkgs, destdir)
download.packages(loadedOnly, destdir)

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
pkglist <- dir(destdir, pattern="*.zip")
oldwd <- getwd()
setwd("C:/Users/jang/Documents/OJT/R/packages")
for(i in 1:length(pkglist)){
  install.packages(pkglist[i])
}
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################


## Install required packages
reqpkgs <- c("RgoogleMaps", "maptools", "rgdal", "ggplot2", "ggmap", "grid",
             "sp", "raster", "rgeos", "spdep","googleVis", "rjson", "shiny", 
             "leafletR","RColorBrewer", "scales", "plotGoogleMaps", "sp",
             "shapefiles", "foreign", "reshape2", "WDI", "datasets")
insedpkgs <- installed.packages()[,1]
inspkg <- reqpkgs[!reqpkgs %in% insedpkgs]
if(length(inspkg) > 0){
  install.packages(inspkg)
}

## Page 003 ---------------------------------------------------------------
library(RgoogleMaps)
ojtpath <- "~/Dropbox/OJT/"
if(.Platform$OS.type=="windows"){
  user <- shell("echo %username%", intern=TRUE)
  ojtpath <- paste("C:/Users/", user, "/Documents/OJT/", sep="")
}
figpath <- paste(ojtpath, "R/figure/", sep="")
MyMap <- GetMap(center=c(0,0), zoom =1,  size = c(640, 410), 
                destfile = paste(figpath, "World1.png", sep=""))
PlotOnStaticMap(MyMap,  size = c(640, 410))


## Page 004 ---------------------------------------------------------------
# cChange the center
MyMap <- GetMap(center=c(37.56654, 126.978), zoom =1,  size = c(640, 390), 
                destfile = paste(figpath, "World2.png", sep=""))
PlotOnStaticMap(MyMap,  size = c(640, 390))


## Page 005 ---------------------------------------------------------------
# Change the zoom variable
par(mfrow=c(1,2))
MyMap <- GetMap(center = c(37.56654, 126.978), zoom =10, size = c(640, 510), 
                destfile = paste(figpath, "World3.png", sep=""))
PlotOnStaticMap(MyMap,  size = c(640, 510))

MyMap <- GetMap(center = c(37.56654, 126.978), zoom =13, size = c(640, 510), 
                destfile = paste(figpath, "World4.png", sep=""))
PlotOnStaticMap(MyMap,  size = c(640, 510))


## Page 006 ---------------------------------------------------------------
# Change the size
par(mfrow=c(1,2))
MyMap <- GetMap(center = c(37.56654, 126.978), zoom =10, size = c(500, 500), 
                destfile = paste(figpath, "World500.png", sep=""))
PlotOnStaticMap(MyMap,  size = c(640, 500))

MyMap <- GetMap(center = c(37.56654, 126.978), zoom =10, size = c(200, 200), 
                destfile = paste(figpath, "World200.png", sep=""))
PlotOnStaticMap(MyMap,  size = c(640, 200))


## Page 007 ---------------------------------------------------------------
# Choose the right maptype
par(mfrow=c(1,3))
MyMap <- GetMap(center = c(37.56654, 126.978), zoom =13, size = c(640, 510), 
                destfile = paste(figpath, "Worlds.png", sep=""), maptype = "satellite")
PlotOnStaticMap(MyMap,  size = c(640, 510))

MyMap <- GetMap(center = c(37.56654, 126.978), zoom =13, size = c(640, 510), 
                destfile = paste(figpath, "Worldr.png", sep=""), maptype = "roadmap")
PlotOnStaticMap(MyMap,  size = c(640, 510))

MyMap <- GetMap(center = c(37.56654, 126.978), zoom =13, size = c(640, 510), 
                destfile = paste(figpath, "Worldt.png", sep=""), maptype = "terrain")
PlotOnStaticMap(MyMap,  size = c(640, 510))


## Page 008 ---------------------------------------------------------------
# Make a plot in gray scale
par(mfrow=c(1,1))
MyMap <- GetMap(center = c(37.56654, 126.978), zoom =13, size = c(640, 510),
                destfile = paste(figpath, "Worldgrey.png", sep=""),
                maptype = "terrain", GRAYSCALE=TRUE)
PlotOnStaticMap(MyMap,  size = c(640, 510))


## Page 009 ---------------------------------------------------------------
# Make a map including several points
# qbbox
bb <- qbbox(c(37.67838, 37.3827, 37.7381),c(126.7712, 127.1189, 127.0337))
print(bb)
MyMap <- GetMap.bbox(bb$lonR, bb$latR, size = c(640, 450), maptype = "roadmap",
                     destfile = paste(figpath, "seoul3.png", sep=""))
PlotOnStaticMap(MyMap,  size = c(640, 450))


## Page 011 ---------------------------------------------------------------
# Make plot on map:
bb <- qbbox(c(37.67838, 37.3827, 37.7381),c(126.7712, 127.1189, 127.0337))
MyMap <- GetMap.bbox(bb$lonR, bb$latR, destfile = paste(figpath, "seoul4.png", sep=""),
                     maptype = "roadmap", size = c(640, 450), hl="ko")
 
#Define the markers:
mymarkers <- cbind.data.frame(lat = c(37.51551, 37.55453), lon = c(126.9078, 126.9707))
 
# adding points
tmp <- PlotOnStaticMap(MyMap, lat = mymarkers[,"lat"], lon = mymarkers[,"lon"],
                       cex=2.5,pch=20,col=c("cyan", "brown1"), add=F)
# adding line
tmp <- PlotOnStaticMap(MyMap, lat = mymarkers[,"lat"], lon = mymarkers[,"lon"],
                       col=c("blueviolet"), add=TRUE, FUN = lines, lwd = 4)


## Page 013 ---------------------------------------------------------------
# Add border(or polygon) of counties
library(maptools)
library(rgdal)
nc <- readShapePoly(paste(ojtpath, "R/Data/시군구_utf8.shp", sep=""),
                    proj4string=CRS("+proj=tmerc +lat_0=38 +lon_0=127 +k=1 +x_0=200000
                                    +y_0=500000 +ellps=bessel +units=m +no_defs"))
nc2 <- spTransform(nc, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
MyMap <- GetMap.bbox(c(124.5782, 130.9558), c(33.16315, 38.64771),
                      destfile = paste(figpath, "skorea.png", sep=""),
                      maptype = "roadmap", size = c(640, 640), hl="ko")
shp <- nc2@polygons[[1]]@Polygons[[1]]@coords[, 2:1]
colnames(shp) <- c("Y", "X")
PlotPolysOnStaticMap(MyMap, shp, lwd=.5, add = F, verbose=.5, col=NA)
for(i in 2:length(nc2@polygons)){
 	for(j in 1:length(nc2@polygons[[i]]@Polygons)){
    	shp <- nc2@polygons[[i]]@Polygons[[j]]@coords[, 2:1]
    	colnames(shp) <- c("Y", "X")
    	PlotPolysOnStaticMap(MyMap, shp, lwd=.5, add = T, verbose=.5, col=NA)
  	}
}

###########################################################################
## ggplot2
###########################################################################
## Page 023 ---------------------------------------------------------------
library(ggplot2)


## Page 024 ---------------------------------------------------------------
head(iris, n = 10) # we can also explicitly set the number of rows to display

table(iris$Species)


## Page 025 ---------------------------------------------------------------
qplot(Sepal.Length, Petal.Length, data = iris)
# Plot Sepal.Length vs. Petal.Length, using data from the `iris` data frame.
# * First argument `Sepal.Length` goes on the x-axis.
# * Second argument `Petal.Length` goes on the y-axis.
# * `data = iris` means to look for this data in the `iris` data frame.    


## Page 026 ---------------------------------------------------------------
qplot(Sepal.Length, Petal.Length, data = iris, color = Species) 


## Page 027 ---------------------------------------------------------------
qplot(Sepal.Length, Petal.Length, data = iris, color = Species, size = Petal.Width)


## Page 028 ---------------------------------------------------------------
qplot(Sepal.Length, Petal.Length, data = iris, color = Species, size = Petal.Width, 
      alpha = I(0.7))


## Page 029 ---------------------------------------------------------------
qplot(Sepal.Length, Petal.Length, data = iris, color = Species,
      xlab = "Sepal Length", ylab = "Petal Length",
      main = "Sepal vs. Petal Length in Fisher's Iris data")


## Page 030 ---------------------------------------------------------------
# # These two invocations are equivalent.
qplot(Sepal.Length, Petal.Length, data = iris, geom = "point")
qplot(Sepal.Length, Petal.Length, data = iris)


## Page 031 ---------------------------------------------------------------
movies = data.frame(
  director = c("spielberg", "spielberg", "spielberg", "jackson", "jackson"),
  movie = c("jaws", "avatar", "schindler's list", "lotr", "king kong"),
  minutes = c(124, 163, 195, 600, 187)
)
# Plot the number of movies each director has.
qplot(director, data = movies, geom = "bar", ylab = "# movies", fill=factor(movie) )
# By default, the height of each bar is simply a count.


## Page 032 ---------------------------------------------------------------
# Here the height of each bar is the total running time of the director's movies.
qplot(director, weight = minutes, data = movies, geom = "bar", ylab = "total length (min.)",
      fill=factor(movie))


## Page 033 ---------------------------------------------------------------
qplot(Sepal.Length, Petal.Length, data = iris, geom = "line", color = Species)
# Using a line geom doesn't really make sense here, but hey.


## Page 034 ---------------------------------------------------------------
# `Orange` is another built-in data frame that describes the growth of orange trees.
qplot(age, circumference, data = Orange, geom = "line",
      colour = Tree, 
      main = "How does orange tree circumference vary with age?")


## Page 035 ---------------------------------------------------------------
# We can also plot both points and lines.
qplot(age, circumference, data = Orange, geom = c("point", "line"), colour = Tree)


## Page 036 ---------------------------------------------------------------
#histogram
qplot(carat, data=diamonds, geom="histogram")
qplot(carat, data=diamonds, geom="histogram", binwidth=0.1)
qplot(carat, data=diamonds, geom="histogram", binwidth=0.01)


## Page 039 ---------------------------------------------------------------
head(mpg)


## Page 040 ---------------------------------------------------------------
p <- ggplot(data=mpg, mapping=aes(x=cty, y=hwy))
p + geom_point()


## Page 041 ---------------------------------------------------------------
summary(p)
summary(p+geom_point())


## Page 042 ---------------------------------------------------------------
p + geom_point(color="red4", size=3)


## Page 043 ---------------------------------------------------------------
p <- ggplot(data=mpg, mapping=aes(x=cty, y=hwy, colour=factor(year)))
p + geom_point()


## Page 044 ---------------------------------------------------------------
p + geom_point() + stat_smooth()


## Page 045 ---------------------------------------------------------------
p <- ggplot(data=mpg, mapping=aes(x=cty, y=hwy))
p + geom_point(aes(colour=factor(year))) + stat_smooth()


## Page 046 ---------------------------------------------------------------
p <- ggplot(data=mpg, mapping=aes(x=cty, y=hwy))
p + geom_point(aes(colour=factor(year))) +
    stat_smooth() +
    scale_color_manual(values= c("steelblue", "red4"))


## Page 047 ---------------------------------------------------------------
p <- ggplot(data=mpg, mapping=aes(x=cty, y=hwy))
p + geom_point(aes(colour=factor(year), size=displ)) +
    stat_smooth() +
    scale_color_manual(values= c("steelblue", "red4"))


## Page 048 ---------------------------------------------------------------
p + geom_point(aes(colour=factor(year),size=displ), alpha=0.5, position = "jitter")+
    stat_smooth()+
    scale_color_manual(values =c("steelblue","red4"))+
    scale_size_continuous(range = c(4, 10))


## Page 049 ---------------------------------------------------------------
p + geom_point(aes(colour=factor(year),size=displ), alpha=0.5, position = "jitter")+
    stat_smooth()+
    scale_color_manual(values =c("steelblue","red4"))+
    scale_size_continuous(range = c(4, 10)) + coord_flip()


## Page 050 ---------------------------------------------------------------
p + geom_point(aes(colour=factor(year),size=displ), alpha=0.5, position = "jitter")+
    stat_smooth()+
    scale_color_manual(values =c("steelblue","red4")) +
    scale_size_continuous(range = c(4, 10))  +
    coord_polar()


## Page 051 ---------------------------------------------------------------
p + geom_point(aes(colour=factor(year),size=displ), alpha=0.5, position = "jitter")+
    stat_smooth()+
    scale_color_manual(values =c("steelblue","red4")) +
    scale_size_continuous(range = c(4, 10))  +
    coord_cartesian(xlim = c(15, 25), ylim=c(15,40))


## Page 052 ---------------------------------------------------------------
p + geom_point(aes(colour=class, size=displ), alpha=0.5, position = "jitter")+
    stat_smooth() +
    scale_size_continuous(range = c(4, 10))  +
    facet_wrap( ~ year, ncol=1)


## Page 053 ---------------------------------------------------------------
p.tmp <- ggplot(mtcars, aes(x=factor(1), fill=factor(cyl))) + geom_bar(width=1)
p.tmp + coord_polar(theta="y")
p.tmp + coord_polar()
ggplot(mtcars, aes(factor(cyl), fill=factor(cyl))) + geom_bar(width=1) + coord_polar()


## Page 054 ---------------------------------------------------------------
ggplot(mpg, aes(class, hwy)) + geom_boxplot()
ggplot(mpg, aes(class, hwy, fill = factor(year))) +
  geom_boxplot()


## Page 055 ---------------------------------------------------------------
#reorder class according to median(hwy)
ggplot(mpg, aes(reorder(class, hwy, median), hwy, fill = factor(year))) +
  geom_boxplot()


## Page 056 ---------------------------------------------------------------
ggplot(mpg, aes(displ, hwy, color = factor(cyl))) +
    geom_point() +
    stat_smooth(method = "lm")

###########################################################################
## ggmap
###########################################################################
## Page 057 ---------------------------------------------------------------
library(ggmap)
koti <- "315 Goyangdaero"
qmap(koti, zoom = 14, source = "osm", scale = 20000)


## Page 061 ---------------------------------------------------------------
#lonlat <-  geocode(c("대화역", "주엽역", "정발산역", "마두역", "백석역"))
locval <- c("대화역", "주엽역", "정발산역", "마두역", "백석역")
if(localeToCharset() == "CP949") locval <- iconv(locval, "CP949", "UTF-8")
lonlat <-  geocode(locval)
qmplot(lon, lat, data =lonlat, colour = I("red"), size = I(5), darken = 0,
       verbose = TRUE, source = "osm", scale = 40000)


## Page 062 ---------------------------------------------------------------
locval <- "대화역"
if(localeToCharset() == "CP949") locval <- iconv(locval, "CP949", "UTF-8")
( gc <- geocode(locval, output = "more"))
revgeocode(as.numeric(gc[1:2]))


## Page 063 ---------------------------------------------------------------
from <- c("houston", "houston", "dallas")
to <- c("waco, texas", "san antonio", "houston")
mapdist(from, to)


## Page 065 ---------------------------------------------------------------
legs_df <- route(
    "marrs mclean science, baylor university",
    "220 south 3rd street, waco, tx 76701a", # National Cancer Center
    alternatives = TRUE)

qmap("424 clay avenue, waco, tx", zoom = 15, maptype = "hybrid",
      base_layer = ggplot(aes(x = startLon, y = startLat), data = legs_df)) +
   geom_leg(
     aes(x = startLon, y = startLat, xend = endLon, yend = endLat, colour = route),
     alpha = 3/4, size = 2, data = legs_df) +
   labs(x = "Longitude", y = "Latitude", colour = "Route") +
   facet_wrap(~ route, ncol = 3) + theme(legend.position = "top")



## Page 066 ---------------------------------------------------------------
data(crime)
head(crime, n=3L)


## Page 067 ---------------------------------------------------------------
qmap("houston", zoom = 13)
gglocator(2)


violent_crimes <- subset(crime, offense != "auto theft" & offense != "theft" &
                                offense != "burglary")

violent_crimes$offense <- factor(violent_crimes$offense,
                                 levels = c('robbery', 'aggravated assault', 'rape', 'murder'))

violent_crimes <- subset(violent_crimes, -95.39681 <= lon & lon <= -95.34188 & 
                                         29.73631 <= lat & lat <= 29.78400)


## Page 069 ---------------------------------------------------------------
library(grid)
theme_set(theme_bw(16))
HoustonMap <- qmap("houston", zoom = 14, color = "bw", legend = "topleft")
crimlabel <- c("Robery","Aggravated Assault","Rape","Murder")
HoustonMap +
   geom_point(aes(x = lon, y = lat, colour = offense, size = offense),
              data = violent_crimes) +
   scale_colour_discrete('Offense', labels = crimlabel) +
   scale_size_discrete('Offense', labels = crimlabel, range = c(1.75,6)) +
   guides(size = guide_legend(override.aes = list(size = 6))) +
   theme(legend.key.size = unit(1.8,'lines'),
         legend.title = element_text(size = 16, face = 'bold'),
         legend.text = element_text(size = 14)) +
   labs(colour = 'Offense', size = 'Offense')


## Page 070 ---------------------------------------------------------------
crimlabel <- c("Robery","Aggravated Assault","Rape","Murder")
HoustonMap +
   stat_bin2d(aes(x = lon, y = lat, colour = offense, fill = offense),
              size = .5, bins = 30, alpha = 1/2, data = violent_crimes) +
   scale_colour_discrete('Offense', labels = crimlabel, guide = FALSE) +
   scale_fill_discrete('Offense', labels = crimlabel) +
   theme(
     legend.text = element_text(size = 15, vjust = .5),
     legend.title = element_text(size = 15, face='bold'),
     legend.key.size = unit(1.8, 'lines')
   )


## Page 072 ---------------------------------------------------------------
houston <- get_map("houston", zoom = 14)
HoustonMap <- ggmap(houston, extent = "device", legend = "topleft")
 
HoustonMap +
   stat_density2d(aes(x = lon, y = lat, fill = ..level.., alpha = ..level..),
                  size = 2, bins = 4, data = violent_crimes, geom = 'polygon') +
   scale_fill_gradient('Violent\nCrime\nDensity') +
   scale_alpha(range = c(.4, .75), guide = FALSE) +
   guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))
 
 overlay <- stat_density2d(
   aes(x = lon, y = lat, fill = ..level..),
   bins = 4, geom = "polygon", alpha = I(.4),
   data = violent_crimes
 )
HoustonMap + 
   overlay +
   scale_fill_gradient('Violent\nCrime\nDensity') + 
   inset(grob = ggplotGrob(ggplot() + overlay + theme_inset()),
         xmin = -95.35836, xmax = Inf, ymin = -Inf, ymax = 29.75062)
 


## Page 074 ---------------------------------------------------------------
houston <- get_map(location = "houston", zoom = 14, color = "bw",
                    source = "osm")
HoustonMap <- ggmap(houston, base_layer = ggplot(aes(x = lon, y = lat),
                                                  data = violent_crimes))
HoustonMap +
   stat_density2d(aes(x = lon, y = lat, fill = ..level..),
                  bins = 5, geom = "polygon", alpha = I(.3),
                  data = violent_crimes) +
   scale_fill_gradient('Violent\nCrime\nDensity',low = "black", high = "red") +
   facet_wrap(~ day, ncol =4)


###########################################################################
## Geospatial Packages
###########################################################################
## Page 076 ---------------------------------------------------------------
library(sp)  # vector data
library(raster)  # raster data
library(rgdal)  # input/output, projections
library(rgeos)  # geometry ops
library(spdep)  # spatial dependenc


## Page 077 ---------------------------------------------------------------
x <- c(11.515, 7.056, 12.945, 12.793, 12.888)
y <- c(24.52, 27.11, 30.09, 24.7, 28.24)
data <- data.frame(Z = c("d", "a", "c", "e", "b"))

# points from scratch
coords = cbind(x, y)

## create SpatialPoints
sp = SpatialPoints(coords)

# make spatial data frame
spdf = SpatialPointsDataFrame(coords, data) 
# spdf = SpatialPointsDataFrame(sp, data)

# promote data frame to spatial
coordinates(data) = cbind(x, y) 
##coordinates(data) = ~ x + y


## Page 079 ---------------------------------------------------------------
# back to data
as.data.frame(data)
data@data
bbox(spdf)


## Page 080 ---------------------------------------------------------------
c1 <- cbind(c(1, 2.2, 3, 4, 6), c(1.5, 1.7, 1.3, 1.4, 1.7))
c2 <- cbind(c(3.5, 5, 6.2), c(.2, 1, 1.1))
c3 <- cbind(c(.8, 2.1, 3.5, 5, 6.2), c(0, 0.3, .1, .12, .15))
# simple line strings
L1 <- Line(c1)
L2 <- Line(c2)
L3 <- Line(c3)


## Page 081 ---------------------------------------------------------------
Ls1 <- Lines(list(L1), ID="a")
Ls2 <- Lines(list(L2,L3), ID="b")


## Page 082 ---------------------------------------------------------------
SL12 <-  SpatialLines(list(Ls1, Ls2))


## Page 083 ---------------------------------------------------------------
SLDF = SpatialLinesDataFrame(SL12, data.frame(Z = c("Road", "River"), row.names = c("a", "b")))
as.data.frame(SLDF)
SpatialLinesLengths(SLDF)




## Page 085 ---------------------------------------------------------------
# single ring feature
x1 <- c(53, 6, 49, 118, 129, 101, 86, 51)
y1 <- c(117, 158, 200, 190, 126, 116, 161, 150)
c1 <- cbind(x1, y1)
r1 <- rbind(c1, c1[1, ])
# join
P1 = Polygon(r1)
Ps1 = Polygons(list(P1), ID = "a")
library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
plot(SpatialPolygons(list(Ps1)), col=cols[1])


## Page 086 ---------------------------------------------------------------
# double ring feature
x2a <- c(155, 173, 267, 274, 263, 230, 220)
y2a <- c(172, 212, 203, 182, 164, 164, 194)
x2b <- c(166, 208, 220, 204, 198, 186)
y2b <- c(164, 178, 156, 154, 163, 158)
# double ring feature
c2a <- cbind(x2a, y2a)
r2a <- rbind(c2a, c2a[1, ])
c2b <- cbind(x2b, y2b)
r2b <- rbind(c2b, c2b[1, ])
P2a <- Polygon(r2a)
P2b <- Polygon(r2b)
Ps2 <- Polygons(list(P2a, P2b), ID = "b")
plot(SpatialPolygons(list(Ps2)), col=cols[2])


## Page 087 ---------------------------------------------------------------
# single ring with hole
x1 <- c(127, 202, 239, 246, 234, 202, 174, 151)
y1 <- c(95, 136, 126, 105,  75,  87,  84,  73)
xh1 <- c(192, 190, 196, 214, 215, 204)
yh1 <- c(96, 107, 116, 113,  93, 100)
c1 <- cbind(x1, y1)
r1 <- rbind(c1, c1[1, ])
P1 <- Polygon(r1)
# single ring with hole
hc1 <- cbind(xh1, yh1)
hr1 <- rbind(hc1, hc1[1, ])
H1 <- Polygon(hr1, hole = TRUE)
P1h <- Polygons(list(P1, H1), ID = "c")
SP1h <- SpatialPolygons(list(P1h))
plot(SP1h, col=cols[3], usePolypath=TRUE)


## Page 088 ---------------------------------------------------------------
# Spatial Polygons Data Frame
SPs = SpatialPolygons(list(Ps1, Ps2, P1h))

SPDF = SpatialPolygonsDataFrame(SPs, 
        data.frame(N = c("one", "two", "three"), 
                   row.names = c("a", "b", "c")))
SPDF@data
plot(SPDF, col=cols[1:3])


## Page 094 ---------------------------------------------------------------
lanepath <- read.csv(paste(ojtpath, "R/Data/bus88.csv", sep=""))
head(lanepath, 4)
L88 <- Line(lanepath)
Ls88 <- Lines(list(L88), ID="bus88")
SL88 <- SpatialLines(list(Ls88))
SLDF88 <-  SpatialLinesDataFrame(SL88, data.frame(Z = c("Busroute88"), row.names = c("bus88")))
proj4string(SLDF88)

proj4string(SLDF88) <- CRS("+init=epsg:5179")
# CRS("+proj=tmerc +lat_0=38 +lon_0=127.5 +k=0.9996 +k=1 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +units=m +no_defs")
# Transform Spatial*
SLDF88wgs84 <- spTransform(SLDF88, CRS("+init=epsg:4326"))
# CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
plot(SLDF88wgs84)


## Page 095 ---------------------------------------------------------------
library(maptools)
spath <- "/Users/jang/Dropbox/OJT/R/Data"

if(.Platform$OS.type=="windows")  spath <-paste(ojtpath, "/R/Data", sep="")

shapes <-  readShapeSpatial(paste(spath, "/bus88.shp", sep=""))
# read/write shapefiles (and others)
#- list formats
head(ogrDrivers(), 4)
shapes <- readOGR(spath, "bus88")
writeOGR(shapes, spath, "busroute88", "ESRI Shapefile")
writeOGR(shapes, paste(spath, "/busroute88.kml", sep=""), "bus88", "KML")


## Page 099 ---------------------------------------------------------------
library("rCharts")


## Page 100 ---------------------------------------------------------------
## Example 1 Facetted Scatterplot
names(iris) = gsub("\\.", "", names(iris))
rPlot(SepalLength ~ SepalWidth | Species, data = iris, 
      color = 'Species', type = 'point')

## Example 2 Facetted Barplot
hair_eye = as.data.frame(HairEyeColor)
rPlot(Freq ~ Hair | Eye, color = 'Eye', data = hair_eye, 
      type = 'bar')


## Page 103 ---------------------------------------------------------------
r1 <- rPlot(mpg ~ wt | am + vs, data = mtcars, type = "point", color = "gear")
r1$print("chart1")
r1


## Page 104 ---------------------------------------------------------------
data(economics, package = "ggplot2")
econ <- transform(economics, date = as.character(date))
m1 <- mPlot(x = "date", y = c("psavert", "uempmed"), type = "Line", data = econ)
m1$set(pointSize = 0, lineWidth = 1)
m1


## Page 105 ---------------------------------------------------------------
hair_eye_male <- subset(as.data.frame(HairEyeColor), Sex == "Male")
n1 <- nPlot(Freq ~ Hair, group = "Eye", data = hair_eye_male, type = "multiBarChart")
n1$print("chart3")
n1


## Page 106 ---------------------------------------------------------------
uspexp <- melt(USPersonalExpenditure)
names(uspexp)[1:2] = c("category", "year")
x1 <- xPlot(value ~ year, group = "category", data = uspexp, type = "line-dotted")
x1$print("chart4")
x1


## Page 107 ---------------------------------------------------------------
h1 <- hPlot(x = "Wr.Hnd", y = "NW.Hnd", data = MASS::survey, 
            type = c("line", "bubble", "scatter"), group = "Clap", size = "Age")
# Warning message:
# In hPlot(x = "Wr.Hnd", y = "NW.Hnd", data = MASS::survey, type = c("line",  :
#   Observations with NA has been removed
h1$print("chart5")
h1


## Page 108 ---------------------------------------------------------------
map3 <- Leaflet$new()
map3$setView(c(51.505, -0.09), zoom = 13)
map3$marker(c(51.5, -0.09), bindPopup = "<p> Hi. I am a popup </p>")
map3$marker(c(51.495, -0.083), bindPopup = "<p> Hi. I am another popup </p>")
map3$print("chart7")
map3


## Page 109 ---------------------------------------------------------------
usp = reshape2::melt(USPersonalExpenditure)
# get the decades into a date Rickshaw likes
usp$Var2 <- as.numeric(as.POSIXct(paste0(usp$Var2, "-01-01")))
p4 <- Rickshaw$new()
p4$layer(value ~ Var2, group = "Var1", data = usp, type = "area", width = 560)
# add a helpful slider this easily; other features TRUE as a default
p4$set(slider = TRUE)
p4$print("chart6")
p4



###########################################################################
## googleVis
###########################################################################

## Page 112 ---------------------------------------------------------------
library(googleVis)
library(rjson)
cat(toJSON(CityPopularity))  ## example data from googleVis


## Page 114 ---------------------------------------------------------------
## demo(WorldBank)
# inds <- c('SP.DYN.TFRT.IN','SP.DYN.LE00.IN', 'SP.POP.TOTL',
#           'NY.GDP.PCAP.CD', 'SE.ADT.1524.LT.FE.ZS')
# indnams <- c("fertility.rate", "life.expectancy", "population",
#             "GDP.per.capita.Current.USD", "15.to.25.yr.female.literacy")
# wdiData <- WDI(country="all", indicator=inds,
#                 start=1960, end=format(Sys.Date(), "%Y"), extra=TRUE)
# colnum <- match(inds, names(wdiData))
# names(wdiData)[colnum] <- indnams
# 
## Create a motion chart
library(googleVis)
# WorldBank <- droplevels(subset(wdiData, !region %in% "Aggregates"))
## save(inds, indnams, wdiData, colnum, WorldBank, file=paste(ojtpath, "R/Data/WorldBank.RData", sep=""))
load(paste(ojtpath, "R/Data/WorldBank.RData", sep=""))
M <- gvisMotionChart(WorldBank,
                     idvar="country", timevar="year",
                     xvar="life.expectancy", yvar="fertility.rate",
                     colorvar="region", sizevar="population",
                     options=list(width=700, height=600))
## Display the chart in the browser
plot(M)
##

## Page 115 ---------------------------------------------------------------
Motion=gvisMotionChart(Fruits, idvar="Fruit", timevar="Year", options=list(height=350, width=400))
# Display chart
plot(Motion)
# Create Google Gadget
cat(createGoogleGadget(Motion), file="motionchart.xml")


## Page 116 ---------------------------------------------------------------
Geo=gvisGeoMap(Exports, locationvar="Country", numvar="Profit", options=list(dataMode='regions'))
# Display chart
plot(Geo)
# Create Google Gadget
cat(createGoogleGadget(Geo), file="geomap.xml")


## Page 117 ---------------------------------------------------------------
AndrewGeo <- gvisGeoMap(Andrew, locationvar="LatLong", numvar="Speed_kt", hovervar="Category",
                        options=list(height=250, width=400, region="US", dataMode='markers'))
# Display chart
plot(AndrewGeo)
# Create Google Gadget


## Page 118 ---------------------------------------------------------------
AndrewMap <- gvisMap(Andrew, "LatLong" , "Tip",
                     options=list(showTip=TRUE, showLine=TRUE, enableScrollWheel=TRUE,
                                  mapType='terrain', useMapTypeControl=TRUE))
# Display chart
plot(AndrewMap)
# Create Google Gadget
cat(createGoogleGadget(AndrewMap), file="andrewmap.xml")


## Page 119 ---------------------------------------------------------------
GeoChart <- gvisGeoChart(Exports, "Country", "Profit",
                         options=list(region="150"))
# Display chart
plot(GeoChart)
# Create Google Gadget
cat(createGoogleGadget(GeoChart), file="geochart.xml")


## Page 120 ---------------------------------------------------------------
Table <- gvisTable(Exports, options=list(width=400, height=270))
# Display chart
plot(Table)
# Create Google Gadget
cat(createGoogleGadget(Table), file="table.xml")


## Page 121 ---------------------------------------------------------------
PopTable <- gvisTable(Population, options=list(width=600, height=300, page='enable'))
# Display chart
plot(PopTable)
# Create Google Gadget
cat(createGoogleGadget(PopTable), file="poptable.xml")


## Page 123 - 124 ----------------------------------------------------------
require(datasets)
states <- data.frame(state.name, state.area)
states3 <- data.frame(state.region, state.division, state.name, state.area)
regions <- aggregate(list(region.area=states3$state.area), list(region=state.region), sum)
 
divisions <- aggregate(list(division.area=states3$state.area),
                       list(division=state.division, region=state.region),
                       sum)
my.states3 <- data.frame(regionid=c("USA",
                                    as.character(regions$region),
                                    as.character(divisions$division),
                                    as.character(states3$state.name)),
                         parentid=c(NA, rep("USA", 4),
                                    as.character(divisions$region),
                                    as.character(states3$state.division)),
                         state.area=c(sum(states3$state.area),
                                      regions$region.area,
                                      divisions$division.area,
                                      states3$state.area))

my.states3$state.area.log=log(my.states3$state.area)

statesTree3 <- gvisTreeMap(my.states3, "regionid", "parentid",
                           "state.area", "state.area.log",
                           options=list(showScale=TRUE, width=400, height=300))
# Display chart
plot(statesTree3)
# Create Google Gadget
cat(createGoogleGadget(statesTree3), file="statestreemap.xml")


## Page 125 ---------------------------------------------------------------
df <- data.frame(country=c("US", "GB", "BR"), val1=c(1,3,4), val2=c(23,12,32))
Line <- gvisLineChart(df, options=list(legend='none', width=300, height=200))
plot(Line)
cat(createGoogleGadget(Line), file="linechart.xml")


## Page 126 ---------------------------------------------------------------
Bar <- gvisBarChart(df, options=list(legend='none', width=300, height=200))
plot(Bar)
cat(createGoogleGadget(Bar), file="barchart.xml")


## Page 127 ---------------------------------------------------------------
Area <- gvisAreaChart(df, options=list(legend='none', width=300, height=300))
plot(Area)
cat(createGoogleGadget(Area), file="areachart.xml")


## Page 128 ---------------------------------------------------------------
Bubble <- gvisBubbleChart(Fruits, idvar="Fruit", xvar="Sales", yvar="Expenses",
                          colorvar="Year", sizevar="Profit",
                          options=list(hAxis='{minValue:75, maxValue:125}',
                                       width=500, height=300))
plot(Bubble)
cat(createGoogleGadget(Bubble), file="bubblechart.xml")


## Page 129 ---------------------------------------------------------------
Pie <- gvisPieChart(CityPopularity, options=list(width=400, height=200))
plot(Pie)
cat(createGoogleGadget(Pie), file="piechart.xml")



###########################################################################
## shiny
###########################################################################
library(shiny)


## Page 129 ---------------------------------------------------------------
runExample("01_hello")



###########################################################################
## plotGoogleMaps
###########################################################################

## Page 144 ---------------------------------------------------------------
if(.Platform$OS.type == "windows"){
  sdata <- paste(ojtpath, "R/Data/sbus_info_win.csv", sep="")
} else{
  sdata <- paste(ojtpath, "R/Data/sbus_info.csv", sep="")
}
tmp <- read.csv(sdata, header=TRUE, colClasses = "character")
head(tmp)
table(tmp[,"typeName"])
busType <- sort(unique(tmp[,"typeName"]))
library(rgdal)
library(scales)
ind <- c(2, 4, 1, 5, 3)
ind <- ind[match(busType,  c("간선", "공항", "광역", "마을", "지선"))]
cols <- brewer.pal(8, "Set1")[ind]
colsa <- alpha(cols, 0.4)



## Page 145 ---------------------------------------------------------------
plot(c(126.4, 127.2), c(37.3, 37.8), type="n", xlab="Longitude", ylab="Latitude", asp=1)
for(i in 1:nrow(tmp)){
  lane <- read.csv(paste(ojtpath, "R/Data/sbus/", tmp[i,"id"], "_path.csv", sep=""),
                   header=TRUE)
  lane <- as.data.frame(lane)
  coordinates(lane) <- c("x", "y")
  proj4string(lane) <- CRS("+proj=tmerc +lat_0=38 +lon_0=127.5 +k=0.9996 
                           +k=1 +x_0=1000000 +y_0=2000000 +ellps=GRS80 
                           +units=m +no_defs")
  lane <- spTransform(lane, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 
                                +units=m +no_defs"))
  lines(coordinates(lane), col = colsa[busType == tmp[i,"typeName"]], lwd=2)
}
legend("topleft", legend=busType, col=cols, lwd=4, lty=1)
title("20% of bus routes in Seoul", cex.main=3)



## Page 147 ---------------------------------------------------------------
lanefull <- NULL
for(iter in 1:nrow(tmp)){
  iter2 <- ifelse(iter < 10, paste("0000", iter, sep=""), 
                  ifelse(iter < 100, paste("000", iter, sep=""), 
                         ifelse(iter < 1000, paste("00", iter, sep=""), 
                                ifelse(iter < 10000, paste("0", iter, sep=""), iter))))
  lane <- read.csv(paste(ojtpath, "R/Data/sbus/", tmp[iter,"id"], "_path.csv",sep=""), header=TRUE)
  lane <- as.data.frame(lane)
  lanefull <- rbind(lanefull, cbind(iter2, lane))   
}
library(plotGoogleMaps)
library(sp)
uid <- unique(lanefull[,1])
n <- length(uid)
lane2 <- list()
for(iter in 1:n){
  lane <- subset(lanefull, lanefull[,1] == uid[iter])[,2:3]
  lane <- apply(lane, 2, as.numeric)
  lane <- Line(lane)
  lane <- Lines(lane, ID=uid[iter])
  lane2[[iter]] <- lane
}


## Page 148 ---------------------------------------------------------------
lane3 <- SpatialLines(lane2, proj4string =  CRS("+proj=tmerc +lat_0=38 +lon_0=127.5 +k=0.9996 +k=1
                                                +x_0=1000000 +y_0=2000000 +ellps=GRS80 +units=m
                                                +no_defs"))
dfv <- as.data.frame(tmp)
rownames(dfv) <- uid
lane4 <- SpatialLinesDataFrame(lane3, dfv, match.ID = TRUE)

#setwd( "/home/jang/")
cols <- brewer.pal(8, "Set1")[c(5, 4, 2, 3, 1)]
mapMeuseCl<- plotGoogleMaps(lane4, filename='myMap2.htm', zcol=4, colPalette=cols,
                           strokeOpacity = .7,
                           strokeWeight = 2,
                           layerName="Bus Type",
                           mapTypeId='ROADMAP')


## Page 150 ---------------------------------------------------------------
library(shapefiles)
write.shapefile2 <- function (shapefile, out.name, arcgis = FALSE, max_nchar=6553){
  write.shp(shapefile$shp, paste(out.name, ".shp", sep = ""))
  write.shx(shapefile$shx, paste(out.name, ".shx", sep = ""))
  write.dbf2(shapefile$dbf, paste(out.name, ".dbf", sep = ""), arcgis, max_nchar)
}

write.dbf2<- function (dbf, out.name, arcgis = FALSE, max_nchar = 6553) {
  library(foreign)
  if (arcgis == T) {
    colnames(dbf$dbf) <- gsub("\\.", "_", colnames(dbf$dbf))
  }
  get("write.dbf", "package:foreign")(dbf$dbf, out.name, max_nchar =6553)
}


## Page 151 ---------------------------------------------------------------
dd <- data.frame(Obs=lanefull[,1], X=as.numeric(lanefull[,2]), Y=as.numeric(lanefull[,3]))
outencoding <- "UTF-8"
if(.Platform$OS.type == "windows") outencoding <- "CP949"
if(.Platform$OS.type == "windows" | outencoding != "CP949"){
  ddTable <- data.frame(Obs=unique(lanefull[,1]),  tmp[1:length(unique(lanefull[,1])), ])
} else { 
  ddTable <- data.frame(Obs=unique(lanefull[,1]), 
                          iconv(as.matrix(tmp[1:length(unique(lanefull[,1])), ]), "UTF-8", "CP949"))
}


ddShapefile <- convert.to.shapefile(dd, ddTable, "Obs", 3)
write.shapefile2(ddShapefile, paste(ojtpath, "R/shp/busroute", sep=""), arcgis=TRUE)


###########################################################################
## leafletR
###########################################################################

## Page 156 ---------------------------------------------------------------
library(leafletR)
data(quakes)
# store data in GeoJSON file (just a subset here)
q.dat <- toGeoJSON(data=quakes[1:99,], dest=tempdir(), name="quakes")
 
# File saved under /tmp/RtmpBcMTre/quakes.geojson
# make style based on quake magnitude
q.style <- styleGrad(prop="mag", breaks=seq(4, 6.5, by=0.5),
                     style.val=rev(heat.colors(5)), leg="Richter Magnitude",
                     fill.alpha=0.7, rad=8)
# create map
q.map <- leaflet(data=q.dat, dest=tempdir(), title="Fiji Earthquakes",
                 base.map="mqsat", style=q.style, popup="mag")
 
# Your leaflet map has been saved under /tmp/RtmpBcMTre/Fiji_Earthquakes/Fiji_Earthquakes.html
browseURL(q.map)


## Page 158 - 160 ----------------------------------------------------------
library(rgdal) #for reading/writing geo files
library(rgeos) #for simplification
library(sp)
# Get the data note that this file is somewhat big so it might take a couple
# of minutes to download
url <- "http://www2.census.gov/geo/tiger/TIGER2010DP1/County_2010Census_DP1.zip"
downloaddir <- paste(ojtpath, "R/Data/tiger", sep="")

destname <- "tiger.zip"
download.file(url, destname)
unzip(destname, exdir=downloaddir, junkpaths=TRUE)

# A little clipping and projecting
filename <- list.files(downloaddir, pattern=".shp", full.names=FALSE)
filename <- gsub(".shp", "", filename)

# ----- Read in shapefile (NAD83 coordinate system)
# ----- this is a fairly big shapefile and takes 1 minute to read
dat <- readOGR(downloaddir, filename) 


# ----- Create a subset of New York counties
subdat <- dat[substring(dat$GEOID10, 1, 2) == "36",]
# ----- Transform to EPSG 4326 - WGS84 (required)
subdat <- spTransform(subdat, CRS("+init=epsg:4326"))
# ----- change name of field we will map
names(subdat)[names(subdat) == "DP0010001"]<-"Population"

#Simplify your shapefile if necessary

# ----- save the data slot
subdat_data <- subdat@data[, c("GEOID10", "Population")]
# ----- simplification yields a SpatialPolygons class
subdat <- gSimplify(subdat, tol=0.01, topologyPreserve=TRUE)
# ----- to write to geojson we need a SpatialPolygonsDataFrame
subdat <- SpatialPolygonsDataFrame(subdat, data=subdat_data)

# Ready to play with LeafletR
# ----- Write data to GeoJSON
leafdat <- paste(downloaddir, "/", filename, ".geojson", sep="") 
writeOGR(subdat, leafdat, layer="", driver="GeoJSON")

# ----- Create the cuts
cuts <- round(quantile(subdat$Population, probs = seq(0, 1, 0.20), na.rm = FALSE), 0)
cuts[1] <- 0 # ----- for this example make first cut zero

# ----- Fields to include in the popup
popup<-c("GEOID10", "Population")

# ----- Gradulated style based on an attribute
sty <- styleGrad(prop="Population", breaks=cuts, right=FALSE, style.par="col",
                 style.val=rev(heat.colors(6)), leg="Population (2010)", lwd=1)

# ----- Create the map and load into browser
map <- leaflet(data=leafdat, dest=downloaddir, style=sty,
               title="index", base.map="osm",
               incl.data=TRUE,  popup=popup)

# Your leaflet map has been saved under d:/Leaflet/index/index.html
# ----- to look at the map you can use this code
browseURL(map)