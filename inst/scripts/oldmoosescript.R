library(FPBKPack2)
library(sp)
library(rgeos)
library(rgdal)

# use Alaska Moose data
data(AKmoose)
# plot data, which are in polygons in Alaska Albers projection
plot(AKmoose)

AKmoose@proj4string
# see http://spatialreference.org/ref/epsg/nad83-alaska-albers/

# transform the projection to lat/lon
latlon = spTransform(AKmoose, CRS("+init=epsg:4326"))
plot(latlon)
latlon@proj4string
# see http://spatialreference.org/ref/epsg/4326/

#get centroids of polygons
centroids = data.frame(ID = latlon@data,
  x = gCentroid(latlon,byid=TRUE)@coords[,'x'],
  y = gCentroid(latlon,byid=TRUE)@coords[,'y'])

# create user-defined UTM coordinates for spatial modeling
xy = LLtoUTM(mean(centroids$x), centroids$y, centroids$x)$xy

# pull the data.frame from the SpatialPolygonsDataFrame object
AKmoose@plotOrder
d1 = AKmoose@data
rownames(d1)
str(d1)
# add new UTM coordinates to data.frame,
d1$x = xy[,'x']
d1$y = xy[,'y']
# check to make sure that ordering didn't get messed up
# by comparing to lat lon in the data frame
d1[1:30, c('CENTRLAT', 'y', 'CENTRLON','x')]
d1$totalnum <- as.numeric(levels(d1$total))[d1$total]
d1$sampind <- as.numeric(levels(d1$surveyed))[d1$surveyed]
d1$totalnum[d1$sampind == 0] <- NA
# fit a model where STRAT is a covariate

detectiondat <- c(rep(1, 42), rep(0, 8))
detectiondf <- data.frame(detectiondat)
detmoose <- get_detection(formula = detectiondat ~ 1,
  data = detectiondf, varmethod = "Delta")

slmfit_out1 = slmfit(formula = totalnum ~ strat, data = d1, xcoordcol = 'x', ycoordcol = 'y',
  CorModel = "Exponential", detectionobj = detmoose)


# predictions
predout1 = predict(slmfit_out1)
predout1$FPBK_Prediction - 1.645 * sqrt(predout1$PredVar)
predout1$FPBK_Prediction + 1.645 * sqrt(predout1$PredVar)

#multiobj <- multistrat(formula = totalnum ~ 1,
#  data = d1, xcoordcol = 'x', ycoordcol = 'y',
#  CorModel = "Exponential", detectionobj = detmoose,
#  stratcol = "strat")


slmfit_out2 = slmfit(formula = totalnum ~ strat, data = d1, xcoordcol = 'x', ycoordcol = 'y',
  CorModel = "Exponential")


# predictions
predout2 = predict(slmfit_out2, detinfo = c(0.84, 0.05))
predout2$FPBK_Prediction

