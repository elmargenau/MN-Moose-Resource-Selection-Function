# Date: April 18, 2022
# Author: Eric Margenau
# This code is designed to read in used and available spatial points and extract data from raster layers. Raster layers are 
# covariate data to be used in statistical analysis. I first read in the covariate rasters and then read in used and available point 
# data and then extract covariate information associated with used and available data. Covariate data are then transformed to remove
# issues with magnitude in values and censored moose are removed prior to writing a .csv file to be ready for RSF analysis. 

# Packages needed 
list.of.packages <- c('dismo', 'rgdal', 'raster', 'rgeos', 'spatstat', 'maptools')
install <- list.of.packages %in% installed.packages()
if(length(list.of.packages[!install]) > 0) install.packages(list.of.packages[!install])
lapply(list.of.packages, require, character.only = TRUE)

# Read in covariate data from this folder
setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection/Shapefiles")

## Raster layers
landcov.typ <- raster('LandCov.tif') # Land cover
recl.tab1 <- read.csv('Landfire_convert.csv') # Reclass land cover class into categories
landcov.mat <- as.matrix(recl.tab1[, c(1:2)])
landcov.recl <- raster::reclassify(landcov.typ, landcov.mat)

landcov.age <- raster('Age.tif') # Age 
recl.tab2 <- as.matrix(read.csv('Reclass_age.csv')) # Reclass age into 5 categories
age.recl <- reclassify(landcov.age, recl.tab2)

sapbiom <- raster('ForageBiomassLANDIS.tif') # Sapling biomass
elev <- raster('Elev.tif') # Elevation
slopeasp <- raster('SlopeAsp.tif') # Slope aspect
slopeperc <- raster('SlopePerc.tif') # Slope percent

# Change working directory to the correct folder, which each used:available ratio has that contains their shapefiles and stuff
setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection/MCP_HomeRanges/MCP_homerange_1to3")
#setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection/MCP_HomeRanges/MCP_homerange_1to1")
#setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection/MCP_HomeRanges/MCP_homerange_1to5")
#setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection/MCP_HomeRanges/MCP_100ptpersqkm")

# Read in xy data
xy.fem <- read.csv('moose_rsf_seasonal_females.csv')
xy.fem$sex <- 'F'
xy.mal <- read.csv('moose_rsf_seasonal_males.csv')
xy.mal$sex <- 'M'
xy.tmp <- data.frame(rbind(xy.fem, xy.mal))
xy.tmp <- xy.tmp[-c(3515961, 5244458), ] # Removed the same rows as 'xy' to match

# Format xy to create point shapefile to talk w/ raster
data.tmp <- cbind(xy.tmp$x, xy.tmp$y)
xy <- matrix(data.tmp, ncol = 2)
which(is.na(xy))
xy <- xy[-c(3515961, 5244458), ]
which(is.na(xy)) # Check to make sure there are no more NAs
xysp <- SpatialPoints(xy, proj4string = CRS("+proj=utm +zone=15 +datum=WGS84"))  ## this projection is the one you also want to use 

## Extract information from raster layers at each point (used and available)
landcov.pt <- raster::extract(landcov.recl, xysp)
landcovage.pt <- raster::extract(age.recl, xysp)
sapbiom.pt <- raster::extract(sapbiom, xysp)
elev.pt <- raster::extract(elev, xysp)
slopeasp.pt <- raster::extract(slopeasp, xysp)
slopeperc.pt <- raster::extract(slopeperc, xysp)

# Combine raster extractions with moose data
# Important to manually check to make sure xy coordinates are aligning with correct values from raster layers
dat.tmp <- data.frame(mooseid = xy.tmp$mooseid, sex = xy.tmp$sex, age = xy.tmp$age, season = xy.tmp$season, yr = xy.tmp$yr,
                      x = xy.tmp$x, y = xy.tmp$y, used = xy.tmp$used, landcovtyp = landcov.pt, landcovage = landcovage.pt,
                      forbiom = sapbiom.pt, elev = elev.pt, aspect = slopeasp.pt, percent = slopeperc.pt)

# Remove moose that are censured from analysis (could do sooner if wanted)
dat.tmp1 <- subset(dat.tmp, mooseid != 71 & mooseid != 17 & mooseid != 34 & mooseid != 60 &
                     mooseid != 102 & mooseid != 119 & mooseid != 153 & mooseid != 154 &
                     mooseid != 168 & mooseid != 170 & mooseid != 185 & mooseid != 193 &
                     mooseid != 213)
summary(dat.tmp1)

# Remove NAs from data frame (this may change depending on the number of available points)
dat.tmp1$season <- as.factor(dat.tmp1$season)
d.na1 <- dat.tmp1[!is.na(dat.tmp1$season), ]
summary(d.na1)
d.na2 <- d.na1[!is.na(d.na1$landcovage), ]
summary(d.na2)
d.na3 <- d.na2[!is.na(d.na2$forbiom), ]
summary(d.na3)

# Transformations, designed to get variables to similar scales
d.na3$forbiom <- scale(d.na3$forbiom)
d.na3$aspect <- (cos(45 - d.na3$aspect) + 1) * 2 + 1  # modified Beers transformation
d.na3$elev <- scale(d.na3$elev)
d.na3$percent <- scale(d.na3$percent)
d.na3$landcovtyp <- as.factor(d.na3$landcovtyp)  # Important for joining
summary(d.na3)

# Read in another file that has the land cover names (using names is easier than anID number)
setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection/Shapefiles")
d.1 <- read.csv('Reclass_landcov.csv')
d.1$landfire <- as.factor(d.1$landfire)
d.1$landcovtyp <- as.factor(d.1$landcovtyp)  # Important for joining
head(d.1)

# Join data frames
d.2 <- plyr::join(d.na3, d.1, by = 'landcovtyp', type = 'left')
df.moose <- d.2[!is.na(d.2$landfire), ]  # Remove any remaining NAs in data frame
which(is.na(df.moose$landfire))
head(df.moose)
df.export <- df.moose[, -16]  # Remove unnecessary columns 
head(df.export)

setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection")
write.csv(df.export, file = 'moose_final_1to3.csv')
