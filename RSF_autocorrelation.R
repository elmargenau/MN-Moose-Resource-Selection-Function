# Date: April 18, 2022
# Author: Eric Margenau
# This code is designed to assess the spatial correlation of used and available points for each individual in each season and each year.

list.of.packages <- c('caret', 'plyr', 'ggplot2', 'glmmTMB', 'groupdata2', 'gstat', 'sp')
install <- list.of.packages %in% installed.packages()
if(length(list.of.packages[!install]) > 0) install.packages(list.of.packages[!install])
lapply(list.of.packages, require, character.only = TRUE)

# Read in moose data
setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection")
d <- read.csv('moose_final_1to3.csv')

cols <- c('sex', 'season', 'landfire', 'mooseid')
d[cols] <- lapply(d[cols], factor)  # Change columns to factors

summary(d)
head(d)

d$mooseid <- sprintf("%03d", as.numeric(as.character(d$mooseid)))  # Make sure all mooseids have three significant digits
d$mooseid <- as.character(d$mooseid)  # Format to character for reading in for-loop

sub.d <- subset(d, used == 0)  # Subset 'used' to either used (1), available (0), or both (not necessary to subset)
summary(sub.d)

mooselist <- unique(sub.0$mooseid) # List all individual animals
nmoose <- length(mooselist)

results.df <- data.frame(mooseid = as.character(), season = as.numeric(), year = as.numeric(), 
                         dist = as.numeric(), gamma = as.numeric())

for (i in 1:nmoose){ # Start moose loop
  
  m.tmp <- sub.0[which(sub.0$mooseid == unique(sub.0$mooseid)[i]), ] # Pulls data for a single moose
  
  seas <- unique(m.tmp$season) # How many seasons for that moose
  yr <- unique(m.tmp$yr) # How many years for that moose
  
  for (z in 1:length(seas)){ # Start season loop
    for (j in 1:length(yr)){ # Start year loop
      
      sub.moose <- m.tmp[m.tmp$season == seas[z] & m.tmp$yr == yr[j], ] 
      if (nrow(sub.moose) == 0) next
      
      d.semivar <- sub.moose[, c(7:15)]
      coordinates(d.semivar) <- ~x+y
      elev.var <- gstat::variogram(elev~1, d.semivar)
      
      tmp.df <- data.frame(mooseid = mooselist[i], season = seas[z], year = yr[j], 
                           dist = elev.var$dist, gamma = elev.var$gamma)
      results.df <- rbind(results.df, tmp.df)
      
    }
  }
}

write.csv(results.df, file = 'SpatialVario_avail.csv')
#write.csv(results.df, file = 'SpatialVario_used.csv')

# Graphing results at three different scales
plot(gamma ~ dist, results.df, xlim = c(0, 10000), ylim = c(0, 3))
plot(gamma ~ dist, results.df, xlim = c(0, 2000), ylim = c(0, 2))
plot(gamma ~ dist, results.df, xlim = c(0, 500), ylim = c(0, 0.5))