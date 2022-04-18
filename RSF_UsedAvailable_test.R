# Date: April 18, 2022
# Author: Eric Margenau
# This code is designed to test the different levels of available points and assess their accuracy in representing the true availability
# within an individual's seasonal home range. For each individual in each season and year, I first calculated the true proportion of land 
# cover within their home range ('prop.cov'). I then summarized the available points for each land cover type ('pc.df') for that individual,
# season, year combination. I then calculated Pearson's correlation coefficient between prop.cov and pc.df. I assumed that the greater the 
# correlation between the true proportion of land cover and number of available points within each land cover would indicate a level of 
# available points that accurately captures what is available to an individual moose. 

list.of.packages <- c('plyr', 'maptools', 'raster', 'landscapemetrics', 'reshape', 'ggplot2')
install <- list.of.packages %in% installed.packages()
if(length(list.of.packages[!install]) > 0) install.packages(list.of.packages[!install])
lapply(list.of.packages, require, character.only = TRUE)

# Read in moose data
setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection")
d <- read.csv('moose_final_100perkm2.csv')  # or one of the below options depending on what level of available points to test
#d <- read.csv('moose_final_1to1.csv')
#d <- read.csv('moose_final_1to3.csv')
#d <- read.csv('moose_final_1to5.csv')

cols <- c('sex', 'season', 'landfire', 'mooseid')
d[cols] <- lapply(d[cols], factor)  # Change columns to factors

summary(d)
head(d)

d$mooseid <- sprintf("%03d", as.numeric(as.character(d$mooseid)))  # Make sure all mooseids have three significant digits
d$mooseid <- as.character(d$mooseid) # Format to character for reading in for-loop

sub.0 <- subset(d, used == 0)  # Subset available points
summary(sub.0)

mooselist <- unique(sub.0$mooseid) # List all individual animals
nmoose <- length(mooselist)

setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection/Shapefiles")
cover.type <- read.csv('Landfire_convert.csv')  # .csv file with land cover names to be joined with clipped rasters in for-loop

results.df <- data.frame(mooseid = as.character(), season = as.numeric(), year = as.numeric(), cor = as.numeric())  # Empty data frame to populated with results

## Do for-loop for each used:available ratio that needs to be tested

for (i in 1:nmoose){ # Start loop through all moose
  
  m.tmp <- sub.0[which(sub.0$mooseid == unique(sub.0$mooseid)[i]), ] # Pulls data for a single moose
  
  seas <- unique(m.tmp$season) # How many seasons for that moose
  yr <- unique(m.tmp$yr) # How many years for that moose
  
  for (z in 1:length(seas)){ # Start season loop
    for (j in 1:length(yr)){ # Start year loop
      
      sub.moose <- m.tmp[m.tmp$season == seas[z] & m.tmp$yr == yr[j], ] 
      
      # Load MCP home range shapefiles
      setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection/MCP_homerange_1to1/Seasonal/Shapefiles")
      shapefile <- readShapePoly(paste('MCP_', mooselist[i], '_', m.tmp$yr[j], '_', m.tmp$seas[z], '.shp', sep=""),
                                 proj4string = CRS("+proj=utm +zone=15 +datum=WGS84"))
      tmp.df <- data.frame(PlotNum = unique(shapefile[["SP_ID"]]), Area = as.numeric(NA))  # Plot label
      
      # Load rasters 
      setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection/Shapefiles")
      raster <- raster::raster('LandCov.tif')  # Load raster of land cover type across entire study region
      
      # For loop to clip the raster based on the home range shapefile
      clip_raster = function(spdata, codename, raster, dsn = getwd()){
        CODES <- unique(spdata[[codename]])
        for (i in 1:length(CODES)) {
          tmp <- spdata[spdata[[codename]] == CODES[i], ]
          raster_tmp <- crop(raster, tmp)
          raster_clip <- mask(raster_tmp, tmp)
          r_name <- paste("Plot_", CODES[i], ".tif", sep = "")
          writeRaster(raster_clip, r_name, overwrite = TRUE)
          plot(raster_clip)
        }
      }
      
      clip_raster(shapefile, "SP_ID", raster)  # Clipped raster of land cover type based on MCP home range
      
      t.rast <- raster('Plot_a.tif')  # Save raster clip as raster for further analysis
      prop.cov <- as.data.frame(lsm_c_pland(t.rast))  # Calculate the proportion of land cover for each home range
      pc.join <- plyr::join(prop.cov, cover.type, by = 'class', match = 'first')  # Join 'prop.cov' with land cover names file (cover.type)
      pc.agg <- aggregate(pc.join$value, by = list(pc.join$Reclass2), FUN = sum)  # Sum proportion for each land cover type
      
      pc.df <- data.frame(summary(sub.moose$landfire))  # Get data frame of within for-loop to compare with clipped raster output
      pc.df$Group.1 <- rownames(pc.df)
      
      join.pc <- plyr::join(pc.agg, pc.df, by = 'Group.1', match = 'first')  # Join clipped raster output with data set with available points
      join.pc[is.na(join.pc)] <- 0
      colnames(join.pc) <- c('landfire', 'propcov', 'available')
      
      cor.df <- cor(join.pc$propcov, join.pc$available, method = 'pearson')  # Pearson's correlation coefficient between clipped raster and available points
      
      tmp.df <- data.frame(mooseid = mooselist[i], season = seas[z], year = yr[j], cor = cor.df)
      results.df <- rbind(results.df, tmp.df)
      
      rm(t.rast)
      rm(prov.cov)
      rm(prop.join)
      rm(pc.agg)
      rm(pc.df)
      rm(join.pc)
      
    }
  }
}

write.csv(results.df, file = 'Correlation_test_100perkm2.csv')  # or one of the below files depending on the level of available points
#write.csv(results.df, file = 'Correlation_test_1to1.csv')
#write.csv(results.df, file = 'Correlation_test_1to3.csv')
#write.csv(results.df, file = 'Correlation_test_1to5.csv')