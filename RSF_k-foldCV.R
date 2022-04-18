# Date: April 18, 2022
# Author: Eric Margenau
# This script is designed to assess the predictive ability of fitted models using methods described by Boyce et al. 2002 (Ecol. Model.).
# For each filtered data set ('subset' function), data are grouped into 1 of 5 folds. All observations from a single individual are
# kept within a single fold. The 5-fold procedure was replicated 5 times. Pearson's correlation coefficient was used to assess the 
# predictive ability of each model, with higher correlations indicating better predictive ability. 

list.of.packages <- c('caret', 'plyr', 'ggplot2', 'glmmTMB', 'groupdata2', 'gstat', 'sp')
install <- list.of.packages %in% installed.packages()
if(length(list.of.packages[!install]) > 0) install.packages(list.of.packages[!install])
lapply(list.of.packages, require, character.only = TRUE)

# Read in moose data
setwd("C:/Users/EricMargenau/OneDrive - USDA/Moose Project/Publications/Resource-selection")
d <- read.csv('moose_final_100perkm2.csv')

cols <- c('sex', 'season', 'landfire', 'mooseid')
d[cols] <- lapply(d[cols], factor)  # Change columns to factors

summary(d)
head(d)

# Subset whatever sex-season combination needed
sex.sub <- subset(d, sex == 'F' & season == 'autumn' & landfire != 'Other' & landfire != 'Agri ' & landfire != 'Water')
summary(sex.sub$landfire)

### K-fold cross validation (following Boyce 2002 & Lehman et al. 2016) ###

# Replicate k-fold cross validation procedure 5 times
kfolds1 <- groupdata2::fold(sex.sub, k = 5, id_col = 'mooseid', method = 'n_last')  # Create folds (folds are based on mooseid
# so that all locations from a given individual is within a single fold, which is why folds are not equal size)
summary(kfolds1$.folds)  # Check to make sure folds are similar sizes

kfolds2 <- groupdata2::fold(sex.sub, k = 5, id_col = 'mooseid', method = 'n_last')
summary(kfolds2$.folds)  # Check to make sure folds are similar sizes

kfolds3 <- groupdata2::fold(sex.sub, k = 5, id_col = 'mooseid', method = 'n_last')
summary(kfolds3$.folds)  # Check to make sure folds are similar sizes

kfolds4 <- groupdata2::fold(sex.sub, k = 5, id_col = 'mooseid', method = 'n_last')
summary(kfolds4$.folds)  # Check to make sure folds are similar sizes

kfolds5 <- groupdata2::fold(sex.sub, k = 5, id_col = 'mooseid', method = 'n_last')
summary(kfolds5$.folds)  # Check to make sure folds are similar sizes


# Empty data frames for loops to populate
out.df <- data.frame(Fold = as.character(), Bin = as.numeric(), AAFreq = as.numeric())
bin.df <- data.frame(Bin = as.numeric(), AAFreq = as.numeric())
fm.df <- data.frame(Fold = as.numeric(), Bin = as.numeric(), AAFreq = as.numeric(), Rep = as.numeric())

# Bins of equal probability width
bin <- data.frame(low = c(0, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91),
                  high = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

## Loop through k-folds(5) and bins(10 for each fold) at 5 replications each
for(r in 1:5){  # Loop through replications
  for(i in 1:5){  # Loop through folds
    training <- subset(get(paste('kfolds', r, sep = '')), .folds != i)  # Training data
    testing <- subset(get(paste('kfolds', r, sep = '')), .folds == i)  # Testing data
    new.lm <- glmmTMB(used ~ landfire * sapbiom + (1 | mooseid) + (1 | yr), 
                      data = training, family = 'binomial')
    new.pred <- predict(new.lm, newdata = testing, allow.new.levels = TRUE)  # Predict new values, keep on link scale as per Northrup et al. 2021
    testing$new_pred <- new.pred  # Add new predictions to testing data frame
    
    for(b in 1:10){  # Loop through bins within each fold
      
      used.pts <- subset(testing, new_pred >= quantile(new_pred, probs = c(bin$low[b], bin$high[b]))[1] & 
                           new_pred <= quantile(new_pred, probs = c(bin$low[b], bin$high[b]))[2] & used == 1)  # Used points in bin
      avail.pts <- subset(testing, new_pred >= quantile(new_pred, probs = c(bin$low[b], bin$high[b]))[1] &
                            new_pred <= quantile(new_pred, probs = c(bin$low[b], bin$high[b]))[2] & used == 0)  # Available points in bin
      
      tab.pts <- data.frame(cbind(summary(used.pts$landfire), summary(avail.pts$landfire)))  # Create data frame with used/available points 
      # for calculating area-adjusted frequency
      tab.rm0 <- data.frame(tab.pts[rowSums(tab.pts[]) > 0, ])  # Removes excess rows with '0' data
      tab.rm0$X2 <- ifelse(tab.rm0$X2 == 0, tab.rm0$X2 + 1, tab.rm0$X2)
      tab.rm0$aafreq <- tab.rm0$X1 / tab.rm0$X2  # Calculate area-adjusted frequency for each variable
      aa.freq <- mean(tab.rm0$aafreq)  # Calculates area-adjusted frequency for the bin
      
      bin.tmp <- data.frame(Bin = b, AAFreq = aa.freq)
      
      bin.df <- rbind(bin.df, bin.tmp)  # Having trouble clearing bin.df after all bins within each fold loop are run, 
      # creates repeats in out.df that need to be removed
      
    }
    
    tmp.df <- data.frame(Fold = i, bin.df)
    out.df <- rbind(out.df, tmp.df)
    
    rm(tmp.df)
    
  }
  
  fm.tmp <- out.df[-c(11:20, 31:50, 61:90, 101:140), ]
  fm.tmp$Rep <- r  # Add replication number to data
  
  fm.df <- rbind(fm.df, fm.tmp)
  
  rm(fm.tmp)
  out.df <- out.df[0, ]  # Reset data frame
  bin.df <- bin.df[0, ]  # Reset data frame
  
}

fm1.results <- fm.df

comb.df <- data.frame(cbind(aggregate(fm1.results$AAFreq, by = list(fm1.results$Bin), FUN = mean), 
                            aggregate(fm1.results$AAFreq, by = list(fm1.results$Bin), FUN = sd)))
cor.test(comb.df$Group.1, comb.df$x, method = 'spearman')

new.df$Bin <- as.factor(new.df$Bin)  # Make Bin and Fold factors for graphing purposes
new.df$Fold <- as.factor(new.df$Fold)

# Graph area-adjusted frequencies
ggplot(comb.df, aes(Group.1, x)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = x - x.1, ymax = x + x.1), width = 0.2) +
  theme_classic() 