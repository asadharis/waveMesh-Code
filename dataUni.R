# In this file we consider the data analysis by only considering the 
# motorcycle dataset. 
# 
# three features which don't allow normal methods being used. 
# 1. Sample size not a power of 2
# 2. unequally distributed covariates
# 3. repeated measurements. 
# Our method can handle all 3 of these issues. 

library(adlift)
library(wavethresh)

source("IRRwd.R")
source("waveMesh.R")

# Be begin with some functions for processing data. 
process.dat <- function(name = "mcycle") {
  if(name == "mcycle") {
    require(MASS)
    dat <- MASS::mcycle
  } else if(name == "RUAndromeda") {
    dat <- read.table(file = "RUAndromeda.txt", 
                      header = TRUE, sep = ",")
    dat <- dat[dat$Magnitude != "<11.6" &
                 dat$Magnitude != "<12.1" &
                 dat$Magnitude != "<12.6", 1:2]
    dat <- droplevels(dat)
    dat$Magnitude <- as.numeric(paste(dat$Magnitude))
  } else if(name == "ipd") {
    data(ipd, package = "wavethresh")
    dat <- ipd; 
    time <- attributes(dat)$tspar["start"]
    time <- time + 0.02 * 0:(length(dat)-1)
    dat <- data.frame("time" = time, "ipd" = as.numeric(dat) )
    
    # Return a 1000 random values of the data set.
    # This makes the sample size not a power of 2
    # makes the covariates irregularly spaced.
    #set.seed(0)
    dat <- dat[sample(1:nrow(dat), size = 1000),]
  }
  # Average values with repeated measurements.
  avg.vals <- tapply(dat[,2], dat[,1], FUN  = mean)
  time <- as.numeric(names(avg.vals))
  timeScale <- (time - min(time))/(diff(range(time)))
  dat2 <- data.frame("time" = time, "time2" = timeScale,
                     "y" = as.numeric(avg.vals))
  
  return(dat2)
}


require(adlift)

filter.number = 8; j0 = 2

dat1 <- process.dat("mcycle")


# We now begin with the data analysis.
# First we implement our proposal to the full
# training set to obtain a working lambda sequence, 
# and also to use later for predictions on test set.
modFull <- waveMesh(dat1$y, dat1$time2, nlam = 50,
                    lam.min.ratio = 1e-6, 
                    nlevel = ceiling(log2(nrow(dat1))),
                    max.iter = 500, tol = 1e-4,
                    theta0 = NULL, filter.number = filter.number, 
                    j0 = j0,
                    family = "DaubExPhase", bc = "symmetric")
lams <- modFull$lams

# CV for our proposal
cvWM <- cv.waveMesh(dat1$y, dat1$time2, K = 5, lams = lams,
                    nlevel = ceiling(log2( nrow(dat1) )), 
                    max.iter = 500, tol = 1e-4,
                    theta0 = NULL, filter.number = filter.number, 
                    j0 = j0,
                    family = "DaubExPhase", bc = "symmetric")
plot(dat1$time, dat1$y, 
     xlab = "Time (ms)", ylab = "Acceleration (g)", 
     cex.lab = 1.4, pch = 16, main = "waveMesh",
     cex.main = 2)
lines(dat1$time, modFull$ans[,cvWM$ind.1se], col = "red", lwd = 3)

IWD <- irr.wd(dat1$time2, dat1$y, filter.number = filter.number, 
              j0 = j0,
              family = "DaubExPhase", bc = "symmetric", 
              policy = "universal")

plot(dat1$time, dat1$y, 
     xlab = "Time (ms)", ylab = "Acceleration (g)", 
     cex.lab = 1.4, pch = 16, 
     main = "Interpolation", cex.main = 2)
lines(dat1$time, IWD, col = "blue", lwd = 3)


ISO <- irr.iso(dat1$time2, dat1$y, filter.number = filter.number, 
               j0 = j0, 
               family = "DaubExPhase", bc = "symmetric", 
               policy = "universal")

plot(dat1$time, dat1$y, 
     xlab = "Time (ms)", ylab = "Acceleration (g)", 
     cex.lab = 1.4, pch = 16, 
     main = "Isometric", cex.main = 2)
lines(dat1$time, ISO, col = "purple", lwd = 3)

AdL <- adlift::denoise(dat1$time2, dat1$y, AdaptNeigh, 
                       1, TRUE, TRUE, 2, rule = "soft")
plot(dat1$time, dat1$y, 
     xlab = "Time (ms)", ylab = "Acceleration (g)", 
     cex.lab = 1.4, pch = 16, 
     main = "Adaptive Lifting", cex.main = 2)
lines(dat1$time, AdL, col = "cyan", lwd = 3)



