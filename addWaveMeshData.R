source("addWaveMesh.R")
source("AMlet.r")
source("wavethresh.r")
# This file consists of the main data analysis for the aditive models. 
# We implement our method we select the lambda value using CV. 

dataAdditive <- function(seed = 1) {
  
  require(wavethresh)
  require(MASS)
  # We begin with loading and processing the dataset.
  data("Boston")
  data <- Boston
  y <- data[, "medv"]
  covars <- data[, c("crim", "indus", "nox", "rm", "age",
                   "dis", "tax", "ptratio", "black", "lstat")]
  # Scale the covariates to lie within the unit interval.
  covars2 <- apply(covars, 2, function(vec) {
    (vec - min(vec))/(diff(range(vec)))
    })
  
  set.seed(seed)
  # Split the data into train/test
  n <- nrow(data)
  n.tr <- 2^8
  train <- sample(1:n, size = n.tr)
  
  # First we implement out proposal via cross validation.
  lam.seq <- 10^(seq(log10(30), log10(0.1), length = 50))
  
  cv.wm <- cv.addWaveMesh(y[train], covars2[train,], 
                          K = 5, max.iter.outer = 1000, 
                          tol.outer = 1e-4, max.iter = 1000, tol = 1e-4,
                          lam.seq = lam.seq, gamma = 0,
                          family = "DaubLeAsymm" , 
                          filter.number=8, bc = "periodic",
                          j0 = 2, nlevel = ceiling(log2(n.tr)))
  
  fit.wm <- fit.additive(y[train], covars2[train,], 
                         max.iter.outer= 1000, tol.outer = 1e-4, 
                         lam.seq = c(cv.wm$lam.min, cv.wm$lam.1se),
                         gamma = 0, family = "DaubLeAsymm" , 
                         filter.number=8, bc = "periodic",
                         j0 = 2, nlevel =  ceiling(log2(n.tr)), 
                         max.iter = 1000, tol = 1e-4)
  
  fhat.wm <- predict.add_mod(fit.wm, new.data = covars2[-train,],
                             type = "response")
  
  mse <- colMeans( (fhat.wm - y[-train] )^2 )

  
  # Followed by Amlet
  fit.amlet1 <- AMlet(y[train], covars2[train,], thresh.value="universal", 
                      sigma=NA, p=1, low.level=0, by.level=F, 
                      filter.number=8, family="DaubLeAsymm",
                      bc="periodic", 
                      conv.thresh=1.e-10, max.iter=ncol(covars2)*1000)
  
  fhat.amlet <- predict.amlet(fit.amlet1, covars2[train,],covars2[-train,])
  mse.amlet1 <- mean((fhat.amlet  - y[-train] )^2)
  
  res <- c("mseMin" = mse[1], "mse1se" = mse[2], "mseAM" = mse.amlet1)
  
  dir <- paste0("simAdditiveData")
  filename <- paste0(dir, "/seed", seed, ".RData")
  
  if(dir.exists(dir)) {
    save(res, file = filename)
  } else {
    dir.create(dir)
    save(res, file = filename)
  }
}

# A helper function for Amlet based on linear interpolation.
predict.amlet <- function(obj, x.tr, new.x) {
  p <- ncol(obj$Shat)
  temp <- sapply(1:p, function(i) {
    approx(x.tr[,i], obj$Shat[,i], new.x[,i], rule = 2)$y
  })
  obj$alpha0+apply(temp,1,sum)
}

args <-  commandArgs(T)
seed <- as.numeric(args[[1]])

# source("genData.R")
# source("waveMesh.R")
# n <- 2^8
# dist.x <- "unif"
# f.name <- "f.poly"

dataAdditive(seed = seed)

q(save = "no")
