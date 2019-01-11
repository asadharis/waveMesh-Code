source("addWaveMesh.R")
source("genData.R")
source("AMlet.r")
source("wavethresh.r")
# This file consists of the main simulation study for the aditive model. 
# We use the simulation setting of Sardy et. al (2004)
# We implement our method we select the lambda value using CV. 

simAdditive <- function(n = 2^12, seed = 1, dist.x = "unif") {
  #n = 500; seed = 1; dist.x = "unif"
  
  require(wavethresh)
  # Set seed to keep fixed design. 
  set.seed(1)
  if(dist.x == "unif") {
    x <- matrix(runif(n*4), ncol = 4)
  } else if(dist == "beta") {
    x <- matrix(rbeta(n*4, shape1 = 2, shape2 = 5), ncol = 4)
  }
  
  f1 <- f.ppoly(x[,1]); f2 <- f.heavi(x[,2]) 
  f3 <- f.poly(x[,3]); f4 <- f.sin(x[,4])
  # f1 <- f1-mean(f1); f2 <- f2-mean(f2); f3 <- f3 - mean(f3)
  # f4 <- f4-mean(f4)
  
  snr <- 10
  var.sd <- var(f1+f2+f4+f3)/snr
  set.seed(seed)
  y <- f1 +f2 + f3 + f4 +  rnorm(n, sd = sqrt(var.sd))
  # First we implement out proposal via cross validation.
  lam.seq <- 10^(seq(log10(30), log10(0.1), length = 50))
  
  cv.wm <- cv.addWaveMesh(y, x, K = 5, max.iter.outer = 1000, 
                          tol.outer = 1e-6, max.iter = 1000, tol = 1e-6,
                          lam.seq = lam.seq, gamma = 0,
                          family = "DaubLeAsymm" , filter.number=8, bc = "periodic",
                          j0 = 0, nlevel = ceiling(log2(n)))
  
  fit.wm <- fit.additive(y, x, max.iter.outer= 1000, tol.outer = 1e-6, 
                         lam.seq = c(cv.wm$lam.min, cv.wm$lam.1se),
                         gamma = 0, family = "DaubLeAsymm" , filter.number=8, bc = "periodic",
                         j0 = 0, nlevel =  ceiling(log2(n)), max.iter = 1000, tol = 1e-6)
  fhat.wm <- predict.add_mod(fit.wm, new.data = x, type = "response")
  
  mse <- colMeans( (fhat.wm - (f1+f2+f3+f4))^2 )
  
  
  ############## Updated to allow for non-powers of 2 ################
  nJ <- ceiling(log2(length(y)))
  nz <- 2^nJ - length(y)
  y.new <- c(rep(NA, floor(nz/2)), y, rep(NA, ceiling(nz/2)))
  ind <- which(!is.na(y.new))
  y.new[is.na(y.new)] <- 0
  x.new <- apply(x, 2, function(col){
    rng <- range(col)
    c(rep(rng[1], floor(nz/2)), col, rep(rng[2], ceiling(nz/2)))
    
  })
    
  ############## Updated to allow for non-powers of 2 ################

    # Followed by Amlet
  fit.amlet1 <- AMlet(y.new, x.new, thresh.value="universal", 
                      sigma=NA, p=1, low.level=0, by.level=F, 
                      filter.number=8, family="DaubLeAsymm",
                      bc="periodic", conv.thresh=1.e-10, max.iter=ncol(x)*1000)
  
  mse.amlet1 <- mean((fit.amlet1$muhat[ind] - (f1+f2+f3+f4))^2)
  
  res <- c("mseMin" = mse[1], "mse1se" = mse[2], "mseAM" = mse.amlet1)
  
  dir <- paste0("simAdditive")
  filename <- paste0(dir, "/n", n, "_distx", dist.x, "_seed", seed, ".RData")
  
  if(dir.exists(dir)) {
    save(res, file = filename)
  } else {
    dir.create(dir)
    save(res, file = filename)
  }
}


args <-  commandArgs(T)
n <- as.numeric(args[[1]])
dist.x <- as.character(args[[2]])
seed <- as.numeric(args[[3]])

# source("genData.R")
# source("waveMesh.R")
# n <- 2^8
# dist.x <- "unif"
# f.name <- "f.poly"

simAdditive(n = n, seed = seed, dist.x = dist.x)

q(save = "no")
