# 8/1/2018: File created as part of additional simulations for 
#           response to reviewer comments. 
#           The methods here modify the Sardy et. al. method to 
#           allow for non powers of 2 data. 
#           
#           This is achieved by padding the response vector by zeros.


# In this file we run a simulation comparing our proposal to existing methods.
# This simulations considers two methods for non-equividistant data with wavelets
# we do this for a sequence of tuning parameter values and compare methods MSE as 
# a function of the degrees of freedom.
# For the second method, in the adlift library, we cannot do this for 
# a sequence of tuning parameters and hence we will represent the results of this method as
# a single point on the graph. 
library(wavethresh)
irr.reg <- function(x, y, lams, filter.number = 8, j0 = 1,
                    family = "DaubExPhase", bc = "symmetric", ...) {
  # This function performs wavlet regression for unequal spaced data 
  # using an interpolation method.
  nlam <- length(lams)
  gd <- makegrid(t = x, y = y)
  fhat <- matrix(NA, ncol = length(lams), nrow = length(y))
  for(i in 1:nlam){
    wty <- irregwd(gd, family= family, filter.number = filter.number,
                   bc = bc)
    th.obj <- threshold(wty, levels = j0:(nlevelsWT(wty)-1), 
                        type = "soft",
                        policy = "manual", 
                        value = lams[i])
    if(bc == "interval") {
      new.val <- wr.int(th.obj)
    } else {
      new.val <- wr(th.obj)
    }
    fhat[, i] <- approx(gd$gridt, new.val, xout = x, rule = 2)$y
  }
  return(fhat)
}

pad.na <- function(x) {
  n <- length(x)
  J <- ceiling(log2(n))
  nz <- 2^J - n
  c(rep(NA, floor(nz/2)), x, rep(NA, ceiling(nz/2)) )
}

irr.reg2 <- function(x, y, lams, filter.number = 8, j0 = 1,
                     family = "DaubExPhase", bc = "symmetric", ...) {
  # This function performs wavlet regression for unequal spaced data 
  # using the isometric wavelet method of Sardy et al.
  # According to paper, this amounts to doing a DWT treated the sorted 
  # data as if it were equally spaced on the grid, i.e. just ignore 
  # irregularities.
  
  # 8/1/2018: One major change here is that we allow for sample size not
  #           a power of 2. 
  
  nlam <- length(lams)
  ny <- length(y)
  fhat <- matrix(NA, ncol = nlam, nrow = ny)
  
  y.pad <- pad.na(y)  
  ind.y <- which(!is.na(y.pad))
  
  y.pad[is.na(y.pad)] <- 0
  for(i in 1:nlam){
    wty <- wd(y.pad, family= family, 
              filter.number = filter.number,
              bc = bc)
    th.obj <- threshold(wty, levels = j0:(nlevelsWT(wty)-1), 
                        type = "soft",
                        policy = "manual", 
                        value = lams[i])
    if(bc == "interval") {
      new.val <- wr.int(th.obj)
    } else {
      new.val <- wr(th.obj)
    }
    fhat[, i] <- new.val[ind.y]
  }
  return(fhat)
}


simulation3 <- function(n = 2^9, nsim = 50, SNR = 5, 
                        dist.x = "unif", f.name = "f.poly", ...) {
  require(wavethresh)
  require(adlift)
  # n = 128; nsim = 100; SNR = 5; f.name = "f.bumps";
  # dist.x = "exp"; 
  f <- get(f.name)
  dat <- genData(nsim, n, SNR, f, dist.x = dist.x)
  m4 <- m5 <- m6 <- mFull <- mIR <- numeric(nsim)
  mAdlift <- mIso <- numeric(nsim)
  for(i in 1:nsim) {
    print(i)
    x <- dat$x
    y <- dat$y[,i]
    f <- dat$f
    
    # We begin with out proposal
    mod4 <- waveMesh(y, x, max.lam = NULL,
                     nlevel = 4,max.iter = 500)
    mod5 <- waveMesh(y, x, max.lam = NULL,
                     nlevel = 5,max.iter = 500)
    mod6 <- waveMesh(y, x, max.lam = NULL,
                     nlevel = 6, max.iter = 500)
    modFull <- waveMesh(y, x, max.lam = NULL,
                        nlevel = floor(log2(n)),
                        max.iter = 500)
    
    m4[i] <- min(colMeans((f - mod4$ans)^2))
    m5[i] <- min(colMeans((f - mod5$ans)^2))
    m6[i] <- min(colMeans((f - mod6$ans)^2))
    mFull[i] <- min(colMeans((f - modFull$ans)^2))
    
    # Now we present results for the wavethres based methods.
    # I.e. interpolation and ignoring scaling.
    mod.irr <- irr.reg(x, y, modFull$lams)
    mIR[i] <- min(colMeans((f - mod.irr)^2))
    
    # Now results for ignoring the unqually spacedness
    mod.irr2 <- irr.reg2(x, y, modFull$lams)
    mIso[i] <- min(colMeans((f - mod.irr2)^2))
    
    # Finally we present results for adaptive lifting technique.
    mod.Adlift <- denoise(x, y, AdaptNeigh,1,TRUE,TRUE, 2, rule = "soft")
    mAdlift[i] <- mean((f - mod.Adlift)^2)
  }
  # y.mat <- dat$y
  # f <- dat$f
  # sigma <- dat$sigma
  res <- data.frame("4 Grid" = m4, "5 Grid" = m5, "6 Grid" = m6,
                    "Full Grid" = mFull, "Interpolaion" = mIR,
                    "Isometric" = mIso, "Adaptive Lifting" = mAdlift)
  dir <- paste0("sim3")
  filename <- paste0("sim3/n", n, "_distx", dist.x, "_f", f.name, ".RData")
  
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
f.name <- as.character(args[[3]])

source("genData.R")
source("waveMesh.R")
# n <- 2^8
# dist.x <- "unif"
# f.name <- "f.poly"
simulation3(n = n, nsim = 100, SNR = 5, dist.x = dist.x, 
            f.name = f.name)

q(save = "no")
