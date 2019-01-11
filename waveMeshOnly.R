
simulation1 <- function(n = 512, nsim = 100, SNR = 5, 
                        dist.x = "unif", f.name = "f.poly",
                        max.lam = NULL,  ...) {
  # n = 2^15; nsim = 50; SNR = 3; f.name = "f.poly";
  # dist.x = "unif"; max.lam = 50; #min.lam =1
  
  f <- get(f.name)
  dat <- genData(nsim, n, SNR, f, dist.x = dist.x)
  
  fhat4 <- matrix(0, nrow = n, ncol = nsim)
  fhat5 <- fhat6 <- fhatFull <- fhat4
  
  mse4 <- mse5 <- mse6 <- mseFull <- numeric(nsim)
  time4 <- time5 <- time6 <- timeFull <- numeric(nsim)
  for(i in 1:nsim) {
    print(i)
    y <- dat$y[,i]
    x <- dat$x
    
    # Fit models.
    time4[i] <- system.time(mod4 <- waveMesh(y, x, max.lam = NULL,
                                             nlevel = 4))[3]
    time5[i] <- system.time(mod5 <- waveMesh(y, x, max.lam = NULL,
                                             nlevel = 5))[3]
    time6[i] <- system.time(mod6 <- waveMesh(y, x, max.lam = NULL,
                                             nlevel = 6))[3]
    timeFull[i] <- system.time(modFull <- waveMesh(y, x, max.lam = NULL,
                                             nlevel = floor(log2(n))))[3]
    
    # Calculate MSEs.
    mse4[i] <- min(colMeans((dat$f - mod4$ans)^2))
    mse5[i] <- min(colMeans((dat$f - mod5$ans)^2))
    mse6[i] <- min(colMeans((dat$f - mod6$ans)^2))
    mseFull[i] <- min(colMeans((dat$f - modFull$ans)^2))
    
    fhat4[,i] <- mod4$ans[,which.min(colMeans((dat$f - mod4$ans)^2))]
    fhat5[,i] <- mod5$ans[,which.min(colMeans((dat$f - mod5$ans)^2))]
    fhat6[,i] <- mod6$ans[,which.min(colMeans((dat$f - mod6$ans)^2))]
    fhatFull[,i] <- modFull$ans[,which.min(colMeans((dat$f - modFull$ans)^2))]
  }

  fits <- cbind("f4" = fhat4[, which.min(mse4)], "f5" = fhat5[, which.min(mse5)],
                "f6" = fhat6[, which.min(mse6)], "ff" = fhatFull[, which.min(mseFull)])
  ans <- cbind("m4" = (mse4), "m5" = (mse5),
         "m6" = (mse6), "mF" = (mseFull), 
         "m4T" = time4, "m5T" = time5, "m6T" = time6,
         "mFT" = timeFull)
  dir <- paste0("sim1")
  filename <- paste0("sim1/n", n, "_distx", dist.x, "_f", f.name, ".RData")
  true.f <- dat$f
  x <- dat$x
  if(dir.exists(dir)) {
    save(ans, fits, true.f, x, file = filename)
  } else {
    dir.create(dir)
    save(ans,fits, true.f, x, file = filename)
  }
  
}


simulation2 <- function(n = 512, nsim = 100, SNR = 5, 
                        dist.x = "unif", f.name = "f.poly",
                        max.lam = NULL,  ...) {
  # This simulation compares wavemesh with adaptive wavemesh.
  
  # n = 2^9; nsim = 50; SNR = 5; f.name = "f.poly";
  # dist.x = "unif";

  f <- get(f.name)
  dat <- genData(nsim, n, SNR, f, dist.x = dist.x)
  
  fhat <- fhatAdap <- matrix(0, nrow = n, ncol = nsim)

  mse <- mseAdap <-  numeric(nsim)
  for(i in 1:nsim) {
    print(i)
    y <- dat$y[,i]
    x <- dat$x
    
    # Fit models.
    modFull <- waveMesh(y, x, max.lam = NULL, nlevel = floor(log2(n)))
    modFullAdap <- waveMesh(y, x, max.lam = NULL, nlevel = floor(log2(n)),
                            adaptive = TRUE)
    
    # Calculate MSEs.
    mse[i] <- min(colMeans((dat$f - modFull$ans)^2))
    mseAdap[i] <- min(colMeans((dat$f - modFullAdap$ans)^2))
    
    fhat[,i] <- modFull$ans[,which.min(colMeans((dat$f - modFull$ans)^2))]
    fhatAdap[, i] <- modFullAdap$ans[,which.min(colMeans((dat$f - modFullAdap$ans)^2))]
  }
  
  fits <- cbind("ff" = fhat[, which.min(mse)], "fa" = fhatAdap[, which.min(mseAdap)])
  ans <- cbind("mF" = (mse), 
               "mA" = mseAdap)
  dir <- paste0("simAdap")
  filename <- paste0("simAdap/n", n, "_distx", dist.x, "_f", f.name, ".RData")
  true.f <- dat$f
  x <- dat$x
  if(dir.exists(dir)) {
    save(ans, fits, true.f, x, file = filename)
  } else {
    dir.create(dir)
    save(ans,fits, true.f, x, file = filename)
  }
  
}


args <-  commandArgs(T)
n <- as.numeric(args[[1]])
dist.x <- as.character(args[[2]])
f.name <- as.character(args[[3]])

# n <- 2^8
# nsim <- 2
# dist.x <- "unif"
# f.name = "f.poly"
source("genData.R")
source("waveMesh.R")

# simulation1(n = n, nsim = 100, SNR = 5, dist.x = dist.x, 
#             f.name = f.name, max.lam = NULL)
simulation2(n = n, nsim = 100, SNR = 5, dist.x = dist.x, 
            f.name = f.name, max.lam = NULL)
q(save = "no")