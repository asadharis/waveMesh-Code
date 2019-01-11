# In this file we aggregate the results and draw appropriate box plots.

gen.res <- function(dist.x = "unif", f.name = "f.poly") {
  #dist.x <- "unif"; f.name = "f.sin"; 
  dir <- paste0("sim1")
  ns <- 2^(5:12)
  files <- paste0("sim1/n",ns, "_distx", dist.x, "_f", f.name, ".RData")
  mydat <- data.frame()
  meth <- c("4 Grid", "5 Grid", "6 Grid", "Full Grid")
  for(i in 1:length(ns)) {
    load(files[i])
    temp.dat <- data.frame("MSE" = c(ans[,1:4]), "Time" = c(ans[,5:8]),
                           "Method" = rep(meth, each = 100), 
                           "n" = ns[i])
    mydat <- rbind(mydat, temp.dat)
  }
  mydat$n <- factor(mydat$n)
  
  # Now we also plot the example plots.
  load(paste0("sim1/n",1024, "_distx", dist.x, "_f", f.name, ".RData"))
  
  mydat2 <- data.frame("x" = rep(x, 3),
                       "f" = c(true.f, fits[, "f6"], 
                               fits[, "ff"]),
                       "Method" = c(rep("True Function", 1024),
                                    rep("6 Grid", 1024), 
                                    rep("Full Grid", 1024)))
  
  require(ggplot2)
  p1 <- ggplot(data = mydat, aes(x=n, y=MSE)) + 
    geom_boxplot(aes(fill=Method), outlier.shape = NA) + 
    xlab("Sample Size") + ylab("Mean Square Error") +
    theme_grey(base_size = 18) + scale_y_log10()+
    ggtitle("Prediction Performance")
  p2 <- ggplot(data = mydat, aes(x=n, y=Time)) + 
    geom_boxplot(aes(fill=Method), outlier.shape = NA) + 
    xlab("Sample Size") + ylab("Evaluation Time (s)") +
    theme_grey(base_size = 18) + ggtitle("Timing Results")
  p3 <- ggplot(data = mydat2, aes(x=x, y=f, color = Method))+
    geom_line(size = 1) + 
    xlab("x") + ylab("f(x)") + 
    guides(color = guide_legend(title=NULL))+
    theme_grey(base_size = 18) + 
    ggtitle("Estimated Functions")
  list(p1,p2, p3)
}


gen.res2 <- function(dist.x = "unif", f.name = "f.poly") {
  #dist.x <- "unif"; f.name = "f.sin"; 
  dir <- paste0("simAdap")
  ns <- 2^(5:12)
  files <- paste0("simAdap/n",ns, "_distx", dist.x, "_f", f.name, ".RData")
  mydat <- data.frame()
  meth <- c("Regular", "Adaptive")
  for(i in 1:length(ns)) {
    load(files[i])
    temp.dat <- data.frame("MSE" = c(ans[,1:2]),
                           "Method" = rep(meth, each = 100), 
                           "n" = ns[i])
    mydat <- rbind(mydat, temp.dat)
  }
  mydat$n <- factor(mydat$n)
  
  # Now we also plot the example plots.
  load(paste0("simAdap/n",1024, "_distx", dist.x, "_f", f.name, ".RData"))
  
  mydat2 <- data.frame("x" = rep(x, 3),
                       "f" = c(true.f, fits[, "ff"], 
                               fits[, "fa"]),
                       "Method" = c(rep("True Function", 1024),
                                    rep("Regular", 1024), 
                                    rep("Adaptive", 1024)))
  mydat2$Method <- factor(mydat2$Method, levels(mydat2$Method)[c(3,2,1)])
  require(ggplot2)
  p1 <- ggplot(data = mydat, aes(x=n, y=MSE)) + 
    geom_boxplot(aes(fill=Method), outlier.shape = NA) + 
    xlab("Sample Size") + ylab("Mean Square Error") +
    theme_grey(base_size = 18) + scale_y_log10()+
    ggtitle("Prediction Performance")
  
  p2 <- ggplot(data = mydat2, aes(x=x, y=f, color = Method))+
    geom_line(size = 1) + 
    xlab("x") + ylab("f(x)") + 
    guides(color = guide_legend(title=NULL))+
    theme_grey(base_size = 18) + 
    ggtitle("Estimated Functions")
  list(p1,p2)
}


results.plot <- function(f.name = "f.poly", 
                         dist.x = "unif",
                         main = "Polynomial") {
  library(ggplot2)
  source("waveletCoefPlots.R")
  obj <- gen.res(dist.x, f.name)
  
  pdf(paste0("plots/", main, ".pdf"))
  print(obj[[1]])
  print(obj[[2]])
  print(obj[[3]])
  p1 <- plot.wavecoef(get(f.name), main = main)
  dev.off()
}

results.plot2 <- function(f.name = "f.poly", 
                         dist.x = "unif",
                         main = "Polynomial") {
  library(ggplot2)
  obj <- gen.res2(dist.x, f.name)
  
  pdf(paste0("plotsAdap/", main, ".pdf"))
  print(obj[[1]])
  print(obj[[2]])
  dev.off()
}


results.plot2("f.poly", "unif", "Polynomial")
results.plot2("f.sin", "unif", "Sine")
#results.plot("f.exp", "unif", "Exponential")

results.plot2("f.bumps", "unif", "Bumps")
#results.plot("f.blocks", "unif", "Blocks")
results.plot2("f.doppler", "unif", "Doppler")
results.plot2("f.heavi", "unif", "Heaviside")
results.plot2("f.ppoly", "unif", "PiecewisePolynomial")

