
# This file contains function(s) for generating the data

genData <- function(nsim  = 1, n = 5e+2, SNR = 5, f = f.blocks,
                    dist.x = "unif", ...) {
  set.seed(1)
  if(dist.x == "unif") {
    x <- runif(n)
  } else if(dist.x == "norm") {
    x <- rnorm(n, ...)
    x <- (x-min(x))/diff(range(x))
  } else if(dist.x == "exp") {
    x <- rexp(n, ...)
    x <- (x-min(x))/diff(range(x))
  } else if(dist.x == "grid") {
    x <- (1:n)/n
  } else if(dist.x == "mix") {
    components <- sample(1:2,prob=c(0.5,0.5),size=n,replace=TRUE)
    mus <- c(-2,2)
    x <- rnorm(n = n, mean = mus)
    x <- (x-min(x))/diff(range(x))
  }
  
  x <- sort(x)
  myf <- f(x)
  var.noise <- var(myf)/SNR
  
  #set.seed(seed)
  y <- matrix(myf, ncol = nsim, nrow = n) + rnorm(n*nsim, sd = sqrt(var.noise))
  return(list("x" = x, "y" = y, 
              "f" = myf, "sigma" = sqrt(var.noise)))
}

f.blocks <- function(x) {
  n <- length(x)
  t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 
         0.78, 0.81)
  h1 <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  blocks <- rep(0, n)
  for (i in seq(1, length(h1))) {
    blocks <- blocks + (h1[i] * (1 + sign(x - t[i])))/2
  }
  blocks/sqrt(var(blocks)) * 7
}
 
f.bumps <- function(x) {
  n <- length(x)
  h2 <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
  t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 
         0.78, 0.81)
  w <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 
         0.005, 0.008, 0.005)
  bumps <- rep(0, n)
  for (i in seq(1, length(h2))) {
    bumps <- bumps + h2[i] * pmax(0, (1 - abs((x - t[i])/w[i])))^4
  }
  bumps/sqrt(var(bumps)) * 7
}

f.heavi <- function(x) {
  heavi <- 4 * sin(4 * pi * x) - 
    sign(x - 0.3) - sign(0.72 - x)
  heavi/sqrt(var(heavi)) * 7
}

f.doppler <- function(x) {
  eps <- 0.05
  doppler <- sqrt(x * (1 - x)) * 
    sin((2 * pi * (1 - eps))/(x + eps))
  doppler/sqrt(var(doppler)) * 7
}

f.ppoly <- function(x)  {
  # x <- seq(0, 1, length = 513)
  # x <- x[1:512]
  y <- rep(0, length(x))
  xsv <- (x <= 0.5)
  y[xsv] <- -16 * x[xsv]^3 + 12 * x[xsv]^2
  xsv <- (x > 0.5) & (x <= 0.75)
  y[xsv] <- (x[xsv] * (16 * x[xsv]^2 - 40 * x[xsv] + 28))/3 - 
    1.5
  xsv <- x > 0.75
  y[xsv] <- (x[xsv] * (16 * x[xsv]^2 - 32 * x[xsv] + 16))/3
  y
}


f.poly <- function(x) {
  poly <- (x - 0.2) * (x - 0.7)
  poly/sqrt(var(poly)) * 7
}

f.sin <- function(x) {
  sf <- 3*sin(6*x)
  sf/sqrt(var(sf)) * 7
}

f.exp <- function(x) {
  ef <- exp(-4*x) - 0.5
  ef/sqrt(var(ef)) * 7
}


# x <- seq(0,1,length = 10000)
# par(mfrow = c(2,3))
# # Easy functions
# plot(x, f.poly(x), type = "l", lwd = 3, col = "green",
#      xlab = "x", ylab = "f(x)")
# 
# # Medium functions
# plot(x, f.ppoly(x), type = "l", lwd = 3, col = "orange",
#      xlab = "x", ylab = "f(x)")
# 
# # Hard functions
# plot(x, f.doppler(x), type = "l", lwd = 3, col = "red",
#      xlab = "x", ylab = "f(x)")
# 
# 
# 
# plot(x, f.sin(x), type = "l", lwd = 3, col = "green",
#      xlab = "x", ylab = "f(x)")
# 
# plot(x, f.heavi(x), type = "l", lwd = 3, col = "orange",
#      xlab = "x", ylab = "f(x)")
# 
# 
# plot(x, f.bumps(x), type = "l", lwd = 3, col = "red",
#      xlab = "x", ylab = "f(x)")
# plot(x, f.blocks(x), type = "l", lwd = 3, col = "red",
#      xlab = "x", ylab = "f(x)")

