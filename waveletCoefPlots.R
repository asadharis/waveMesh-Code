source("genData.R")
library(wavethresh)

plot.wavecoef <- function(f = f.poly, ...) {
  j <- 10
  n <- 2^j
  x <- (1:n)/n
  wdobj <- wd(f(x), filter.number = 8, family = "DaubExPhase",
               bc = "symmetric")
  plot(wdobj, scaling = "global", ...)
  
}

# plot.wavecoef(f.poly, main = "Polynomial Function")
# plot.wavecoef(f.exp, main = "Exponential Function")
# plot.wavecoef(f.sin, main = "Sine Function")
# 
# plot.wavecoef(f.ppoly, main = "Piecewise Polynomial")
# plot.wavecoef(f.doppler, main = "Doppler")
# plot.wavecoef(f.heavi, main = "Heaviside")
# plot.wavecoef(f.bumps, main = "Bumps")
# plot.wavecoef(f.blocks, main = "Blocks")
