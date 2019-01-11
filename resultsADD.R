process.res <- function(n = 128, dist = "unif") {
  dir <- "simAdditive"
  
  files <- paste0(dir, "/n", n,"_distx", 
                  dist, "_seed", 1:100, ".RData")
  nf <- length(files)
  ans <- c()
  for(i in 1:nf){
    #print(i)
    load(files[i])
    ans <- rbind(ans, res)
  }
  
  mu <- colMeans(ans)
  sd <- apply(ans, 2, function(x){
    sd(x)/sqrt(nf)
  })
  return(list("mu" = round(mu[-1],2), "sd" = round(sd[-1],2)))
}

process.res(n = 64)
process.res(n = 100)
process.res(n = 128)
process.res(n = 256)
process.res(n = 500)
process.res(n = 512)
