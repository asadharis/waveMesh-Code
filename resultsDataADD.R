process.res <- function() {
  dir <- "simAdditiveData"
  files <- paste0(dir, "/",list.files(path = dir))
  n <- length(files)
  ans <- c()
  for(i in 1:n){
    load(files[i])
    ans <- rbind(ans, res)
  }
  
  colMeans(ans)
  
}
