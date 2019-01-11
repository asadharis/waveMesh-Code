
# 8/1/2018: File modified a bit to allow for sample sizes not a powr of 2
#           Done in response to reveiwer comments.

gen.res <- function(dist.x = "unif") {
  dat <- data.frame()
  # By default this file generates Table 1 of main manuscript.
  # for other tables please comment/uncomment code accordingly.
  
  for(n in 2^(6:9)) { 
    #for(n in c(75,100,300, 500)){ 
    for(f in c("f.poly", "f.sin", "f.ppoly","f.heavi","f.bumps","f.doppler")) {
      
      file <- paste0("sim2/n", n, "_distx", dist.x, "_f", f, ".RData")
      #file <- paste0("sim3/n", n, "_distx", dist.x, "_f", f, ".RData")
      
      load(file)
      names(res) <- c("4 Grid", "5 Grid", "6 Grid", "Full Grid",
                      "Interpolation", "Isometric", "Adaptive Lifting")
      res <- cbind(res, "n" = rep(paste("n =",n), nrow(res)),
                   "Function" = rep(f, nrow(res)))
      
      dat <- rbind(dat, res)
      
    }
  }
  
  return(dat)
}

#md <- gen.res("norm")
md <- gen.res("unif")
temp <- apply(md[,1:7], 1, function(x){
  x/x[4]
})

md2 <- md
md2[, 1:7] <- t(temp)

fmt <- function(x, digits, ...) {
  s <- format(x, digits=digits, ...)
  is_stderr <- (1:length(s)) > length(s) %/% 2
  s[is_stderr] <- sprintf("(%s)", s[is_stderr])
  s[!is_stderr] <- latexNumeric(s[!is_stderr])
  s
}
se <- function(x) {
  round(sd(x)/sqrt(length(x)) * 100,2)
}
mu <- function(x){
  paste0(format(round(mean(x),2), nsmall=2), " (", format(round(se(x), 2), nsmall=2), ")")
}

library(tables)
tab <- tabular(Function*n ~ (`5 Grid` + `6 Grid` + 
                               Interpolation + Isometric + 
                               `Adaptive Lifting`) * (mu)
               , data = md2)
latex(tab)
