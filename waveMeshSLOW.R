genR <- function(x, nlevel = 5) {
  #require(Matrix)
  # This function creates the interpolation matrix.
  
  # Get total number of intervals.
  K <- 2^nlevel
  # Create vector keeping track of intervals.
  k.vec <- as.numeric(cut(x, breaks = (0:K)/K, labels = 1:K)) - 1
  
  n <- length(x)
  # There HAS to be a better way to do this!
  R <- matrix(0, ncol = K, nrow = n)
  for(i in 1:n) {
    if(k.vec[i] == 0) {
      R[i, k.vec[i] + 1] <- 1
    } else {
      R[i, k.vec[i]] <- (1 + k.vec[i]- K*x[i])
      R[i, k.vec[i]+1] <- K*x[i] - k.vec[i]
    }
  }
  return(R)
}

next.iter <- function(obj1, obj2, w.theta, lambda, step,
                      filter.number = 8, family = "DaubExPhase",
                      bc = "interval", min.lev = 1) {
  # Obj1: (I - step*t(R)%*%R)
  # Obj2: step * t(R) %*% y
  # w.theta: The object W %*% theta^(l-1), i.e. the fitted function 
  #           from the previous iteration.
  
  z <- obj1 %*% w.theta + obj2
  wty <- wd(data = z, filter.number = filter.number, 
            family = family, bc = bc)
  th.obj <- threshold(wty, levels = min.lev:(nlevelsWT(wty)-1), 
                      type = "soft",
                      policy = "manual", 
                      value = step*lambda)
  if(bc == "interval") {
    return(wr.int(th.obj))
  } else {
    return(wr(th.obj))
  }
  
}

waveMesh <- function(y, x, lams = NULL, nlam = 50,
                     max.lam = 1, lam.min.ratio = 1e-3, 
                     nlevel = 7, 
                     max.iter = 100, tol = 1e-4,
                     theta0 = NULL, filter.number = 8, j0 = 1,
                     family = "DaubExPhase", bc = "interval") {
  require(wavethresh)
  nlevel <- min(nlevel, floor(log2(length(y))))
  # First we calculate the R and R^TR matrix and get the maximum eigen-value.
  # The step size 1/Max.eigen can be used.
  R <- genR(x, nlevel = nlevel)
  RtR <- t(R) %*% R
  step <- 1/max(eigen(RtR, symmetric = TRUE, only.values = TRUE)$values)
  
  # Some fixed vectors and matrices used throughout the algorithm
  #
  # One quantity of interest is Z^tilde = (I - t*t(R)*R)W\theta + t*t(R)*y
  #
  # First we have the first term (I - t*t(R)*R)
  obj1 <- -step*RtR
  diag(obj1) <- diag(obj1) + 1
  # Second we have t*t(R) %*% y
  obj2 <- step*t(R) %*% y
  
  # Initialize the theta vector if not intialized via warm starts.
  if(is.null(theta0)) {
    theta0 <- t(R) %*% y
  }
  if(is.null(lams)) {
    lams <- 10^seq(log10(max.lam), log10(max.lam*lam.min.ratio), 
                   length = nlam)
  } else {
    nlam <- length(lams)
  }
  
  K <- 2^nlevel
  w.theta.hat <- matrix(NA, ncol = nlam, nrow = K)
  old.val <- theta0
  for(i in 1:nlam){
    conv <- FALSE
    count <- 1
    while(!conv & count <= max.iter) {
      new.val <- next.iter(obj1, obj2, old.val, 
                           lambda = lams[i],
                           step = step, filter.number = filter.number,
                           bc = bc, family = family, min.lev = j0)
      delta <- mean((new.val - old.val)^2)/mean(old.val^2)
      cat("Lambda NUmber:", i, "; Delta:", delta,"\n")
      if(delta <= tol) {
        conv <- TRUE
        old.val <- new.val
        w.theta.hat[, i] <- new.val
      } else {
        old.val <- new.val
      }
    }
  }
  
  return(R %*% w.theta.hat)
}
