
next.iter2 <- function(obj1, obj2, 
                       w.theta.inter, w.theta, 
                       count, lambda, step,
                       filter.number = 8, family = "DaubExPhase",
                       bc = "interval", min.lev = 1) {
  # Obj1: (I - step*t(R)%*%R)
  # Obj2: step * t(R) %*% y
  # w.theta: The object W %*% theta^(l-1), i.e. the fitted function 
  #           from the previous iteration.
  K <- length(w.theta)
  
  z <- tridiag_y(obj1, w.theta.inter, K) + obj2
  wty <- wd(data = z, filter.number = filter.number, 
            family = family, bc = bc)
  th.obj <- threshold(wty, levels = min.lev:(nlevelsWT(wty)-1), 
                      type = "soft",
                      policy = "manual", 
                      value = step*lambda)
  
  if(bc == "interval") {
    new.w.theta <- wr.int(th.obj)
  } else {
    new.w.theta <- wr(th.obj)
  }
  
  new.w.theta.inter <- new.w.theta + (new.w.theta - w.theta)*((count - 1)/(count + 2))
  list("new.inter" = new.w.theta.inter, "new.val" = new.w.theta)
}

waveMesh2 <- function(y, x, lams = NULL, nlam = 50,
                      max.lam = NULL, lam.min.ratio = 1e-3, 
                      nlevel = 7, 
                      max.iter = 100, tol = 1e-4,
                      theta0 = NULL, filter.number = 8, j0 = 1,
                      family = "DaubExPhase", bc = "interval") {
  
  # lams = NULL; nlam = 50;
  # max.lam = NULL; lam.min.ratio = 1e-3;
  # nlevel = 8;
  # max.iter = 100; tol = 1e-4;
  # theta0 = NULL; filter.number = 8; j0 = 1;
  # family = "DaubExPhase"; bc = "interval"
  
  require(wavethresh)
  K <- 2^nlevel
  n <- length(y)
  nlevel <- min(nlevel, floor(log2(length(y))))
  
  # First we calculate the R and R^TR matrix and get the maximum eigen-value.
  # The step size 1/Max.eigen can be used.
  # First evaluate vector of cluster membership
  k.vec <- gen_k_vec(x, K)
  k.ind <- gen_k_ind(x, K)
  spR <- gen_r(k.vec, x, K)
  spRtR <- rtr(spR, k.ind, K)
  step <- as.numeric(1/max_eigen(spRtR, K))
  
  # One quantity of interest is Z^tilde = (I - t*t(R)*R)W\theta + t*t(R)*y
  #
  # First we have the first term (I - t*t(R)*R)
  # This will be saved as a sparse matrix, i.e. will be
  # a matrix with two columns representing the diagonal and off-diagnal entries
  obj1 <- -step*spRtR
  obj1[,1] <- obj1[,1] + 1
  
  # Second we have t*t(R) %*% y
  obj2 <- step * rty(spR, y, k.ind, K)
  
  # Initialize the W%*%theta vector if not specified.
  if(is.null(theta0)) {
    theta0 <- obj2*0
  }
  
  if(is.null(lams)) {
    if(is.null(max.lam)) {
      temp <- wd(obj2, filter.number = filter.number, 
                 family = family, bc = bc)
      coefs <- sapply(j0:(nlevel-1), function(k){accessD(temp, k)})
      max.lam <- max(abs(do.call("c", coefs)))
      lams <- 10^seq(log10(max.lam), log10(max.lam*lam.min.ratio), 
                     length = nlam)
    } else {
      lams <- 10^seq(log10(max.lam), log10(max.lam*lam.min.ratio), 
                     length = nlam)
    }
  } else {
    nlam <- length(lams)
  }
  
  w.theta.hat <- matrix(NA, ncol = nlam, nrow = K)
  old.val <- theta0
  for(i in 1:nlam){
    conv <- FALSE
    count <- 1
    old.iter <- old.val
    while(!conv & count <= max.iter) {
      new.val <- next.iter2(obj1, obj2, old.iter, old.val,
                            count, lambda = lams[i],
                            step = step, filter.number = filter.number,
                            bc = bc, family = family, min.lev = j0)
      
      count <- count + 1
      delta <- mean((new.val$new.val - old.val)^2)/mean(old.val^2 + 1e-15)
      cat("Lambda NUmber:", i, "; Delta:", delta,"\n")
      if(delta <= tol | count == max.iter) {
        conv <- TRUE
        old.val <- new.val$new.val
        w.theta.hat[, i] <- new.val$new.val
      } else {
        old.val <- new.val$new.val
        old.iter <- new.val$new.inter
      }
    }
  }
  
  return(rtheta(spR, w.theta.hat, k.ind, K, n)) 
}

