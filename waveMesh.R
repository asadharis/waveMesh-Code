# The main file for implementing the Meshy wavelet optimization.
# Our algorithm is based on the proximal gradient descent, by 
# using an appropriate majorizing surrogate. 

# A few points:
#   1. We mostly calculate a discrete wavelet transform. I.e. no big matrix multiplication
#   2. We do matrix multiplication with a sparse matrix R and (R^TR), the latter is tridiagonal
#   3. We can calculate the Lipchitz constant exactly, so no need to do a line search.

Rcpp::sourceCpp('helpers.cpp')
next.iter <- function(obj1, obj2, w.theta, lambda, step,
                      filter.number = 8, family = "DaubExPhase",
                      bc = "interval", min.lev = 1) {
  # Obj1: (I - step*t(R)%*%R)
  # Obj2: step * t(R) %*% y
  # w.theta: The object W %*% theta^(l-1), i.e. the fitted function 
  #           from the previous iteration.
  K <- length(w.theta)
  
  z <- tridiag_y(obj1, w.theta, K) + obj2
  wty <- wd(data = z, filter.number = filter.number, 
            family = family, bc = bc)
  th.obj <- threshold(wty, levels = min.lev:(nlevelsWT(wty)-1), 
                      type = "soft",
                      policy = "manual", 
                      value = step*lambda)
  mother.coef <- sapply(min.lev:(nlevelsWT(th.obj)-1), function(j){accessD(th.obj, j)})
  mother.coef <- do.call("c", mother.coef)
  #father.coef <- accessC(th.obj, min.lev)

  l1.norm <- sum(abs(mother.coef))
    
  if(bc == "interval") {
    new.val <- wr.int(th.obj)
  } else {
    new.val <- wr(th.obj)
  }
  list("new.val" = new.val, "l1" = l1.norm)
}

get.z <- function(w.theta, rty_vec, spRtR, step, lambda,
                  filter.number = 8, family = "DaubExPhase",
                  bc = "interval", min.lev = 1) {
  K <- length(w.theta)
  obj1 <- -step*spRtR
  obj1[,1] <- obj1[,1] + 1
  
  temp_y <- tridiag_y(obj1, w.theta, K) + step*rty_vec
  wty <- wd(data = temp_y, filter.number = filter.number, 
            family = family, bc = bc)
  th.obj <- threshold(wty, levels = min.lev:(nlevelsWT(wty)-1), 
                      type = "soft",
                      policy = "manual", 
                      value = step*lambda)
  
  mother.coef <- sapply(min.lev:(nlevelsWT(th.obj)-1), function(j){accessD(th.obj, j)})
  mother.coef <- do.call("c", mother.coef)
  #father.coef <- accessC(th.obj, min.lev)
  
  l1.norm <- sum(abs(mother.coef))
  
  if(bc == "interval") {
    new.val <- wr.int(th.obj)
  } else {
    new.val <- wr(th.obj)
  }
  return(list("new.val" = matrix(new.val), "l1" = l1.norm))
}

line.search <- function(w.theta, rty_vec, spRtR,  
                        spR, k.ind, step, lambda,
                        filter.number = 8, family = "DaubExPhase",
                        bc = "interval", min.lev = 1,
                        alpha = 0.5) {
  step0 <- step
  K <- length(w.theta)
  n <- nrow(spR)
  conv <- FALSE
  res.old <- y - rtheta(spR, w.theta, k.ind, K, n)
  loss.old <- 0.5*sum(res.old^2)
  while(!conv) {
    z <- get.z(w.theta, rty_vec, spRtR, step, lambda,
               filter.number = filter.number, 
               family = family,
               bc = bc, min.lev = min.lev)
    lhs <- 0.5*sum((y - rtheta(spR, z$new.val, k.ind, K, n))^2)
    rhs <- loss.old - t(res.old) %*% rtheta(spR, z$new.val - w.theta, k.ind, K, n) + 
      (1/(2*step))*sum((z$new.val - w.theta)^2)
    #print(c(lhs, rhs))
    if(lhs <= rhs) {
      conv <- TRUE
    } else {
      step <- alpha * step
    }
  }
  #print(c(step0, step))
  return(z)
}


waveMesh <- function(y, x, lams = NULL, nlam = 50,
                     max.lam = NULL, lam.min.ratio = 1e-3, 
                     nlevel = 7,
                     max.iter = 100, tol = 1e-4,
                     theta0 = NULL, filter.number = 8, j0 = 1,
                     family = "DaubExPhase", bc = "symmetric",
                     ls.alpha = 0.5) {
  
  # lams = NULL; nlam = 50;
  # max.lam = 3; lam.min.ratio = 1e-3;
  # nlevel = 8;
  # max.iter = 100; tol = 1e-4;
  # theta0 = NULL; filter.number = 8; j0 = 1;
  # family = "DaubExPhase"; bc = "periodic"
  
  require(wavethresh)
  K <- 2^nlevel
  n <- length(y)
  #nlevel <- min(nlevel, floor(log2(length(y))))

  # First we calculate the R and R^TR matrix and get the maximum eigen-value.
  # The step size 1/Max.eigen can be used.
  # First evaluate vector of cluster membership
  k.vec <- gen_k_vec(x, K)
  k.ind <- gen_k_ind(x, K)
  spR <- gen_r(k.vec, x, K)
  spRtR <- rtr(spR, k.ind, K)
  if(length(x)==K){
    if(all(sort(x) == (1:K)/K)) {
      step <- 1
    } else {
      step <- as.numeric(1/max_eigen(spRtR, K))
    }
  }else {
    step <- as.numeric(1/max_eigen(spRtR, K))
  }
  #step <- 10
  
  # One quantity of interest is Z^tilde = (I - t*t(R)*R)W\theta + t*t(R)*y
  #
  # First we have the first term (I - t*t(R)*R)
  # This will be saved as a sparse matrix, i.e. will be
  # a matrix with two columns representing the diagonal and off-diagnal entries
  obj1 <- -step*spRtR
  obj1[,1] <- obj1[,1] + 1
  
  # Second we have t*t(R) %*% y
  rty_vec <- rty(spR, y, k.ind, K)
  obj2 <- step * rty_vec
  
  # Initialize the W%*%theta vector if not specified.
  if(is.null(theta0)) {
    theta0 <- obj2*0
  }
  
  if(is.null(lams)) {
    if(is.null(max.lam)) {
      temp <- wd(obj2/step, filter.number = filter.number, 
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
  #print(lams[1])
  
  w.theta.hat <- matrix(NA, ncol = nlam, nrow = K)
  old.val <- theta0
  objective.val <- vector("list", nlam)
  conv.v <- c()
  for(i in 1:nlam){
    #cat("Lambda:", i, "\n")
    conv <- FALSE
    count <- 1
    while(!conv & count <= max.iter) {
      # new.val <- line.search(old.val, rty_vec, spRtR, spR, k.ind,
      #                        step = step, lambda = lams[i],
      #                        filter.number = filter.number,
      #                        family = family,
      #                        bc = bc, min.lev = j0,
      #                        alpha = ls.alpha)
      new.val <- next.iter(obj1, obj2, old.val, lambda = lams[i],
                            step = step, filter.number = filter.number,
                            bc = bc, family = family, min.lev = j0)
      

      obj <- 0.5 * sum((y - rtheta(spR, matrix(new.val$new.val), k.ind, K, n))^2) +
        lams[i] * new.val$l1
      objective.val[[i]] <- c(objective.val[[i]], obj)
      count <- count + 1
      delta <- mean((new.val$new.val - old.val)^2)/mean(old.val^2 + 1e-15)
      #cat("Lambda NUmber:", i, "; Delta:", delta,"\n")
      if(delta <= tol | count == max.iter) {
        conv <- TRUE
        old.val <- new.val$new.val
        w.theta.hat[, i] <- new.val$new.val
      } else {
        old.val <- new.val$new.val
      }
    }
    conv.v <- c(conv, conv.v)
  }
  
#  return(rtheta(spR, w.theta.hat, k.ind, K, n))
  return(list("ans" = rtheta(spR, w.theta.hat, k.ind, K, n),
              "obj" = objective.val,
              "lams" = lams,"conv" = conv.v))
}

predict.waveMesh <- function(obj, x.tr, new.x) {
  nlams <- length(obj$lams)
  apply(obj$ans, 2, function(vec) {
    approx(x.tr, vec, xout = new.x, rule = 2)$y
  })
}

solve.prox.waveMesh <- function(y, x, lambda1, lambda2, ...) {
  temp <- waveMesh(y, x, lams = lambda2, ...)
  #print(temp$conv)
  f_hat <- temp$ans
  max(1 - lambda1/sqrt(sum(f_hat^2)), 0)*f_hat
}


# A function to do cross validation for univariate waveMesh.
# By default it does 5-fold cross validation.
cv.waveMesh <- function(y, x, K = 5, lams = NULL, nlam = 50,
                     max.lam = NULL, lam.min.ratio = 1e-3, 
                     nlevel = 7,
                     max.iter = 100, tol = 1e-4,
                     theta0 = NULL, filter.number = 8, j0 = 1,
                     family = "DaubExPhase", bc = "symmetric",
                     ls.alpha = 0.5) {
  
  if(is.null(lams)) {
    lams <- 10^seq(log10(max.lam), log10(max.lam*lam.min.ratio), 
                     length = nlam)
  } else {
    nlam <- length(lams)
    max.lam <- max(lams)
    lam.min.ratio <- min(lams)/max(lams)
  }
  
  n <- length(y)
  # Sample the bins in case x is ordered. 
  set.seed(1)
  folds <- sample(cut(1:n, breaks = K, labels = F))
  cv.mse <- matrix(NA, ncol = nlam, nrow = K)
  
  for(j in 1:K) {
    cat("CV Fold:", j,"\n")
    x.train <- x[folds!=j]
    y.train <- y[folds!=j]
    
    x.tr <- (x.train - min(x.train))/(diff(range(x.train)))

    ans <- waveMesh(y.train, x.tr, lams = lams, nlam = nlam,
                    max.lam = max.lam, lam.min.ratio = lam.min.ratio, 
                    nlevel = nlevel,
                    max.iter = max.iter, tol = tol,
                    theta0 = theta0, 
                    filter.number = filter.number, j0 = j0,
                    family = family, bc = bc, ls.alpha = ls.alpha)
    # Obtain predicted values on test set.
    x.test <- x[folds==j]
    x.te <- (x.test - min(x.train))/(diff(range(x.train)))
    
    fhats <- predict.waveMesh(ans, x.tr = x.tr, new.x = x.te)
    cv.mse[j, ] <- colMeans((fhats - y[folds==j])^2)
  }
  
  cv.mse.mu <- apply(cv.mse,2,mean, na.rm = TRUE)
  cv.mse.se <- apply(cv.mse,2,function(x,...){sd(x, ...)/sqrt(length(x[!is.na(x)]))},
                     na.rm = TRUE)
  ind.min <- which.min(cv.mse.mu)
  lam.min <- lams[ind.min]
  
  ind.1se <- head(which(cv.mse.mu <= min(cv.mse.mu) + 
                          cv.mse.se[ind.min]),1)
  
  lam.1se <- lams[ind.1se]
  list("mse" = cv.mse.mu, "mse.se" = cv.mse.se, 
       "lam.min" = lam.min, "lams" = lams, 
       "lam.1se" = lam.1se, "ind.min" = ind.min, 
       "ind.1se" = ind.1se)
}

# ################ Testing my functions #####################
# set.seed(1)
# n <- 2^6
# #x <- pmin(sort(rexp(n, rate = 4)),1)
# #x <- (x-min(x))/diff(range(x))
# x <- sort(runif(n))
# #x <- (1:n)/n
# 
# f <- sin(5*x+0.3)
# f <- (x-0.5)^2
# SNR <- 10
# var.noise <- var(f)/SNR
# y <- f + rnorm(n, sd = sqrt(var.noise))
# 
# plot(x,y, col = "red")
# #
# lams <- 10^seq(log10(2), log10(0.001), length = 50)
# #
# ans <- waveMesh(y - mean(y), x, max.lam = NULL, nlam = 50, 
#                 nlevel = ceiling(log2(n)), lam.min.ratio = 1e-3,
#                 lams = lams,
#             j0 = 4,  bc = "periodic", family = "DaubLeAsymm",
#                               max.iter = 10000, tol = 1e-30)
# 
# plot(lams, colMeans((ans$ans - f)^2), log = "y")
# ind <- which.min(colMeans((ans$ans - f)^2))
# plot(x,f)
# lines(x, ans$ans[,ind], col = "red", lwd = 3)
# 
# 
# library(glmnet)
# j <- 10
# w <- GenW(n = 2^j, filter.number = 8, family = "DaubExPhase", bc = "periodic")
# #r <- genR(x, nlevel = j)
# # 
# xmat <- w
# las <- glmnet(xmat, y, intercept = FALSE, alpha = 1, lambda = ans$lams/n,
#        penalty.factor = c(0, rep(1, ncol(xmat)-1)) , standardize = FALSE,
#        thresh = 1e-30)
# 
# fhat.las <- predict(las, newx = xmat)
# 
# mse <- colMeans((f - ans$ans)^2)
# plot(mse, log = "y", type = "l")
# # 
# mse2 <- colMeans((f - predict(las, newx = xmat))^2)
# lines(mse2, col = "blue")
# fhat1 <- ans$ans
# fhat2 <- predict(las, newx= xmat)
# 
# max(abs(fhat1 - fhat2))
# 
# plot(fhat1[,10]-fhat2[,10])
# 
# 
# library(glmnet)
# w <- myw(n = n, filter.number = 8, family = "DaubExPhase", bc = "interval")
# 
# 
# 
# plot(x,y, col = "red")
# lines(x,f, lwd = 3)
# lines(x, ans$ans[,which.min(mse)], col = "blue", lwd = 3)
# min(mse)
