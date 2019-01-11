source("waveMesh.R")

# This function fits a sparse additive model for one
# value of the tuning parameters for the case of a least squares loss.
# This is done via a block-coordinate descent method. This function is
# used when loss="gaussian".
#
# Args:
#   y: Response vector of size n.
#   x: The design matrix of size n*p.
#   ord: A n*p matrix of the ordering for the covariate matrix x.
#   lambda1, lambda2: Tuning parameters.
#   init_fhat: Initial value of the estimated functions.
#   k: order of the Trend filtering problem, ignored for sobolev.
#   max_iter: Maximum iterations of the algorithm.
#   tol: Tolerance for stopping criteria.
#   method: Charecter vector indicating method to use. 
#           Redundant, we use waveMesh only.
#
# Returns:
#   A list with the following components:
#   fhat: The estimate of the component functions, an n*p matrix.
#   intercept: The estimate of the intercept term, a scalar.
#   conv: A boolean indicator of convergence status.

blockCoord_one <- function(y, x, init_fhat,
                           max_iter_outer = 100, tol_outer = 1e-4,
                           lambda1, lambda2,...) {
  
  n <- length(y)
  p <- ncol(x)
  
  counter <- 1
  converged <- FALSE
  
  old.fhat <- init_fhat
  
  # Begin loop for block coordinate descent
  while(counter < max_iter_outer & !converged) {
    for(j in 1:p) {
      res <- y - apply(init_fhat[, -j], 1, sum)- mean(y)
      init_fhat[, j] <- solve.prox.waveMesh(res - mean(res), x[,j], 
                                            lambda1, lambda2, ...)
      # Center each component to have mean 0.
      #init_fhat[, j] <- init_fhat[, j] - mean(init_fhat[, j])
    }
    
    temp_res1 <- mean((init_fhat - old.fhat)^2)
    temp_res2 <- mean((old.fhat)^2)
    
    #cat("Iteration: ", counter,", Tol:", temp_res1/(temp_res2+1e-30),"\n")
    # Compare the relative change in parameter values and compare to
    # tolerance to check convergence. Else update iteration.
    #cat("counter: ", counter, ", Del: ", temp_res1/(temp_res2+1e-30), "\n")
    if(temp_res1/(temp_res2+1e-30) <= tol_outer) {
      converged <- TRUE
    } else {
      old.fhat <- init_fhat
      counter <- counter+1
    }
  }
  list("fhat" = init_fhat,
       "intercept" = mean(y),
       "conv" = converged)
}


# This function solves the full problem for a sequence of lambda
# values.
#
# Args:
#   y: Response vector of size n.
#   x: An n*p design matrix.
#   max.iter: Maximum number of iterations for the algorithm
#   tol: The tolerance for the algorithm
#   initpars: Initial parameter values, defaults to NULL f^hat = 0.
#   lambda.max: maximum lambda value
#   lambda.min.ratio: Ratio between largest and smallest lambda value.
#   nlam: Number of lambda values.
#   loss: "Gaussian" for least squares norm, and "binomial" for logistic loss.
#   ininitintercept: Initial value for the intercept term.
#   coord.desc: A boolean, indicating if block coordinate descent
#               should be used for loss=="gaussian"
#   gamma: For alternate parametrization, i.e. lambda_1 = lambda*gamma 
#          and lambda2 = lambda(1-gamma). Gamma = 0 fits non-sparse additive model. 

fit.additive <- function(y, x, max.iter.outer = 100, tol.outer = 1e-4,
                         initpars = NULL, lambda.max = 1, lambda.min.ratio = 1e-3,
                         nlam = 50, lam.seq = NULL, loss = "gaussian",
                         initintercept = NULL,
                         coord.desc = TRUE,
                         gamma = 0,...) {
  
  # Initialize some values.
  x <- as.matrix(x)
  n <- length(y)
  p <- ncol(x)
  
  if(is.null(lam.seq)) {
    lam.seq <- seq(log10(lambda.max), log10(lambda.max*lambda.min.ratio),
                   length = nlam)
    lam.seq <- 10^lam.seq
  }
  
  nlam <- length(lam.seq)
  
  if(is.null(initpars)) {
    initpars <- matrix(0, ncol = p, nrow = n)
  }
  
  if(is.null(initintercept) & loss == "binomial") {
    mp <- mean(y)
    initintercept <- log(mp/(1-mp))
  }
  if(is.null(initintercept) & loss == "gaussian") {
    initintercept <- mean(y)
  }
  
  if(loss == "binomial"){
    y[y==0] <- -1
  }
  
  if(is.null(gamma)) {
    lam1 <- lam.seq
    lam2 <- lam.seq^2
  } else {
    lam1 <- gamma * lam.seq
    lam2 <- (1-gamma) * lam.seq
  }
  
  ans <- array(0, dim = c(n, p, nlam))
  ans.inters <- numeric(nlam)
  conv.vec <- c()
  
  # print(loss)
  # print(coord.desc)
  for(i in 1:nlam) {
    cat("Lambda value: ", i, "\n")
    if(loss=="gaussian" & coord.desc) {
      temp <- blockCoord_one(y, x, init_fhat = initpars,
                             max_iter_outer = max.iter.outer, 
                             tol_outer = tol.outer,
                             lambda1 = lam1[i],
                             lambda2 = lam2[i], ...)
    } else {
      #print(i)
      temp <- proxGrad_one(y, x_ord, lam1[i], lam2[i],
                           init_fhat = initpars, init_intercept = initintercept,
                            max_iter = max.iter, tol = tol,
                           step_size = step, alpha = alpha,
                           loss = loss, ...)
    }
    
    ans[, , i] <- temp$fhat
    ans.inters[i] <- temp$intercept
    
    initintercept <- ans.inters[i]
    initpars <- ans[, , i]
    conv.vec <- c(conv.vec, temp$conv)
  }
  
  obj <- list("f_hat" = ans,
              "intercept" = ans.inters,
              "x" = x,
              "lam" = lam.seq,
              "gamma" = gamma,
              "loss" = loss,
              "conv" = conv.vec)
  
  class(obj) <- "add_mod"
  return(obj)
}

# The Generic predict function our framework based on linear
# interpolation.
#
# Args:
#   obj: A fitted additive model, object of type "add_mod"
#   new.data: A new x_matrix which we wish to fit.
#   type: "function" will return the fitted components f_1,...,f_p.

predict.add_mod <- function(obj, new.data, type = "function") {
  nlam <- length(obj$lam)
  p <- dim(obj$f_hat)[2]
  
  ans <- array(0, dim = c(dim(new.data), nlam) )
  for(i in 1:nlam) {
    for(j in 1:p) {
      ans[,j,i] <- approx(obj$x[, j], obj$f_hat[,j,i], new.data[, j],
                          rule = 2)$y
    }
  }
  
  if(type == "function") {
    return(ans)
  } else {
    if(obj$loss == "gaussian"){
      temp <- apply(ans, c(1,3),sum)
      temp <- t(apply(temp, 1, "+", obj$intercept))
      
      return(temp)
    } else {
      temp <- apply(ans, c(1,3),sum)
      temp <- t(apply(temp, 1, "+", obj$intercept))
      
      return(1/(1 + exp(-temp)))
    }
  }
  return(ans)
}

# A function to do cross validation. 
# By default it does 5-fold cross validation.
cv.addWaveMesh <- function(y, x, K = 5, max.iter.outer = 100, tol.outer = 1e-4,
                           initpars = NULL, lambda.max = 1, lambda.min.ratio = 1e-3,
                           nlam = 50, lam.seq = NULL, loss = "gaussian",
                           initintercept = NULL,
                           coord.desc = TRUE,
                           gamma = 0, ...) {

  if(is.null(lam.seq)) {
    lam.seq <- seq(log10(lambda.max), log10(lambda.max*lambda.min.ratio),
                   length = nlam)
    lam.seq <- 10^lam.seq
  }
  
  n <- length(y)
  folds <- cut(1:n, breaks = 5, labels = F)
  cv.mse <- matrix(NA, ncol = nlam, nrow = K)
  for(j in 1:K) {
    cat("CV Fold:", j,"\n")
    x.train <- x[folds!=j, ]
    y.train <- y[folds!=j]
    
    x.mins <- apply(x.train, 2, min)
    x.ranges <- apply(x.train, 2, function(x){diff(range(x))})
    x.tr <- apply(x.train, 1, "-", x.mins)
    x.tr <- t(apply(x.tr, 2, "/", x.ranges))
    
#    ans <- fit.additive(y.train, x.tr)
    ans <- fit.additive(y.train, x.tr, max.iter.outer = max.iter.outer, 
                        tol.outer = tol.outer,
                        initpars = initpars, lambda.max = lambda.max,
                        lambda.min.ratio = lambda.min.ratio,
                        nlam = nlam, lam.seq = lam.seq, loss = loss,
                        initintercept = initintercept,
                        coord.desc = coord.desc,
                        gamma = gamma,...)
    # Obtain predicted values on test set.
    x.test <- x[folds==j, ]
    x.te <- apply(x.test, 1, "-", x.mins)
    x.te <- t(apply(x.te, 2, "/", x.ranges))
    
    fhats <- predict.add_mod(ans, new.data = x.te,type = "response")
    cv.mse[j, ] <- colMeans((fhats - y[folds==j])^2)
  }
    
  cv.mse.mu <- apply(cv.mse,2,mean, na.rm = TRUE)
  cv.mse.se <- apply(cv.mse,2,function(x,...){sd(x, ...)/sqrt(length(x[!is.na(x)]))}, na.rm = TRUE)
  ind.min <- which.min(cv.mse.mu)
  lam.min <- lam.seq[ind.min]
  
  ind.1se <- head(which(cv.mse.mu <= min(cv.mse.mu) + 
                          cv.mse.se[ind.min]),1)
  
  lam.1se <- lam.seq[ind.1se]
  list("mse" = cv.mse.mu, "mse.se" = cv.mse.se, 
       "lam.min" = lam.min, "lams" = lam.seq, 
       "lam.1se" = lam.1se)
}
