# This file consists of functions which take the competitor methods
# of Sardy et al. and Koval and Silverman (2000) and implements a  
# K-fold cross validation scheme for each method. 

# First for the interpolation scheme of Koval and Silvermal (2000).
irr.wd <- function(x, y, lams, filter.number = 8, j0 = 1,
                    family = "DaubExPhase", bc = "symmetric", ...) {
  # This function performs wavlet regression for unequal spaced data 
  # using an interpolation method.
  nlam <- length(lams)
  gd <- makegrid(t = x, y = y)
  fhat <- matrix(NA, ncol = length(lams), nrow = length(y))
  for(i in 1:nlam){
    wty <- irregwd(gd, family= family, filter.number = filter.number,
                   bc = bc)
    th.obj <- threshold(wty, levels = j0:(nlevelsWT(wty)-1), 
                        type = "soft", #policy = "cv")
                         policy = "manual", 
                         value = lams[i])
    if(bc == "interval") {
      new.val <- wr.int(th.obj)
    } else {
      new.val <- wr(th.obj)
    }
    fhat[, i] <- approx(gd$gridt, new.val, xout = x, rule = 2)$y
  }
  return(fhat)
}

# Now we write a CV for the isometric wavelet methods. 
irr.iso <- function(x, y, lams, filter.number = 8, j0 = 1,
                   family = "DaubExPhase", bc = "symmetric", ...) {
  # This function performs wavlet regression for unequal spaced data 
  # using the isometric wavelet method of Sardy et al.
  # According to paper, this amounts to doing a DWT treated the sorted 
  # data as if it were equally spaced on the grid, i.e. just ignore 
  # irregularities.
  nlam <- length(lams)
  fhat <- matrix(NA, ncol = length(lams), nrow = length(y))
  for(i in 1:nlam){
    wty <- wd(y, family= family, filter.number = filter.number,
              bc = bc)
    th.obj <- threshold(wty, levels = j0:(nlevelsWT(wty)-1), 
                        type = "soft",
                        policy = "manual", 
                        value = lams[i])
    if(bc == "interval") {
      new.val <- wr.int(th.obj)
    } else {
      new.val <- wr(th.obj)
    }
    fhat[, i] <- new.val
  }
  return(fhat)
}

predict.irr <- function(fh, x.tr, new.x) {
  apply(fh, 2, function(vec) {
    approx(x.tr, vec, xout = new.x, rule = 2)$y
  })
}


#################################################################
#################################################################
#################################################################
#################################################################

cv.irr <- function(x, y, lams, method = irr.wd, 
                   K = 4, filter.number = 8, j0 = 1,
                   family = "DaubExPhase", bc = "symmetric") {
  
  nlam <- length(lams)
  max.lam <- max(lams)
  lam.min.ratio <- min(lams)/max(lams)
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
    
    ans <- method(x.tr, y.train, lams = lams,
                  filter.number = filter.number, j0 = j0,
                  family = family, bc = bc)
    # Obtain predicted values on test set.
    x.test <- x[folds==j]
    x.te <- (x.test - min(x.train))/(diff(range(x.train)))
    
    fhats <- predict.irr(ans, x.tr = x.tr, new.x = x.te)
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
