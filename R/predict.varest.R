"predict.varest" <-
function(object, ..., n.ahead = 10, ci = 0.95){
  K <- object$K
  p <- object$p
  obs <- object$obs
  type <- object$type
  data.all <- object$datamat
  ynames <- colnames(object$y)
  Z <- object$datamat[, -c(1 : K)]
  B <- B(object)
  ## Deterministic and lagged y's
  ## Retrieval of A in matrix (whole)
  if(type == "const"){
      Zdet <- matrix(rep(1, n.ahead), nrow = n.ahead, ncol = 1)
      Zy <- as.matrix(data.all[, -c(K + 1)])
    }else if(type == "trend"){
      trdstart <- Z[nrow(Z), 1] + 1 
      Zdet <- matrix(seq(trdstart, length = n.ahead), nrow = n.ahead, ncol = 1)
      Zy <- as.matrix(data.all[, -c(K + 1)])
    }else if(type == "both"){
      trdstart <- Z[nrow(Z), 2] + 1
      Zdet <- matrix(c(rep(1, n.ahead), seq(trdstart, length = n.ahead)), nrow = n.ahead, ncol = 2)
      Zy <- as.matrix(data.all[, -c(K + 1, K + 2)])
    }else if(type == "none"){
      Zdet <- NULL
      Zy <- as.matrix(data.all)
    }
  yse <- matrix(NA, nrow = n.ahead, ncol = K)
  sig.y <- .fecov(x = object, n.ahead = n.ahead)
  for(i in 1 : n.ahead){
    yse[i, ] <- sqrt(diag(sig.y[, , i]))
  }
  yse <- -1 * qnorm((1 - ci) / 2) * yse
  colnames(yse) <- paste(ci, "of", ynames)
  ## forecast recursion
  forecast <- matrix(NA, ncol = K, nrow = n.ahead)
  lasty <- c(Zy[nrow(Zy), ])
  for(i in 1 : n.ahead){
    lasty <- lasty[1 : (K * p)]
    Z <- c(Zdet[i, ], lasty)
    forecast[i, ] <- B %*% Z
    temp <- forecast[i, ]
    lasty <- c(temp, lasty)
  }
  colnames(forecast) <- paste(ynames, ".fcst", sep="")
  lower <- forecast - yse
  colnames(lower) <- paste(ynames, ".lower", sep="")
  upper <- forecast + yse
  colnames(upper) <- paste(ynames, ".upper", sep="")
  forecasts <- list()
  for(i in 1 : K){
    forecasts[[i]] <- cbind(forecast[, i], lower[, i], upper[, i], yse[, i])
    colnames(forecasts[[i]]) <- c("fcst", "lower", "upper", "CI")
  }
  names(forecasts) <- ynames
  result <- list(fcst = forecasts, endog = object$y, model = object) 
  class(result) <- "varprd"
  return(result)
}
