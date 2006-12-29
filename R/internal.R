require(MASS)
require(strucchange)
## Forecast variance-covariance matrix
".fecov" <-
function(x, n.ahead) {
  sigma.u <- crossprod(x$resid)/(x$obs - ncol(x$datamat[, -c(1:x$K)]))
  Sigma.yh <- array(NA, dim = c(x$K, x$K, n.ahead))
  Sigma.yh[, , 1] <- sigma.u
  Phi <- Phi(x, nstep = n.ahead)
  if (n.ahead > 1) {
    for (i in 2:n.ahead) {
      temp <- matrix(0, nrow = x$K, ncol = x$K)
      for (j in 2:i) {
        temp <- temp + Phi[, , j] %*% sigma.u %*% t(Phi[, , j])
      }
      Sigma.yh[, , i] <- temp + Sigma.yh[, , 1]
    }
  }
  return(Sigma.yh)
}
## Forecast variance-covariance matrix (SVAR)
".fecovsvar" <-
function(x, n.ahead) {
  Sigma.yh <- array(NA, dim = c(x$var$K, x$var$K, n.ahead))
  Phi <- Phi(x, nstep = n.ahead)
  Sigma.yh[, , 1] <- Phi[, , 1]%*%t(Phi[, , 1])
  if (n.ahead > 1) {
    for (i in 2:n.ahead) {
      temp <- matrix(0, nrow = x$var$K, ncol = x$var$K)
      for (j in 2:i) {
        temp <- temp + Phi[, , j]%*%t(Phi[, , j])
      }
      Sigma.yh[, , i] <- temp + Sigma.yh[, , 1]
    }
  }
  return(Sigma.yh)
}
## irf (internal)
".irf" <-
function(x, impulse, response, y.names, n.ahead, ortho, cumulative){
  if(class(x) == "varest"){
    if(ortho){
      irf <- Psi(x, nstep = n.ahead)
    } else {
      irf <- Phi(x, nstep = n.ahead)
    }
  } else if(class(x) == "svarest"){
    irf <- Phi(x, nstep = n.ahead)
  }
  dimnames(irf) <- list(y.names, y.names, NULL)
  idx <- length(impulse)
  irs <- list()
  for(i in 1 : idx){
    irs[[i]] <- matrix(t(irf[response , impulse[i], 1 : (n.ahead + 1)]), nrow = n.ahead+1)
    colnames(irs[[i]]) <- response
    if(cumulative){
      if(length(response) > 1) irs[[i]] <- apply(irs[[i]], 2, cumsum)
      if(length(response) == 1){
        tmp <- matrix(cumsum(irs[[1]]))
        colnames(tmp) <- response
        irs[[1]] <- tmp
      }
    }
  }
  names(irs) <- impulse
  result <- irs
  return(result)    
}
## bootstrapping irf
".boot" <-
function(x, n.ahead, runs, ortho, cumulative, impulse, response, ci, seed, y.names){
  if(!(is.null(seed))) set.seed(abs(as.integer(seed)))
  ifelse(class(x) == "varest", VAR <- x, VAR <- x$var)
  p <- VAR$p
  K <- VAR$K
  obs <- VAR$obs
  total <- VAR$totobs
  type <- VAR$type
  B <- B(VAR)
  BOOT <- list()
  ysampled <- matrix(0, nrow = total, ncol = K)
  colnames(ysampled) <- colnames(VAR$y)
  Zdet <- switch(type,
                 "const" = matrix(rep(1, obs), nrow = obs, ncol = 1),
                 "trend" = matrix(seq(1 : obs), nrow = obs, ncol = 1),
                 "both" = matrix(rep(1, obs), seq(1 : obs), nrow = obs, ncol = 2),
                 "none" = NULL)
  resorig <- scale(VAR$resid, scale = FALSE)
  B <- B(VAR)
  for(i in 1:runs){
    booted <- sample(c(1 : obs), replace=TRUE)
    resid <- resorig[booted, ]
    lasty <- c(t(VAR$y[p : 1, ]))
    ysampled[c(1 : p), ] <- VAR$y[c(1 : p), ]
    for(j in 1 : obs){
      lasty <- lasty[1 : (K * p)]
      Z <- c(Zdet[j, ], lasty)
      ysampled[j + p, ] <- B %*% Z + resid[j, ]
      lasty <- c(ysampled[j + p, ], lasty) 
    }
    varboot <- update(VAR, y = ysampled)
    if(class(x) == "svarest"){
      varboot <- update(x, x = varboot)
    }
    BOOT[[i]] <- .irf(x = varboot, n.ahead = n.ahead, ortho = ortho, cumulative = cumulative, impulse = impulse, response = response, y.names=y.names)
  }
  lower <- ci / 2
  upper <- 1 - ci / 2
  mat.l <- matrix(NA, nrow = n.ahead + 1, ncol = length(response))
  mat.u <- matrix(NA, nrow = n.ahead + 1, ncol = length(response))
  Lower <- list()
  Upper <- list()
  idx1 <- length(impulse)
  idx2 <- length(response)
  idx3 <- n.ahead + 1
  temp <- rep(NA, runs)
  for(j in 1 : idx1){
    for(m in 1 : idx2){
      for(l in 1 : idx3){
        for(i in 1 : runs){
          if(idx2 > 1){
            temp[i] <- BOOT[[i]][[j]][l, m]
          } else {
            temp[i] <- matrix(BOOT[[i]][[j]])[l, m]
          }
        }
        mat.l[l, m] <- quantile(temp, lower)
        mat.u[l, m] <- quantile(temp, upper)
      }
    }
    colnames(mat.l) <- response
    colnames(mat.u) <- response
    Lower[[j]] <- mat.l
    Upper[[j]] <- mat.u
  }
  names(Lower) <- impulse
  names(Upper) <- impulse
  result <- list(Lower = Lower, Upper = Upper)
  return(result)
}
## Duplication matrix
".duplicate" <-
function(n){
  D <- matrix(0, nrow = n^2, ncol = n * (n + 1) / 2)
  count <- 0
  for(j in 1 : n){
    D[(j - 1) * n + j, count + j] <- 1
    if((j + 1) <= n){
      for(i in (j + 1):n){
        D[(j - 1) * n + i, count + i] <- 1
        D[(i - 1) * n + j, count + i] <- 1
      }
    }
    count <- count + n - j
  }
  return(D)
}
## univariate ARCH test
".arch.uni" <-
function(x, lags.single){
  lags.single <- lags.single + 1
  mat <- embed(scale(x)^2, lags.single)
  arch.lm <- summary(lm(mat[, 1] ~ mat[, -1]))
  STATISTIC <- arch.lm$r.squared*length(resid(arch.lm))
  names(STATISTIC) <- "Chi^2"
  PARAMETER <- lags.single - 1
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "ARCH test (univariate)"
  result <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = deparse(substitute(x)))
  class(result) <- "htest"
  return(result)
}
## multivariate ARCH test
".arch.multi" <-
function(x, lags.multi, obj.name, K, obs){
  col.arch.df <- 0.5 * K * (K + 1)
  arch.df <- matrix(NA, nrow = obs, ncol = col.arch.df)
  for( i in 1 : obs){
    temp <- outer(x[i,], x[i,])
    arch.df[i,] <- temp[lower.tri(temp, diag=TRUE)]
  }
  lags.multi <- lags.multi + 1
  arch.df <- embed(arch.df, lags.multi)
  archm.lm0 <- lm(arch.df[ , 1:col.arch.df] ~ 1)
  archm.lm0.resids <- resid(archm.lm0)
  omega0 <- cov(archm.lm0.resids)
  archm.lm1 <- lm(arch.df[ , 1 : col.arch.df] ~ arch.df[ , -(1 : col.arch.df)])
  archm.lm1.resids <- resid(archm.lm1)
  omega1 <- cov(archm.lm1.resids)
  R2m <- 1 - (2 / (K * (K + 1))) * sum(diag(omega1 %*% solve(omega0)))
  n <- nrow(archm.lm1.resids)
  STATISTIC <- 0.5 * n * K * (K+1) * R2m
  names(STATISTIC) <- "Chi^2"
  lags.multi <- lags.multi - 1
  PARAMETER <- lags.multi * K^2 * (K + 1)^2 / 4
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "ARCH (multivariate)"
  result <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object", obj.name))
  class(result) <- "htest"
  return(result)
}
## univariate normality test
".jb.uni" <-
function(x, obs){
  x <- as.vector(x)
  m1 <- sum(x) / obs
  m2 <- sum((x - m1)^2) / obs
  m3 <- sum((x - m1)^3) / obs
  m4 <- sum((x - m1)^4) / obs
  b1 <- (m3 / m2^(3 / 2))^2
  b2 <- (m4/m2^2)
  STATISTIC <- obs * b1 / 6 + obs * (b2 - 3)^2 / 24
  names(STATISTIC) <- "Chi^2"
  PARAMETER <- 2
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = 2)
  METHOD <- "JB-Test (univariate)"
  result <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = deparse(substitute(x)))
  class(result) <- "htest"
  return(result)
}
## multivariate normality test
".jb.multi" <-
function(x, obs, K, obj.name){
  P <- chol(crossprod(x) / obs)
  resids.std <- x %*% solve(P)
  b1 <- apply(resids.std, 2, function(x) sum(x^3) / obs)
  b2 <- apply(resids.std, 2, function(x) sum(x^4) / obs)
  s3 <- obs * t(b1) %*% b1 / 6
  s4 <- obs * t(b2 - rep(3, K)) %*% (b2 - rep(3, K)) / 24
  STATISTIC <- s3 + s4
  names(STATISTIC) <- "Chi^2"
  PARAMETER <- 2 * K
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "JB-Test (multivariate)"
  result1 <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object", obj.name))
  class(result1) <- "htest"
  STATISTIC <- s3
  names(STATISTIC) <- "Chi^2"
  PARAMETER <- K
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "Skewness only (multivariate)"
  result2 <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object", obj.name))
  class(result2) <- "htest"
  STATISTIC <- s4
  names(STATISTIC) <- "Chi^2"
  PARAMETER <- K
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "Kurtosis only (multivariate)"
  result3 <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object", obj.name))
  class(result3) <- "htest"
  result <- list(JB = result1, Skewness = result2, Kurtosis = result3)
  return(result)
}
## Convenience function for computing lagged x
".matlag1" <-
function(x, lag = 1){
  totcols <- ncol(x)
  nas <- matrix(NA, nrow = lag, ncol = totcols)
  x <- rbind(nas, x)
  totrows <- nrow(x)
  x <- x[-c((totrows - lag + 1) : totrows), ]
  return(x)
}
## Multivariate Portmanteau Statistik
".pt.multi" <-
function(x, K, obs, lags.pt, obj.name, resids){
  C0 <- crossprod(resids) / obs
  C0inv <- solve(C0)
  tracesum <- rep(NA, lags.pt)
  for(i in 1 : lags.pt){
    Ut.minus.i <- .matlag1(resids, lag = i)[-c(1 : i), ]
    Ut <- resids[-c(1 : i), ]
    Ci <- crossprod(Ut, Ut.minus.i) / obs
    tracesum[i] <- sum(diag(t(Ci) %*% C0inv %*% Ci %*% C0inv))
  }
  vec.adj <- obs - (1 : lags.pt)
  Qh <- obs * sum(tracesum)
  Qh.star <- obs^2 * sum(tracesum / vec.adj)
  nstar <- K^2 * x$p
  ## htest objects for Qh and Qh.star
  STATISTIC <- Qh
  names(STATISTIC) <- "Chi^2"
  PARAMETER <- (K^2 * lags.pt - nstar)
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "Portmanteau Test (asymptotic)"
  PT1 <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object", obj.name))
  class(PT1) <- "htest"
  STATISTIC <- Qh.star
  names(STATISTIC) <- "Chi^2"
  PARAMETER <- (K^2 * lags.pt - nstar)
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "Portmanteau Test (adjusted)"
  PT2 <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object", obj.name))
  class(PT2) <- "htest"
  result <- list(PT1 = PT1, PT2 = PT2)
  return(result)
}
## Convenience function for computing lagged residuals
".matlag2" <-
function(x, lag = 1){
  K <- ncol(x)
  obs <- nrow(x)
  zeromat <- matrix(0, nrow = obs, ncol = K * lag)
  idx1 <- seq(1, K * lag, K)
  idx2 <- seq(K, K * lag, K)
  for(i in 1:lag){
    lag <- i + 1
    res.tmp <- embed(x, lag)[, -c(1 : (K * i))]
    zeromat[-c(1 : i), idx1[i] : idx2[i]] <- res.tmp
  }
  resids.l <- zeromat
  return(resids.l)
}
## Breusch-Godfrey and Edgerton-Shukur Test
".bgserial" <-
function(x, K, obs, lags.bg, obj.name, resids){
  ylagged <- as.matrix(x$datamat[, -c(1 : K)])
  resids.l <- .matlag2(resids, lag = lags.bg)
  if(is.null(x$restrictions)){
    regressors <- as.matrix(cbind(ylagged, resids.l))
    lm0 <- lm(resids ~ -1 + regressors)
    lm1 <- lm(resids ~ -1 + ylagged)
    sigma.1 <- crossprod(resid(lm1)) / obs
    sigma.0 <- crossprod(resid(lm0)) / obs
  } else {
    resid0 <- matrix(NA, ncol = K, nrow = obs)
    resid1 <- matrix(NA, ncol = K, nrow = obs)
    for(i in 1 : K){
      datares <- data.frame(ylagged[, which(x$restrictions[i, ] == 1)])
      regressors <- data.frame(cbind(datares, resids.l))
      lm0 <- lm(resids[, i] ~ -1 + ., data=regressors)
      lm1 <- lm(resids[, i] ~ -1 + ., data=datares)
      resid0[, i] <- resid(lm0)
      resid1[, i] <- resid(lm1)
      sigma.0 <- crossprod(resid0) / obs
      sigma.1 <- crossprod(resid1) / obs
    }
  }
  LMh.stat <- obs * (K - sum(diag(crossprod(solve(sigma.1), sigma.0))))
  STATISTIC <- LMh.stat
  names(STATISTIC) <- "Chi^2"
  PARAMETER <- lags.bg * K^2
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "Breusch-Godfrey LM test"
  LMh <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object", obj.name))
  class(LMh) <- "htest"
  ## small sample correction of Edgerton Shukur
  R2r <- 1 - det(sigma.0) / det(sigma.1)
  m <- K * lags.bg
  q <- 0.5 * K * m - 1
  n <- ncol(x$datamat) - K
  N <- obs - n - m - 0.5 * (K - m + 1)
  r <- sqrt((K^2 * m^2 - 4)/(K^2 + m^2 - 5))
  LMFh.stat <- (1 - (1 - R2r)^(1 / r))/(1 - R2r)^(1 / r) * (N * r - q) / (K * m)
  STATISTIC <- LMFh.stat
  names(STATISTIC) <- "F statistic"
  PARAMETER1 <- lags.bg * K^2
  names(PARAMETER1) <- "df1"
  PARAMETER2 <- floor(N * r - q)
  names(PARAMETER2) <- "df2"
  PVAL <-   1 - pf(LMFh.stat, PARAMETER1, PARAMETER2)
  METHOD <- "Edgerton-Shukur F test"
  LMFh <- list(statistic = STATISTIC, parameter = c(PARAMETER1, PARAMETER2), p.value = PVAL, method = METHOD, data.name = paste("Residuals of VAR object", obj.name))
  class(LMFh) <- "htest"
  return(list(LMh = LMh, LMFh = LMFh))
}

