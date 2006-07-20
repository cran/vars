"VARselect" <-
function(y, lag.max=10, type = c("const", "trend", "both", "none")){
  y <- as.matrix(y)
  if (any(is.na(y)))
    stop("\nNAs in y.\n")
  colnames(y) <- make.names(colnames(y))
  K <- ncol(y)  
  lag.max <- abs(as.integer(lag.max))
  type <- match.arg(type)
  lag <- lag.max + 1
  ylagged <- embed(y, lag)[, -c(1 : K)]
  yendog <- y[-c(1 : lag.max), ]
  sample <- nrow(ylagged)
  rhs <- switch(type,
                "const" = rep(1, sample),
                "trend" = seq(1, sample),
                "both" = cbind(rep(1, sample), seq(1, sample)),
                "none" = NULL)
  idx <- seq(K, K*lag.max, K)
  criteria <- matrix(NA, nrow = 4, ncol = lag.max)
  rownames(criteria) <- c("AIC(n)", "HQ(n)", "SC(n)", "FPE(n)")
  colnames(criteria) <- paste(seq(1 : lag.max))
  for(i in 1 : lag.max){
    ys.lagged <- cbind(ylagged[, c(1 : idx[i])], rhs)
    sampletot <- nrow(y)
    nstar <- ncol(ys.lagged)
    resids <- resid(lm(yendog ~ -1 + ys.lagged))
    sigma.det <- det(crossprod(resids)/sample)
    criteria[1, i] <- log(sigma.det) + (2 / sample) * i * K^2 
    criteria[2, i] <- log(sigma.det) + (2 * log(log(sample)) / sample) * i * K^2
    criteria[3, i] <- log(sigma.det) + (log(sample) / sample)* i* K^2
    criteria[4, i] <- ((sample + nstar)/(sample - nstar))^K * sigma.det    
  }
  order <- apply(criteria, 1, which.min)
  return(list(selection=order, criteria=criteria))
}
