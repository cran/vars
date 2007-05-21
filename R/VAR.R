"VAR" <-
function(y, p = 1, type = c("const", "trend", "both", "none")){
  y <- as.matrix(y)
  if (any(is.na(y)))
    stop("\nNAs in y.\n")
  if(ncol(y) < 2)
    stop("The matrix 'y' should contain at least two variables. For univariate analysis consider ar() and arima() in package stats.\n")
  if(is.null(colnames(y))){
    colnames(y) <- paste("y", 1:ncol(y), sep = "")
    warning(paste("No columne names supplied in y, using:", paste(colnames(y), collapse = ", "), ", instead.\n"))
  }
  colnames(y) <- make.names(colnames(y))
  y.orig <- y
  type <- match.arg(type)
  obs <- dim(y)[1]
  K <- dim(y)[2]
  sample <- obs - p
  ylags <- embed(y, dimension = p + 1)[ , -(1:K)]
  temp1 <- NULL
  for (i in 1:p) {
    temp <- paste(colnames(y), ".l", i, sep = "")
    temp1 <- c(temp1, temp)
  }
  colnames(ylags) <- temp1
  yend <- y[-c(1:p), ]
  if(type == "const"){
    rhs <- cbind(rep(1, sample), ylags)
    colnames(rhs) <- c("const", colnames(ylags))
  }else if(type == "trend"){
    rhs <- cbind(seq(1, sample), ylags)
    colnames(rhs) <- c("trend", colnames(ylags))
  }else if(type == "both"){
    rhs <- cbind(rep(1, sample), seq(1, sample), ylags)
    colnames(rhs) <- c("const", "trend", colnames(ylags))
  }else if(type == "none"){
    rhs <- ylags
    colnames(rhs) <- colnames(ylags)
  }
  datamat <- as.data.frame(rhs)
  colnames(datamat) <- colnames(rhs)
  resid <- matrix(NA, nrow=nrow(datamat), ncol=K)
  equation <- list()
  for(i in 1:K){
    y <- yend[, i]
    equation[[colnames(yend)[i]]] <- lm(y ~ -1 + ., data=datamat)
    resid[, i] <- resid(equation[[i]])
  }
  colnames(resid) <- paste("resids of", colnames(yend))
  result <- list(varresult=equation, resid=resid, datamat=data.frame(cbind(yend, rhs)), y=y.orig, type=type, p=p, K=K, obs=sample, totobs=sample+p, restrictions=NULL, call=match.call())
  class(result) <- "varest"
  return(result)
}

