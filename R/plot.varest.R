"plot.varest" <-
function(x, ...){
  K <- x$K
  resids <- resid(x)
  fitted <- fitted(x)
  op <- par(no.readonly = TRUE)
  for(i in 1 : K){
    par(ask = TRUE)
    layout(matrix(c(1, 1, 2, 2, 3, 4), nrow = 3, ncol = 2, byrow = TRUE))
    plot.ts(x$datamat[, i], main = paste("Diagram of fit for", colnames(x$datamat)[i], sep=" "), ylim = c(min(c(x$datamat[, i], fitted[, i])), max(c(x$datamat[, i], fitted[, i]))), ylab = "", lty = 1)
    lines(fitted[, i], col = "blue", lty = 2)
    plot.ts(resids[, i], main = "Residuals", ylab = "", lty = 1)
    abline(h = 0, col = "red")
    acf(resids[, i], main = "ACF Residuals", ylab = "")
    pacf(resids[, i], main = "PACF Residuals", ylab = "")
  }
  on.exit(par(op)) 
}
