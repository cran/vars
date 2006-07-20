"plot.varcheck" <-
function(x, ...){
  K <- ncol(x$resid)
  resids <- x$resid
  for(i in 1 : K){
    layout(matrix(c(1, 2, 3, 4, 5, 6), nrow=3, ncol=2, byrow=TRUE))
    plot.ts(resids[, i], main = paste("Diagram of fit for", colnames(resids)[i], "residuals", sep=" "), ylim = c(min(resids[, i]), max(resids[, i])), ylab = "", lty = 1)
    abline(h = 0, col = "red")
    hist(resids[, i], main = "Histogram and EDF", freq = FALSE)
    lines(density(resids[, i]), col = "blue")
    acf(resids[, i], main = "ACF Residuals", ylab = "")
    pacf(resids[, i], main = "PACF Residuals", ylab = "")
    acf(resids[, i]^2, main = "ACF of squared Residuals", ylab = "")
    pacf(resids[, i]^2, main = "PACF of squared Residuals", ylab = "")
    if (interactive()){
      cat("\nType <Return> to continue: ")
      readline()
    }
  }
}
