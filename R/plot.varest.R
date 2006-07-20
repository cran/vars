"plot.varest" <-
function(x, ...){
  x.sum <- summary(x)
  K <- x$K
  for(i in 1 : K){
    layout(matrix(c(1, 1, 2, 2, 3, 4), nrow = 3, ncol = 2, byrow = TRUE))
    plot.ts(x$datamat[, i], main = paste("Diagram of fit for", colnames(x$datamat)[i], sep=" "), ylim = c(min(c(x$datamat[, i], x$varresult[[i]]$fitted.values)), max(c(x$datamat[, i], x$varresult[[i]]$fitted.values))), ylab = "", lty = 1)
    lines(x$varresult[[i]]$fitted.values, col = "blue", lty = 2)
    plot.ts(x$varresult[[i]]$residuals, main = "Residuals", ylab = "", lty = 1)
    abline(h = 0, col = "red")
    acf(x$varresult[[i]]$residuals, main = "ACF Residuals", ylab = "")
    pacf(x$varresult[[i]]$residuals, main = "PACF Residuals", ylab = "")
    if (interactive()){
      cat("\nType <Return> to continue: ")
      readline()
    }
  }
}
