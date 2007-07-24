"plot.varprd" <-
function(x, ...){
  K <- ncol(x$endog)
  smpl <- nrow(x$endog)
  ynames <- colnames(x$endog)
  op <- par(no.readonly = TRUE)
  for(i in 1 : K){
    par(ask = TRUE)
    fcsty <- c(rep(NA, smpl - 1), x$endog[smpl, i], x$fcst[[i]][, 1])
    fcstl <- c(rep(NA, smpl - 1), x$endog[smpl, i], x$fcst[[i]][, 2])
    fcstu <- c(rep(NA, smpl - 1), x$endog[smpl, i], x$fcst[[i]][, 3])
    smply <- c(x$endog[, i], rep(NA, length(x$fcst[[i]][, 1])))
    min.y <- min(na.omit(c(fcsty, fcstl, fcstu, smply)))
    max.y <- max(na.omit(c(fcsty, fcstl, fcstu, smply)))               
    plot.ts(fcsty, ylab = "", xlab = "", ylim = c(min.y, max.y), main = paste("Forecast of series", ynames[i]), col = "blue", lty = 2)
    lines(smply, col = "black", lty = 1)
    lines(fcstl, col = "red", lty = 3)
    lines(fcstu, col = "red", lty = 3)
    abline(v = smpl, col = "grey", lty = 4)
  }
  on.exit(par(op))
}
