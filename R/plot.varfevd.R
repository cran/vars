"plot.varfevd" <-
function(x, ...){
  K <- length(x)
  ynames <- names(x)
  op <- par(no.readonly = TRUE)
  for(i in 1 : K){
    barplot(t(x[[i]]), main = paste("FEVD for", ynames[i]), col = palette()[1 : K], ylab = "Percentage", xlab = "Horizon", names.arg = paste(1 : nrow(x[[i]])), ylim = c(0, 1.2))
    legend("top", legend = ynames, fill = palette()[1 : K], ncol = K)
    if (interactive()){
      cat("\nType <Return> to continue: ")
      readline()
    }
  }
  par(op)
}
