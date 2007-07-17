"plot.varirf" <-
function(x, ...){
  idx1 <- length(x$impulse)
  idx2 <- length(x$response)
  op <- par(no.readonly = TRUE)
  for(i in 1 : idx1){
    layout(matrix(1 : idx2, nrow = idx2, ncol = 1))
    for(j in 1 : idx2){
      limit1 <- min(c(x$irf[[i]][, j], x$Lower[[i]][, j], x$Upper[[i]][, j]))
      limit2 <- max(c(x$irf[[i]][, j], x$Lower[[i]][, j], x$Upper[[i]][, j]))
      if((x$model == "varest") || (x$model == "vec2var")){
        if(x$ortho){
          text <- paste("Orthogonal Impulse Response from", x$impulse[i], "to", x$response[j], sep = " ")
        } else {
          text <- paste("Forecast Error Impulse Response from", x$impulse[i], "to", x$response[j], sep = " ")
        }
      } else if(x$model == "svarest"){
        text <- paste("SVAR Impulse Response from", x$impulse[i], "to", x$response[j], sep = " ")
      }
      else if (x$model == "svecest") {
        text <- paste("SVECM Impulse Response from", x$impulse[i], "to", x$response[j], sep = " ")
      }
      if(x$cumulative) text <- paste(text, "(cumulative)", sep=" ")
      plot.ts(x$irf[[i]][, j], ylab = "", xlab = "", ylim = c(limit1, limit2), main = text)
      abline(h = 0, col = "gray")
      if(x$boot){
        lines(x$Lower[[i]][, j], col = "red", lty = 2)
        lines(x$Upper[[i]][, j], col = "red", lty = 2)
        mtext(paste((1-x$ci)*100, "% Bootstrap CI, ", x$runs, "runs"), side = 1, line = 2, outer = FALSE)
      }
    }
    if(idx1 > 1){
      if (interactive()){
        cat("\nType <Return> to continue: ")
        readline()
      }
    }
  }
  par(op)
}
