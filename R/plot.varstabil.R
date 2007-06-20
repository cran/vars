"plot.varstabil" <-
function(x, ...){
  K <- x$K
  op <- par(no.readonly = TRUE)
  for(i in 1 : K){
    title <- paste(x[[1]][[i]]$type, "of equation", x[[2]][i])
    plot(x[[1]][[i]], main = title)
    if (interactive()){
      cat("\nType <Return> to continue: ")
      readline()
    }
  }
  par(op)
}
