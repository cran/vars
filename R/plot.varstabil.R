"plot.varstabil" <-
function(x, ...){
  K <- x$K
  op <- par(no.readonly = TRUE)
  for(i in 1 : K){
    par(ask = TRUE)
    title <- paste(x[[1]][[i]]$type, "of equation", x[[2]][i])
    plot(x[[1]][[i]], main = title)
  }
  on.exit(par(op))
}
