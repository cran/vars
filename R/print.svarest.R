"print.svarest" <-
function(x, ...){
  cat("\nEstimated A matrix:\n")
  print(x$A)
  cat("\nEstimated standard errors for A matrix:\n")  
  print(x$Ase)
  cat("\nEstimated B matrix:\n")
  print(x$B)
  cat("\nEstimated standard errors for B matrix:\n")  
  print(x$Bse)
  cat("\nLR overidentification test:\n")
  print(x$LR)
  invisible(x)                                                               
}
