"print.svarest" <-
function(x, ...){
  title <- paste("# SVAR:", x$type, "#", sep=" ")
  row <- paste(rep("#", nchar(title)), collapse="")
  cat("\n")
  cat(row, "\n")
  cat(title, "\n")
  cat(row, "\n")
  cat("\n")
  if(identical(x$type, "Blanchard-Quah")){
    cat("\nEstimated contemporaneous impact matrix:\n")
    print(x$B)
    cat("\nEstimated identified long run impact matrix:\n")
    print(x$LRIM)
    invisible(x)
  } else {
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
}
