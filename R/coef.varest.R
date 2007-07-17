"coef.varest" <-
function(object, ...){
  return(sapply(object$varresult, coef))
}
