"predict" <-
function(x, ..., n.ahead, ci, dumvar){
  UseMethod("predict", x)
}
