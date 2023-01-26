
#' @export
fdp <- function(selected,H0){
  if(length(selected)==0){
    fdp <- 0
  }else{
    fdp <- sum(H0[selected]==1)/length(selected)}
  return(fdp)
}

#' @export
Pow <- function(selected,H0){
  if(sum(1-H0)==0){
    Pow <- 0
  }else{
    Pow <- sum(H0[selected]==0)/sum(1-H0)
  }
  return(Pow)
}
