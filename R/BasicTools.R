#' Title
#' FDP Calculate
#'
#' Calculate the false discovery proportion.
#' @param selected a vector indicating the indexes of the selected hypotheses.
#' @param H0 a vector indicating the underlying states of the hypotheses.
#'
#' @return a numerical value representing the false discovery proportion.
#' @export
fdp <- function(selected,H0){
  if(length(selected)==0){
    fdp <- 0
  }else{
    fdp <- sum(H0[selected]==1)/length(selected)}
  return(fdp)
}

#' Title
#' Power Calculate
#'
#' Calculate the true positive percentage.
#' @param selected a vector indicating the indexes of the selected hypotheses.
#' @param H0 a vector indicating the underlying states of the hypotheses.
#'
#' @return a numerical value representing the true positive percentage.
#' @export
Pow <- function(selected,H0){
  if(sum(1-H0)==0){
    Pow <- 0
  }else{
    Pow <- sum(H0[selected]==0)/sum(1-H0)
  }
  return(Pow)
}
