#' JMirror: A Package for Joint Mirror Procedure
#'
#' The JMirror package provides two important functions for implementing joint mirror procedure: JointMirror.R and JointMirror.Qvalue.
#'
#'
#'
#'
#' @name JMirror
#' @useDynLib JMirror, .registration=TRUE
#' @importFrom Rcpp loadModule
#' @rawNamespace export(JointMirror)
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  Rcpp::loadModule("JOINTMIRRORMODE", what = TRUE)
}
