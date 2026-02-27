#include "JointMirror.h"


// [[Rcpp::depends(RcppArmadillo)]]
RCPP_MODULE(JOINTMIRRORMODE){
  class_<JointMirror>("JointMirror")
  .constructor<arma::mat, arma::mat>()
  .method("print", &JointMirror::print)
  .method("InitPara", &JointMirror::InitPara)
  .method("runJM", &JointMirror::runJM)
  .method("getMaskInd", &JointMirror::getMaskInd)
  .method("getRejInd", &JointMirror::getRejInd)
  .method("getRootInd", &JointMirror::getRootInd)
  .method("getFDPest", &JointMirror::getFDPest)
  .method("getProbInRej", &JointMirror::getProbInRej)
  .method("getunMaskInd", &JointMirror::getunMaskInd)
  .method("excludeInd", &JointMirror::getexcludeInd)
  .method("getTime", &JointMirror::getTime)
  ;

}
