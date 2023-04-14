
//#include "GadditiveModel_FullRank.hpp"
#include "JointMirror.h"


// [[Rcpp::depends(RcppArmadillo)]]
//RCPP_EXPOSED_CLASS(JOINTMIRROR)
RCPP_MODULE(JOINTMIRROR){
  class_<JointMirror>("JointMirror")
  .constructor< mat, mat>()
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

