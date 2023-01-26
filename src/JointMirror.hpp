#ifndef FUNCTIONAL_JOINTMIRROR_ALGORITHM
#define FUNCTIONAL_JOINTMIRROR_ALGORITHM

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include <vector>
//#include "OptionDefine.cpp"

using namespace Rcpp;
using namespace arma;
using namespace std;

class JointMirror{
public:
  JointMirror( mat Pval_, mat HMat_);
  void InitPara(double offset_, double fdr_level_,
                double init_p_cut_,
                int rank_Mode_);
  
  void print();
  
  void runJM();
  vec getProbInRej(){
    return(ProbInRej);
  };
  uvec getRejInd(){
    return(find(maskLogi==1 && RejSide==1));
  };
  
  uvec getMaskInd(){
    return(find(maskLogi==1));
  };
  
  uvec getunMaskInd(){
    return(unmaskInd);
  };
  
  uvec getexcludeInd(){
    return(excludeInd);
  };
  
  IntegerVector getRootInd(){
    return(rootInd);
  };
  
  vec getFDPest(){
    return(FDP_est_seq);
  };
  NumericVector getTime(){
    NumericVector  res(timer);   // 
    return(res);
  };
private:
  int rank_Mode;
  double fdr_level,init_p_cut;
  // start the timer
  Timer timer;
  
  int m,K;
   mat PVal,ProjPval,DistPval;
  //mat DistMat;
   mat HMat;
  sp_mat AdjMat_sp;
  vector<uvec> childInd;
  vector<uvec> parentInd;
  vec out_degree;
  rowvec in_degree;
  vec ProbInRej,ProbInRej_de,ProbInRej_en;
  
  double RejCount,ControlCount;
  double offset;
  double FDP_est;
  vec FDP_est_seq;
  
  vec RejSide;
  
  int rootNum=0;
  IntegerVector rootInd;
  
  int unmaskNum;
  uvec unmaskInd,excludeInd,maskInd;
  vec maskLogi;
  
  void createDAG();
  
  void InitData();
  void InitRejControl();
  void InitProbInRej(int startiter);
  void UpdateProbInRej(int oldrootNum,int removerootind, int RejInd);
  
  int UpdateRoot(int rm_rootind);
  
};

#endif







