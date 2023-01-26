#include "JointMirror.hpp"

// [[Rcpp::depends(RcppArmadillo)]]

JointMirror::JointMirror( mat Pval_, mat HMat_){
  PVal=Pval_;
  ProjPval=min(PVal,1-PVal);
  HMat=HMat_;
  
  DistPval=ProjPval*HMat;
  
  K = PVal.n_cols;
  m = PVal.n_rows;
  
}

void JointMirror::InitPara(double offset_, double fdr_level_,
                           double init_p_cut_,
                           int rank_Mode_){
  fdr_level=fdr_level_;
  offset=offset_;
  rank_Mode=rank_Mode_;
  init_p_cut=init_p_cut_;
}

void JointMirror::InitData(){
  //Rcpp::Rcout << "InitData" <<  std::endl;
  
  //Rcpp::Rcout << "11"<<  std::endl;
  ProbInRej = zeros(m);
  ProbInRej_de = zeros(m);
  ProbInRej_en = zeros(m);
  
  out_degree = zeros(m);
  in_degree = zeros<rowvec>(m);
  unmaskInd = zeros<uvec>(m);
  unmaskInd.fill(m);
  rootNum=0;
  
  
  //Rcpp::Rcout << "22" <<  std::endl;
  childInd.resize(m);
  parentInd.resize(m);
  
  RejSide = zeros(m);
  for(int i=0;i<m;i++){
    RejSide(i) = all(PVal.row(i)<0.5);
  }
  RejCount=0;
  ControlCount=0;
    
  // Initialize exclude and unmask index
  vec excludeLogi = zeros(m);
  vec unmaskLogi = zeros(m);
  for(int i=0;i<m;i++){
    excludeLogi(i) = sum(PVal.row(i)>=0.5)>1;
    unmaskLogi(i) = max(0.0,
               (max(ProjPval.row(i))>=init_p_cut)-excludeLogi(i));
  }
//  unmaskLogi <- max(zeros(m),unmaskLogi-excludeLogi);
  
  unmaskNum = sum(unmaskLogi);
  unmaskInd(regspace<uvec>(0,unmaskNum-1)) = find(unmaskLogi==1);
  
  excludeInd = find(excludeLogi==1);
  
  maskInd = find(excludeLogi==0&&unmaskLogi==0);
  
  maskLogi.ones(m);
  maskLogi(find(excludeLogi==1||unmaskLogi==1)).zeros();
  
  
  // Init Transitive Reductive Matrix
  // Generate Closure firstly
  
  timer.step("GenRedMat");
  // (x1,...,xk)<(y1,...,yk) iff x1<y1,...,xk<yk
  
  timer.step("Create Transitive Closure");
  int i,j;
  if(rank_Mode!=3){
    int maskNum  = maskInd.n_elem;
    // We just need to create AdjMat for the masked projected p-values to save space
    mat AdjMat=mat(maskNum, maskNum, fill::zeros);
  if(rank_Mode==1){
    //for(auto &i:maskInd){
    //for(auto &j:maskInd){
    for(int isub=0;isub<maskNum;isub++){
      for(int jsub=0;jsub<maskNum;jsub++){
        // Avoid cycle
        i = maskInd(isub);
        j = maskInd(jsub);
        if(all(ProjPval.row(i)>=ProjPval.row(j))&&any(ProjPval.row(i)!=ProjPval.row(j))){
          AdjMat(isub,jsub)=1;
        }
      }
    }
  }
  
  // (x1,...,xk)<(y1,...,yk) iff max(x1,...,xk)<max(y1,...,yk)
  if(rank_Mode==2){
    for(int isub=0;isub<maskNum;isub++){
      for(int jsub=0;jsub<maskNum;jsub++){
        // Avoid cycle
        i = maskInd(isub);
        j = maskInd(jsub);
        if(max(ProjPval.row(i))>max(ProjPval.row(j))){
          AdjMat(isub,jsub)=1;
        }
      }
    }
  }
  
  // Otherwise (x1,...,xk) are not comparable; Turn to rankmode 3
  // Do not suggest
  // Re_Initilize
  out_degree=vec(maskNum, fill::zeros);
  in_degree=rowvec(maskNum, fill::zeros);
  
  out_degree=sum(AdjMat,1);
  in_degree=sum(AdjMat,0);
  
  timer.step("Create Transitive Reuction");
  //mat AdjMatSub = AdjMat()
  //AdjMat = AdjMat - ((AdjMat*AdjMat)!=0);
  
  // Try to accelerate
  uvec out_sort_ind = sort_index(out_degree,"descend");
  uvec out_sort_ind_sort = sort_index(out_sort_ind,"ascend");
  
  AdjMat = AdjMat.cols(out_sort_ind);
  
  vec in_degree_col=conv_to< vec >::from(in_degree);
  in_degree_col = in_degree_col(out_sort_ind);
  //out_degree = out_degree(out_sort_ind);
  
  uvec withTransparentInd = find(in_degree_col>0);
  uvec out_sort_ind_withTransparent = out_sort_ind(withTransparentInd);
  mat AdjMatTrunc = AdjMat.cols(withTransparentInd);
  mat AdjMatCopy=AdjMatTrunc;
  // Only vertex has at least two children will have a transition path go through another vertex
  uvec withTranschildInd = find(out_degree>1);
  uvec iTransparentInd;
  
  timer.step("Create Transitive Reduction(Reordering)");
  //uvec withcpInd = find(out_degree!=0&&in_degree!=0);
  // Our establishment will be more efficient after reordering
  for(auto &i:withTranschildInd){
    // iTransparentInd = find(AdjMatTrunc.row(i)==1);
    // for(auto &j:iTransparentInd){
    //   for(auto &k:iTransparentInd){
    //     // i->k (must stand) and k->j(only thing we need to check) 
    //     // With stored matrix, we can jump from the loop quickly
    //     if(AdjMatCopy(k,j)==1){
    //       AdjMatTrunc(i,j)=0;
    //       break;
    //     }
    //   }
    iTransparentInd = find(AdjMatTrunc.row(i)==1);
    //for(auto j:iTransparentInd){
    for(auto j=iTransparentInd.begin()+1;j<iTransparentInd.end();j++){
      // If there is a path i->k->j, outdegree(k)>outdegree(j)
    for(auto k=iTransparentInd.begin();k<j;k++){
      //for(auto k:iTransparentInd){
        // There is a path i->k->j ( so that i->j comes from transition)
        if(AdjMatCopy(out_sort_ind_withTransparent(*k),(*j))==1){
          //cout<<i<<" "<<k<<" "<<j<<endl;
          AdjMatTrunc(i,(*j))=0;
          break;
        }
      }
    }
  }
  
  AdjMat.cols(withTransparentInd) = AdjMatTrunc;
  AdjMat = AdjMat.cols(out_sort_ind_sort); // Recover the previous index
  
  timer.step("Create Transitive Reduction Degree");
  // Recalculate outdegree and in degree
  out_degree=vec(m, fill::zeros);
  in_degree=rowvec(m, fill::zeros);
  
  out_degree(maskInd)=sum(AdjMat,1);
  in_degree(maskInd)=sum(AdjMat,0);
  
  AdjMat_sp = sp_mat(AdjMat); // The row and column name is maskInd
  }else{
  // If rank mode is 3, adjmat is zero, i.e., no connection;
    out_degree=vec(m, fill::zeros);
    in_degree=rowvec(m, fill::zeros);
    //AdjMat_sp = sp_mat(m,m);
  }
  timer.step("GenDAG");
  //Rcpp::Rcout << "CreateDAG" <<  std::endl;
  createDAG();
  
  timer.step("InitRej");
  //Rcpp::Rcout << "InitProbInRej" <<  std::endl;
  InitProbInRej(0);
  //Rcpp::Rcout << "InitRejControl" <<  std::endl;
  InitRejControl();
}


void JointMirror::InitRejControl(){
  for(auto&iterI:maskInd){
    RejCount += RejSide(iterI);
    ControlCount += (1-RejSide(iterI));
  }
  FDP_est = (ControlCount+offset)/max(RejCount,1.0);
  
  FDP_est_seq = ones(m);//*FDP_est;
}

void JointMirror::createDAG(){
  uvec withchildInd = find(out_degree!=0);
  uvec withparentInd = find(in_degree!=0);
  // Initlize according to the out degree and in degree
  timer.step("Create ChildSet");
  for(auto &i:withchildInd){
    childInd.at(i) = uvec(out_degree(i));
  }
  
  for(auto &i:withparentInd){
    parentInd.at(i) = uvec(in_degree(i));
  }
  
  // Find the root index for DAG
  timer.step("Create RootIndex");
  vec in_degree_col=conv_to< vec >::from(in_degree);
  rootInd = find((maskLogi==1)&&(in_degree_col==0));
  rootNum = sum((maskLogi==1)&&(in_degree_col==0));
  //Rcout<<sum(in_degree_col)<<endl;
  //Create DAG
  uvec child_index(m,fill::zeros);
  uvec parent_index(m,fill::zeros);
  
  timer.step("Create Sparse copy");
  sp_mat::const_iterator it     = AdjMat_sp.begin();
  sp_mat::const_iterator it_end = AdjMat_sp.end();
  uword row_cur_mask,col_cur_mask;
  int row_cur,col_cur;
  timer.step("Start Constructing DAG");
  for(; it != it_end; ++it){
    // Since our matrix is a 0-1 matrix & sparse;
    // It has value means there is a edge;
    row_cur_mask=it.row();
    col_cur_mask=it.col();
    row_cur=maskInd(row_cur_mask);
    col_cur=maskInd(col_cur_mask);
    //cout << "ind: " << row_cur<<" "<<col_cur<<", val:"<<  (*it)    << endl;
    childInd.at(row_cur)(child_index(row_cur)) = col_cur;
    parentInd.at(col_cur)(parent_index(col_cur)) = row_cur;
    child_index(row_cur)=child_index(row_cur)+1;
    parent_index(col_cur)=parent_index(col_cur)+1;
  }
  
  // for(auto &i:maskInd){
  //   //Rcpp::Rcout << i <<  std::endl;
  //   //Rcpp::Rcout << find(AdjMat.row(i)!=0) <<  std::endl;
  //   if(out_degree(i)!=0){
  //     childInd.at(i)=find(AdjMat_sp.row(i)!=0);
  //   }
  //   //Rcpp::Rcout << find(AdjMat.col(i)!=0) <<  std::endl;
  //   
  //   if(in_degree(i)!=0){
  //   parentInd.at(i)=find(AdjMat_sp.col(i)!=0);
  //   }
  //   if(in_degree(i)==0){
  //     rootInd.push_back(i);
  //     rootNum = rootNum+1;
  //   }
  // }
}

//Detele a root from root set
//We may include new roots
int JointMirror::UpdateRoot(int rm_rootind){
  //Rcpp::Rcout << regspace<uvec>(0,unmaskNum-1) <<  std::endl;
  //Rcpp::Rcout << rootInd <<  std::endl;
  int newNum=0;
  int root_unmask_all = rootInd(rm_rootind);
  //Rcpp::Rcout << "Before" <<rootInd<<  std::endl;
  //Rcpp::Rcout << "Remove" <<rm_rootind<<  std::endl;
  rootInd.erase(rm_rootind); // Remove from the root set
  // Unmask
  //Rcpp::Rcout << "Unmask" <<root_unmask_all<<  std::endl;
  maskLogi(root_unmask_all) = 0; 
  unmaskInd(unmaskNum) =root_unmask_all;
  unmaskNum = unmaskNum+1;
  
  rootNum-=1;
  //Rcpp::Rcout << "After" <<rootInd<<  std::endl;
  
  // Add new roots
  uvec childind = childInd.at(root_unmask_all);
  //Rcpp::Rcout << "Add new" <<  childind<<std::endl;
  if(!childind.is_empty()){
    for(auto &iterI:childind){
      in_degree(iterI) -= 1;
      if(in_degree(iterI)==0){
        rootNum+=1;
        newNum+=1;
        rootInd.push_back(iterI);
      }
    }
  }
  return(newNum);
}



//Waiting:ChangeTheRangeOfIter
void JointMirror::InitProbInRej(int startiter){
  //Rcpp::Rcout << regspace<uvec>(0,unmaskNum-1) <<  std::endl;
  //Rcpp::Rcout << rootInd <<  std::endl;
  uvec unmaskind = unmaskInd(regspace<uvec>(0,unmaskNum-1));
  int rootind;
  vec Distvec;
  vec Kervec;
  rowvec rootDistPval;
  mat DiffDistPval;
  for(int iterI=startiter;iterI<rootNum;iterI++){
    //Rcpp::Rcout << "iter "<<iterI<<" rootNum"<<rootNum<<"rootInd "<<rootInd <<  std::endl;
    rootind = rootInd.at(iterI);
    rootDistPval = DistPval.row(rootind);
    //Rcpp::Rcout << "DistMat "<<iterI<<" rootind" <<  std::endl;
    //X.cols(0,3).each_col() -= v;
    //XX.each_col( [](const vec& b){ b.print(); } )
    //Distvec = DistPval.rows(unmaskind).each_row( [](const vec& b){ norm(rootDistPval-b,2); } ); //DistMat.col(rootind);
    DiffDistPval = DistPval.rows(unmaskind);
    DiffDistPval.each_row()-=rootDistPval;
    Distvec = sqrt(sum(DiffDistPval%DiffDistPval,1)); //DistMat.col(rootind);
    //Distvec = sqrt(DistPval.row(rootind) * DistPval.row(unmaskind).t()); //DistMat.col(rootind);
    Kervec = normpdf(Distvec);//normpdf(Distvec(unmaskind));
    //Rcpp::Rcout << normpdf(Distvec(unmaskind)) <<  std::endl;
    ProbInRej_de(rootind) = sum(Kervec);
    //Rcpp::Rcout << "RejSide"<<RejSide(unmaskind) <<  std::endl;
    ProbInRej_en(rootind) = dot(Kervec,RejSide(unmaskind));
    ProbInRej(rootind) = ProbInRej_en(rootind)/ProbInRej_de(rootind);
  }
}


void JointMirror::UpdateProbInRej(int oldrootNum,int removerootind, int RejInd){
  //vec Distvec = DistMat.col(removerootind);
  int orootind;
  double w_val;
  
  for(int oldrootind=0;oldrootind<oldrootNum;oldrootind++){
    orootind = rootInd.at(oldrootind);
    w_val = normpdf(norm(DistPval.row(orootind)-DistPval.row(removerootind),2));
    ProbInRej_de(orootind) += w_val;//normpdf(Distvec(orootind));
    ProbInRej_en(orootind) += w_val*RejInd;//normpdf(Distvec(orootind))*RejInd;
    ProbInRej(orootind) = ProbInRej_en(orootind)/ProbInRej_de(orootind);
  }
}

void JointMirror::runJM(){
  int root_unmask_in,root_unmask_all;
  int is_Rej;
  int newrootNum;
  int SearchInd=0;
  timer.step("InitData");
  InitData();
  //Rcpp::Rcout << "StartRun" <<  std::endl;
  
  timer.step("StartRun");
  while(FDP_est>=fdr_level && RejCount>0){
    //Rcpp::Rcout << "Start Searching" << SearchInd<< std::endl;
    //Rcpp::Rcout << "rootInd: " << rootInd<< " "<<as<arma::uvec>(rootInd) << std::endl;
    root_unmask_in = ProbInRej(as<arma::uvec>(rootInd)).index_min();
    root_unmask_all = rootInd(root_unmask_in);
    
    //Rcpp::Rcout <<"MinProj: "<<ProbInRej(as<arma::uvec>(rootInd)).min()<<" Pval: "<<PVal.row(root_unmask_all)<<" RejSide: "<<RejSide(root_unmask_all)<< std::endl;
    
    RejCount -= RejSide(root_unmask_all);
    ControlCount -= (1-RejSide(root_unmask_all));
    
    FDP_est = (ControlCount+offset)/max(RejCount,1.0);
    FDP_est_seq(root_unmask_all) = FDP_est;
    //Rcpp::Rcout<<"Root: "<<root_unmask_all << " RejCount: "<<RejCount<< " ControlCount: "<<ControlCount<< " FDP_est:"<<FDP_est<<  std::endl;
    
    //Rcpp::Rcout << "Update Root Set" <<  std::endl;
    // Update Root Set
    newrootNum = UpdateRoot(root_unmask_in);
    //Rcpp::Rcout << "Update InitProbInRej Set" <<  std::endl;
    // Update InitProbInRej Set
    if(newrootNum!=0){
      InitProbInRej(rootNum-newrootNum);
    }
    if(rootNum!=newrootNum){
      UpdateProbInRej(rootNum-newrootNum,root_unmask_all,RejSide(root_unmask_all));
    }
    SearchInd++;
    
    //Rcpp::Rcout <<"FDP: "<<FDP_est<<" fdr_level:"<<fdr_level<< " UnmaskIn: " << root_unmask_all<<" RootNum:"<<rootNum<<" unmaskNum:"<<unmaskNum<<" NewRootNum:"<<newrootNum<< std::endl;
    
  }
  
  timer.step("EndRun");
  if(RejCount!=0){
    FDP_est_seq(find(maskLogi==1)).fill(FDP_est);
  }
}

void JointMirror::print(){
  Rcpp::Rcout << "We reject "<<RejCount<<" hypothesises with the FDP estimates "<< 
    FDP_est<<" with offset "<< offset<< ". There are "<<
    ControlCount<<" hypothesis locating into the control region."<<  std::endl;
}

