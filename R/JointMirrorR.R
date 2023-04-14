#' Title
#' Joint Mirror Procedure
#'
#' A multiple testing procedure for testing simultaneous signals.
#' @param Pval a m*K matrix for p-values;
#' @param init.thred a value between 0 to 0.5 to determine the initial state;
#' @param offset a positive value at the enumerator of FDP estimator to guarantee finite sample FDR control;
#' @param trgt.fdr.level a value between 0 to 1 representing the target fdr level;
#' @param rank.Mode a character representing the selected partial order, including "Product","Max"/"Pmax","NoShape"/"EmptyPoset".
#' @param Ker.base a character indicating which part of the data to use to adjust the bandwidth, including "Pval","ProjPval","InProjPval","InPval".
#' @param Ker.BW the method for tuning bandwidth, including "Silverman", "Scott" and "Hpi".
#' @param is.gc a bool indicating whether we should free some spaces.
#'
#' @return list
#' selected is the indexes of the selected hypotheses.
#' Hker.inv.chol is the bandwidth matrix for kerne density estimation.
#' JMirrorS is a module storing the results of joint mirror procedure.
#'
#' @import kernelboot ks
#' @importFrom  methods new
#'
#' @export
#' @examples
#' set.seed(10)
#' m <- 1000
#' mu1H10 <- 1.5
#' mu1H11 <- 2
#' mu2H01 <- 2.5
#' mu2H11 <- 3
#' pi.seq <- c(0.4,0.2,0.2,0.2)
#' H <- rmultinom(m,1,pi.seq)
#' H0 <- colSums(H[1:3,])
#' mu1 <- mu2 <- rep(0,m)
#'
#' mu1[H[2,]==1] <- mu1H10
#' mu1[H[4,]==1] <- mu1H11
#' mu2[H[3,]==1] <- mu2H01
#' mu2[H[4,]==1] <- mu2H11
#'
#' X1 <- rnorm(m,mu1)
#' X2 <- rnorm(m,mu2)
#'
#' p1 <- 2*(pnorm(-abs(X1)))
#' p2 <- 2*(pnorm(-abs(X2)))
#'
#' trgt.fdr.level <- 0.2
#' Pval <- cbind(p1,p2)
#' JM.Product.Res <- JointMirror.R(Pval = cbind(p1,p2),rank.Mode = "Product",
#'                                 trgt.fdr.level=trgt.fdr.level)
#' JM.Max.Res <- JointMirror.R(Pval = cbind(p1,p2),rank.Mode = "Max",
#'                             trgt.fdr.level=trgt.fdr.level)
#' JM.EmptyPoset.Res <- JointMirror.R(Pval = cbind(p1,p2),rank.Mode = "EmptyPoset",
#'                                    trgt.fdr.level=trgt.fdr.level)
#'
#' c(fdp(JM.Product.Res$selected,H0),fdp(JM.Max.Res$selected,H0),fdp(JM.EmptyPoset.Res$selected,H0))
#' c(Pow(JM.Product.Res$selected,H0),Pow(JM.Max.Res$selected,H0),Pow(JM.EmptyPoset.Res$selected,H0))
JointMirror.R <- function(Pval,init.thred=0.2,
                          offset=1,trgt.fdr.level=0.1,
                          rank.Mode = "Product",
                          #h=1,#is.adapt.kernel =F,
                          Ker.base = "InProjPval",
                          Ker.BW = "Silverman", # Rule of Thumb
                          is.gc = T){
  ProjPval <- pmin(Pval,1-Pval)

  if(Ker.base=="Pval"){
    BinP <- Pval
  }else if(Ker.base=="ProjPval"){
    BinP <- ProjPval
  }else if(Ker.base=="InProjPval"){
    ind.Ex <- which(rowSums(Pval>=0.5)>1)
    BinP <- ProjPval[-ind.Ex,]
  }else if(Ker.base=="InPval"){
    ind.Ex <- which(rowSums(Pval>=0.5)>1)
    BinP <- Pval[-ind.Ex,]
  }

  if(Ker.BW == "Hpi"){
    # Wand & Jones (1995) and Chacon, J.E. & Duong, T. (2010)
    Ker.mat <- Hpi(BinP)
    Hker.inv.chol <- t(solve(chol(Ker.mat)))
    if(is.gc){
      gc()
    }
    #print(Hker.inv.chol)
  }else{
    n <- dim(BinP)[1]
    d <- dim(BinP)[2]
    if( Ker.BW == "Silverman"){
      #
      #SqrtKer.diag <-(4/(d+2))^(1/(d+4))*n^(-1/(d+4)) * apply(BinP,2,sd)
      # library(kernelboot)
       Ker.mat <- bw.silv(Pval, na.rm = FALSE)
    }
    if(Ker.BW =="Scott"){
      #SqrtKer.diag <-n^(-1/(d+4)) * apply(BinP,2,sd)
      Ker.mat <-  bw.scott(Pval,na.rm = FALSE)
    }
    #Hker.inv.chol <- diag(1/SqrtKer.diag)
    Hker.inv.chol <- t(solve(chol(Ker.mat)))
  }
  #DistMat <- as.matrix(dist(ProjPval %*% Hker.inv.chol))
  rank_Mode <- switch(rank.Mode,"Product"=1,
                      "Pmax"=2,"Max"=2,3)
  JMirrorSolver <- new(JointMirror,Pval,Hker.inv.chol)
  #cat("Run JM procedure...... \n")
  JMirrorSolver$InitPara(offset_=offset,fdr_level_=trgt.fdr.level,
                         init_p_cut_=init.thred,rank_Mode_=rank_Mode)
  JMirrorSolver$runJM()

  # Add 1 to make the index consistent
  selected <- JMirrorSolver$getRejInd()+1
  return(list(selected=selected,
              Hker.inv.chol = Hker.inv.chol,
              JMirrorS=JMirrorSolver))
}

#' Title
#' Closed FDP estimate
#'
#' If there are several target fdr levels to compare, we can set trgt.fdr.level=0 in JointMirror.R() and input result here. We can directly compare the FDP estimates and target fdr level to obtain discoveries.
#' @param Pval a m*K matrix for p-values;
#' @param JMirrorS is a module storing the results of joint mirror procedure; usually, this is the result of JMirrorR() with trgt.fdr.level=0.
#'
#' @return list
#' Unmask.Order is the unmasking order of each hypothesis.
#' is.RejSide represents whether the hypothesis is in the rejection side.
#' FDP.est.Total represents the adjusted FDP estimates that is decreasing with respect to the unmasking order.
#' @export
#'
#' @examples
#' set.seed(10)
#' m <- 1000
#' mu1H10 <- 1.5
#' mu1H11 <- 2
#' mu2H01 <- 2.5
#' mu2H11 <- 3
#' pi.seq <- c(0.4,0.2,0.2,0.2)
#' H <- rmultinom(m,1,pi.seq)
#' H0 <- colSums(H[1:3,])
#' mu1 <- mu2 <- rep(0,m)
#'
#' mu1[H[2,]==1] <- mu1H10
#' mu1[H[4,]==1] <- mu1H11
#' mu2[H[3,]==1] <- mu2H01
#' mu2[H[4,]==1] <- mu2H11
#'
#' X1 <- rnorm(m,mu1)
#' X2 <- rnorm(m,mu2)
#'
#' p1 <- 2*(pnorm(-abs(X1)))
#' p2 <- 2*(pnorm(-abs(X2)))
#'
#' Pval = cbind(p1,p2)
#' JM.Product.Res <- JointMirror.R(Pval,rank.Mode = "Product",trgt.fdr.level=0)
#' JM.Product.Qval <- JointMirror.Qvalue(Pval,JM.Product.Res$JMirrorS)
#' trgt.fdr.level <- 0.2
#' FDP.est <- JM.Product.Qval$FDP.est.Total
#' is.RejSide <- JM.Product.Qval$is.RejSide
#' JM.Product.sel <- which(FDP.est<=trgt.fdr.level & is.RejSide)
#' c(fdp(JM.Product.sel,H0),Pow(JM.Product.sel,H0))
JointMirror.Qvalue <- function(Pval,JMirrorS){
  ## Provide ``Qvalue"
  m <- dim(Pval)[1]
  is.RejSide <- apply(Pval,1,function(x){min(x<0.5)})

  FDP.est.val <- JMirrorS$getFDPest()
  # In.Ind is the indexes of unmasked hypothesis
  # In.Ind also includes the order of unmasking
  In.Ind <- JMirrorS$getunMaskInd()+1
  In.Ind <- In.Ind[In.Ind!=m+1]
  Ex.Ind <- (1:m)[!(1:m %in% In.Ind)]

  # To ensure FDP decreases as rejection region shrinks
  FDP.est.In <- FDP.est.val[In.Ind]
  FDP.est.In <- cummin(FDP.est.In)
  # FDP.est.In <- sapply(1:length(In.Ind),
  #                      function(i){
  #                        min(FDP.est.val[In.Ind[1:i]])})

  FDP.est.Total <- rep(1,m)
  FDP.est.Total[In.Ind] <-  FDP.est.In
  FDP.est.Total <- pmin(1,FDP.est.Total)
  return(list(Unmask.Order = rev(In.Ind),
              is.RejSide=is.RejSide,
              FDP.est.Total=FDP.est.Total))
}
