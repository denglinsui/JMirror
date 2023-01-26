# JMirror

This R package is for large-scale inference for testing signals in several experiments simultaneously. 

# Installation

```
library(devtools)
devtools::install_github("https://github.com/denglinsui/JMirror")
```

# Example 
```
library(JMirror)

set.seed(10)
m <- 1000
mu1H10 <- 1.5
mu1H11 <- 2
mu2H01 <- 2.5
mu2H11 <- 3
pi.seq <- c(0.4,0.2,0.2,0.2)
H <- rmultinom(m,1,pi.seq)
H0 <- colSums(H[1:3,])
mu1 <- mu2 <- rep(0,m)

mu1[H[2,]==1] <- mu1H10
mu1[H[4,]==1] <- mu1H11
mu2[H[3,]==1] <- mu2H01
mu2[H[4,]==1] <- mu2H11

X1 <- rnorm(m,mu1)
X2 <- rnorm(m,mu2)

p1 <- 2*(pnorm(-abs(X1)))
p2 <- 2*(pnorm(-abs(X2)))

Pval <- cbind(p1,p2)
trgt.fdr.level <- 0.2

# Implementation suggested for single target fdr level
JM.Product.Res <- JointMirror.R(Pval,rank.Mode = "Product",
                                trgt.fdr.level=trgt.fdr.level)
JM.Max.Res <- JointMirror.R(Pval),rank.Mode = "Max",
                            trgt.fdr.level=trgt.fdr.level)
JM.EmptyPoset.Res <- JointMirror.R(Pval,rank.Mode = "EmptyPoset",
                                   trgt.fdr.level=trgt.fdr.level)

c(fdp(JM.Product.Res$selected,H0),fdp(JM.Max.Res$selected,H0),fdp(JM.EmptyPoset.Res$selected,H0))
c(Pow(JM.Product.Res$selected,H0),Pow(JM.Max.Res$selected,H0),Pow(JM.EmptyPoset.Res$selected,H0))

# Implementation suggested for several target fdr levels
JM.Product.Res <- JointMirror.R(Pval,rank.Mode = "Product",trgt.fdr.level=0)
JM.Product.Qval <- JointMirror.Qvalue(Pval,JM.Product.Res$JMirrorS)
trgt.fdr.level <- 0.2
FDP.est <- JM.Product.Qval$FDP.est.Total
is.RejSide <- JM.Product.Qval$is.RejSide
JM.Product.sel <- which(FDP.est<=trgt.fdr.level & is.RejSide==T)
c(fdp(JM.Product.sel,H0),Pow(JM.Product.sel,H0))
```
