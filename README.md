---
title: "sg"
output: rmarkdown::github_document
---



# Description

This package implements methods for methods for conducting subgroup analyses and estimating the impact of implementing the optimal individualized treatment strategy in the population.

# Installation

To install the package via [devtools](https://www.rstudio.com/products/rpackages/devtools/), use


```r
devtools::install_github("alexluedtke12/sg")
```

# Example

We now present an example on a sample of size n=1000. The SuperLearner library is small so that the example does not take long to run -- in practice, we recommend using a larger SuperLearner library. 




```r
library(sg)

sim.data = function(nn){
	W = data.frame(W1=rnorm(nn),W2=rnorm(nn),W3=rnorm(nn),W4=rnorm(nn))
	A = rbinom(nn,1,1/2)
	Y = rbinom(nn,1,Qbar(A,W))
	return(list(W=W,A=A,Y=Y))
}

SL.library = c('SL.mean','SL.glm')

n = 1000
Qbar = function(a,w){plogis(-1 + 2*a*w$W1)} # function(a,w){plogis(a*(w$W1>0)*w$W1)}

# data to run methods on
dat = sim.data(n)
W = dat$W
A = dat$A
Y = dat$Y

# data to evaluate true parameter values
n.mc = 2e4
dat = sim.data(n.mc)
W.mc = dat$W
A.mc = dat$A
Y.mc = dat$Y
blip.mc = Qbar(1,W.mc)-Qbar(0,W.mc)
```

Below each example, we also print the true impact of implementing the rule, i.e. the quantity that sg.cvtmle is estimating.

```r
# Contrast optimal treatment strategy against treating with probability 1/2
sg.cvtmle(W,A,Y,SL.library=SL.library,txs=c(0,1),baseline.probs=c(1/2,1/2),CATE.SL.library=SL.library,sig.trunc=0.001,family=binomial(),kappa=1,num.SL.rep=2,num.est.rep=2,lib.ests=TRUE,verbose=FALSE)
```

```
## $est
## SuperLearner      SL.mean       SL.glm 
##   0.13394596   0.06408995   0.13412034 
## 
## $ci
##                      lb         ub
## SuperLearner 0.10831312 0.15957880
## SL.mean      0.03743888 0.09074101
## SL.glm       0.10849102 0.15974966
## 
## $est.mat
##              SuperLearner    SL.mean    SL.glm
## Repetition 1    0.1357762 0.06413537 0.1357762
## Repetition 2    0.1321157 0.06404452 0.1324645
```

```r
# truth
mean(Qbar(1,W.mc)*((blip.mc>0)-1/2) + Qbar(0,W.mc)*((blip.mc<=0)-1/2))
```

```
## [1] 0.1242482
```

```r
# Contrast optimal treatment strategy against treating everyone with tx 0
sg.cvtmle(W,A,Y,SL.library=SL.library,txs=c(0,1),baseline.probs=c(1,0),CATE.SL.library=SL.library,sig.trunc=0.001,family=binomial(),kappa=1,num.SL.rep=2,num.est.rep=2,lib.ests=TRUE,verbose=FALSE)
```

```
## $est
## SuperLearner      SL.mean       SL.glm 
##    0.1991847    0.1279397    0.1990086 
## 
## $ci
##                      lb        ub
## SuperLearner 0.15403755 0.2443319
## SL.mean      0.07459551 0.1812839
## SL.glm       0.15386105 0.2441562
## 
## $est.mat
##              SuperLearner   SL.mean    SL.glm
## Repetition 1    0.1994651 0.1273042 0.1991129
## Repetition 2    0.1989044 0.1285752 0.1989044
```

```r
# truth
mean((blip.mc)*(blip.mc>=0))
```

```
## [1] 0.1655005
```

```r
# Resource constraint: at most 25% can be treated. Contrast against treating everyone with tx 0
sg.cvtmle(W,A,Y,SL.library=SL.library,txs=c(0,1),baseline.probs=c(1,0),CATE.SL.library=SL.library,sig.trunc=0.001,family=binomial(),kappa=0.25,num.SL.rep=2,num.est.rep=2,lib.ests=TRUE,verbose=FALSE)
```

```
## $est
## SuperLearner      SL.mean       SL.glm 
##    0.1863364    0.1280433    0.1863364 
## 
## $ci
##                      lb        ub
## SuperLearner 0.13782888 0.2348439
## SL.mean      0.07467649 0.1814101
## SL.glm       0.13782509 0.2348477
## 
## $est.mat
##              SuperLearner   SL.mean    SL.glm
## Repetition 1    0.1874353 0.1293991 0.1874353
## Repetition 2    0.1852375 0.1266875 0.1852375
```

```r
# truth
mean((blip.mc)*(blip.mc>=max(quantile(blip.mc,1-0.25),0)))
```

```
## [1] 0.1290834
```
