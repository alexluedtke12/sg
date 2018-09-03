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


```
## Super Learner
```

```
## Version: 2.0-24
```

```
## Package created on 2018-08-10
```


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
sg.cvtmle(W,A,Y,txs=c(0,1),baseline.probs=c(1/2,1/2),SL.library=SL.library,sig.trunc=0.001,family=binomial(),kappa=1,num.SL.rep=2,num.est.rep=2,lib.ests=TRUE,verbose=FALSE)
```

```
## $est
## SuperLearner  SL.mean_All   SL.glm_All 
##   0.11957252   0.04102641   0.11957252 
## 
## $ci
##                      lb         ub
## SuperLearner 0.09299205 0.14615300
## SL.mean_All  0.01354829 0.06850453
## SL.glm_All   0.09299205 0.14615300
## 
## $est.mat
##              SuperLearner SL.mean_All SL.glm_All
## Repetition 1    0.1180312  0.04098487  0.1180312
## Repetition 2    0.1211139  0.04106795  0.1211139
```

```r
# truth
mean(Qbar(1,W.mc)*((blip.mc>0)-1/2) + Qbar(0,W.mc)*((blip.mc<=0)-1/2))
```

```
## [1] 0.1247774
```

```r
# Contrast optimal treatment strategy against treating everyone with tx 0
sg.cvtmle(W,A,Y,txs=c(0,1),baseline.probs=c(1,0),SL.library=SL.library,sig.trunc=0.001,family=binomial(),kappa=1,num.SL.rep=2,num.est.rep=2,lib.ests=TRUE,verbose=FALSE)
```

```
## $est
## SuperLearner  SL.mean_All   SL.glm_All 
##   0.15672328   0.08229836   0.15672328 
## 
## $ci
##                      lb        ub
## SuperLearner 0.11095379 0.2024928
## SL.mean_All  0.02744129 0.1371554
## SL.glm_All   0.11095379 0.2024928
## 
## $est.mat
##              SuperLearner SL.mean_All SL.glm_All
## Repetition 1    0.1542718  0.08180770  0.1542718
## Repetition 2    0.1591747  0.08278901  0.1591747
```

```r
# truth
mean((blip.mc)*(blip.mc>=0))
```

```
## [1] 0.16659
```

```r
# Resource constraint: at most 25% can be treated. Contrast against treating everyone with tx 0
sg.cvtmle(W,A,Y,txs=c(0,1),baseline.probs=c(1,0),SL.library=SL.library,sig.trunc=0.001,family=binomial(),kappa=0.25,num.SL.rep=2,num.est.rep=2,lib.ests=TRUE,verbose=FALSE)
```

```
## $est
## SuperLearner  SL.mean_All   SL.glm_All 
##  0.122349451 -0.001525606  0.121357394 
## 
## $ci
##                       lb         ub
## SuperLearner  0.09466465 0.15003425
## SL.mean_All  -0.02794915 0.02489793
## SL.glm_All    0.09362072 0.14909407
## 
## $est.mat
##              SuperLearner  SL.mean_All SL.glm_All
## Repetition 1    0.1232119 -0.006544996  0.1212211
## Repetition 2    0.1214870  0.003493783  0.1214936
```

```r
# truth
mean((blip.mc)*(blip.mc>=max(quantile(blip.mc,1-0.25),0)))
```

```
## [1] 0.1295534
```
