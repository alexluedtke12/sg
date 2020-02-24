sg
================

Description
===========

This package implements methods for methods for conducting subgroup analyses and estimating the impact of implementing the optimal individualized treatment strategy in the population.

Installation
============

To install the package via [devtools](https://www.rstudio.com/products/rpackages/devtools/), use

``` r
devtools::install_github("alexluedtke12/sg")
```

Example
=======

We now present an example on a sample of size n=1000. The SuperLearner library is small so that the example does not take long to run -- in practice, we recommend using a larger SuperLearner library.

    ## Warning: package 'SuperLearner' was built under R version 3.5.2

    ## Super Learner

    ## Version: 2.0-25

    ## Package created on 2019-08-05

``` r
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

``` r
# Contrast optimal treatment strategy against treating with probability 1/2
sg.cvtmle(W,A,Y,txs=c(0,1),baseline.probs=c(1/2,1/2),CATE.SL.library=SL.library,sig.trunc=0.001,family=binomial(),kappa=1,num.SL.rep=2,num.est.rep=2,lib.ests=TRUE,verbose=FALSE)
```

    ## $est
    ## SuperLearner      SL.mean       SL.glm 
    ##    0.1002347    0.0481884    0.1002347 
    ## 
    ## $ci
    ##                      lb         ub
    ## SuperLearner 0.07432583 0.12614355
    ## SL.mean      0.02170676 0.07467005
    ## SL.glm       0.07432583 0.12614355
    ## 
    ## $est.mat
    ##              SuperLearner    SL.mean     SL.glm
    ## Repetition 1   0.09570899 0.04835484 0.09570899
    ## Repetition 2   0.10476039 0.04802197 0.10476039

``` r
# truth
mean(Qbar(1,W.mc)*((blip.mc>0)-1/2) + Qbar(0,W.mc)*((blip.mc<=0)-1/2))
```

    ## [1] 0.1251412

``` r
# Contrast optimal treatment strategy against treating everyone with tx 0
sg.cvtmle(W,A,Y,txs=c(0,1),baseline.probs=c(1,0),CATE.SL.library=SL.library,sig.trunc=0.001,family=binomial(),kappa=1,num.SL.rep=2,num.est.rep=2,lib.ests=TRUE,verbose=FALSE)
```

    ## $est
    ## SuperLearner      SL.mean       SL.glm 
    ##   0.15120948   0.09685982   0.15040278 
    ## 
    ## $ci
    ##                      lb        ub
    ## SuperLearner 0.10614544 0.1962735
    ## SL.mean      0.04391033 0.1498093
    ## SL.glm       0.10538518 0.1954204
    ## 
    ## $est.mat
    ##              SuperLearner    SL.mean    SL.glm
    ## Repetition 1    0.1531282 0.09701104 0.1515148
    ## Repetition 2    0.1492907 0.09670860 0.1492907

``` r
# truth
mean((blip.mc)*(blip.mc>=0))
```

    ## [1] 0.1674317

``` r
# Resource constraint: at most 25% can be treated. Contrast against treating everyone with tx 0
sg.cvtmle(W,A,Y,txs=c(0,1),baseline.probs=c(1,0),CATE.SL.library=SL.library,sig.trunc=0.001,family=binomial(),kappa=0.25,num.SL.rep=2,num.est.rep=2,lib.ests=TRUE,verbose=FALSE)
```

    ## $est
    ## SuperLearner      SL.mean       SL.glm 
    ##   0.14105268   0.09691537   0.14105268 
    ## 
    ## $ci
    ##                      lb        ub
    ## SuperLearner 0.09308925 0.1890161
    ## SL.mean      0.04399252 0.1498382
    ## SL.glm       0.09308832 0.1890170
    ## 
    ## $est.mat
    ##              SuperLearner    SL.mean    SL.glm
    ## Repetition 1    0.1389333 0.09686683 0.1389333
    ## Repetition 2    0.1431721 0.09696391 0.1431721

``` r
# truth
mean((blip.mc)*(blip.mc>=max(quantile(blip.mc,1-0.25),0)))
```

    ## [1] 0.1297942
