
#' SuperLearner for Estimating the Conditional Average Treatment Effect
#' 
#' @description
#' This function estimates the average treatment effect conditional on baseline covariates.
#' @details
#' If outcome is bounded in [0,1], then this functions respects that fact when estimating the outcome regression but not when estimating the conditional average treatment effect using the double robust loss presented in the below cited paper. Loss functions exist which do respect the bounds on the outcome (see the below cited paper), but we have not implemented them here.
#' @keywords SuperLearner, conditional average treatment effect, blip function
#' @param W data frame with observations in the rows and baseline covariates used to form the subgroup in columns.
#' @param A binary treatment vector, where 1 indicates that an individual was treated.
#' @param Y real-valued outcome.
#' @param SL.library SuperLearner library (see documentation for \code{SuperLearner} in the corresponding package) used to estimate outcome regression, treatment mechanism (if unknown), and conditional average treatment effect function.
#' @param txs A vector of length one or two indicating the two treatments of interest in A that will be used for the treatment assignment problem. Individuals with treatment \code{A} equal to a value not in \code{txs} will be treated as censored at treatment.
#' @param g0 if known (as in a randomized controlled trial), a vector of probabilities of treatment A=1 given covariates. If \code{NULL}, \code{SuperLearner} will be used to estimate these probabilities.
#' @param family \code{binomial()} if outcome bounded in [0,1], or \code{gaussian()} otherwise. See \code{Details}.
#' @param num.SL.folds number of folds to use in SuperLearner.
#' @param project if \code{TRUE} and family is \code{binomial()}, then projects the initial estimate of the blip function into its attainable range [-1,1].
#' @param num.SL.rep final output is an average of num.SL.rep super-learner fits (repetition ensures minimal reliance on initial choice of folds)
#' @param SL.method method that the SuperLearner function uses to select a convex combination of learners
#' @param id optional cluster identification variable
#' @param obsWeights observation weights
#' @param stratifyCV stratify validation folds by event counts (does this for estimation of outcome regression, treatment mechanism, and conditional average treatment effect function). Useful for rare outcomes
#' @param ipcw if TRUE, then does not estimate outcome regression (just sets it to zero for all covariate-treatment combinations)
#' @param lib.ests Also return the candidate optimal rule estimates in the super-learner library
#' @return a list containing
#' \item{est}{Vector containing an estimate of the conditional average treatment effect function for each individual in the data set (conditional on the covariate strata they belong to). If \code{txs} is of length two, then conditional average treatment effect is defined as the difference in conditional mean outcome if receiving the second treatment in \code{txs} versus receiving the first treatment in \code{txs}. If \code{txs} is of length one, then conditional average treatment effect is defined as the difference in conditional mean outcome if receiving the treatment in \code{txs} versus the expected outcome for a treatment randomly drawn according to the observed distribution (conditional on covariates).}
#' \item{SL}{\code{SuperLearner} object used to generate this estimate.}
#' @references A. R. Luedtke and M. J. van der Laan, ``Super-learning of an optimal dynamic treatment rule,'' \emph{International Journal of Biostatistics} (to appear), 2014.
#' @examples
#' SL.library = c('SL.mean','SL.glm')
#' Qbar = function(a,w){plogis(a*w$W1)}
#' n = 1000
#' W = data.frame(W1=rnorm(n),W2=rnorm(n),W3=rnorm(n),W4=rnorm(n))
#' A = rbinom(n,1,1/2)
#' Y = rbinom(n,1,Qbar(A,W))
#' 
#' cate.est = sg.SL(W,A,Y,SL.library=SL.library,family=binomial())$est
#' plot(W$W1,cate.est)
#' 
#' # compare to the truth
#' cate.truth = Qbar(1,W) - Qbar(0,W)
#' plot(cate.est,cate.truth)
#' @export


sg.SL = function(W,A,Y,SL.library,txs=c(0,1),g0=NULL,family=binomial(),num.SL.folds=10,project=TRUE,num.SL.rep=5,SL.method="method.NNLS2",id=NULL,obsWeights=NULL,stratifyCV=FALSE,ipcw=FALSE,lib.ests=FALSE,...){
	require(SuperLearner)

	if(any(names(list(...))=="separate.reg")) warning("The separate.reg option is deprecated. All outcome regressions are now stratified based on treatment status, i.e. the deprecated separate.reg always equals TRUE.")

	# Recode A
	A.new = rep(NA,length(A))
	if(length(txs)==2){
		# Recode A so that first value in txs is coded as 0, second is coded as 1
		A.new[A==txs[1]] = 0
		A.new[A==txs[2]] = 1
	} else if(length(txs)==1){
		# Recode A so that the value in txs is coded as 1
		A.new[A==txs[1]] = 1
	} else stop("txs must be of length 1 or 2")
	A.new[!(A%in%txs)] = 1e10
	A = A.new

	if((family$family=='binomial') & stratifyCV) {
		blip.cvControl = list(V=num.SL.folds,validRows=CVFolds(length(Y),id,Y,SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.SL.folds)))
	} else {
		if(sum(!(Y%in%c(0,1)))>0 & stratifyCV){
			warning('Can only stratify on binary outcome. Will stratify on estimation of treatment mechanism only.')
		}
		blip.cvControl = list(V=num.SL.folds)
	}
	n = nrow(W)
	if(ipcw){
		Qbar.1W <- Qbar.0W <- rep(0,nrow(W))
	} else {
		Qbar.1W = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(Y[A==1],W[A==1,],newX=W,family=family,SL.library=SL.library,cvControl=list(V=num.SL.folds,stratifyCV=(family$family=='binomial') & stratifyCV),id=id[A==1],obsWeights=obsWeights[A==1])$SL.predict[,1]}))/num.SL.rep
		if(length(txs)==2){
			Qbar.0W = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(Y[A==0],W[A==0,],newX=W,family=family,SL.library=SL.library,cvControl=list(V=num.SL.folds,stratifyCV=(family$family=='binomial') & stratifyCV),id=id[A==0],obsWeights=obsWeights[A==0])$SL.predict[,1]}))/num.SL.rep
		}
	}
	if(length(txs)==2) Qbar.AW = (A==1)*Qbar.1W + (A==0)*Qbar.0W
	if(length(g0)==0){
		g1 = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(as.numeric(A==1),W,family=binomial(),SL.library=SL.library,cvControl=list(V=num.SL.folds,stratifyCV=stratifyCV),id=id,obsWeights=obsWeights)$SL.predict[,1]}))/num.SL.rep
		if(length(unique(A))>2){
			g0 = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(as.numeric(A==0),W,family=binomial(),SL.library=SL.library,cvControl=list(V=num.SL.folds,stratifyCV=stratifyCV),id=id,obsWeights=obsWeights)$SL.predict[,1]}))/num.SL.rep
			if(sum(g0+g1>1)>0){
				inds = which(g0+g1>1)
				sm = g0 + g1
				g0[inds] = g0[inds]/sm[inds]
				g1[inds] = g1[inds]/sm[inds]
			}
		} else {
			g0 = 1-g1
		}
		g = cbind(g0,g1)
	} else g=g0
	if(length(txs)==2){
		Z = ((A==1)-(A==0))/((A==1)*g[,2] + (A==0)*g[,2] + (A!=0 & A!=1)) * (Y-Qbar.AW) + Qbar.1W - Qbar.0W
	} else {
		Z = (A==1)/((A==1)*g[,2] + (A!=1)) * (Y-Qbar.1W) + Qbar.1W - Y
	}
	

	SL.obj = SuperLearner(Z,W,SL.library=SL.library,cvControl=blip.cvControl,id=id,obsWeights=obsWeights,method=SL.method)
	blip = Reduce('+',lapply(1:num.SL.rep,function(i){SL.obj$SL.predict[,1]}))/num.SL.rep
	if(project & family$family=='binomial') blip = pmin(pmax(blip,-1),1)

	if(lib.ests){
		est.mat = Reduce('+',lapply(1:num.SL.rep,function(i){SL.obj$library.predict}))/num.SL.rep
		if(project & family$family=='binomial') est.mat = pmin(pmax(est.mat,-1),1)
		if(length(txs)==2) {
			return(list(est=blip,SL=SL.obj,Qbar.1W=Qbar.1W,Qbar.0W=Qbar.0W,Qbar.AW=Qbar.AW,Z=Z,est.mat=est.mat))
		} else {
			return(list(est=blip,SL=SL.obj,Qbar.1W=Qbar.1W,Z=Z,est.mat=est.mat))
		}
	} else {
		if(length(txs)==2) {
			return(list(est=blip,SL=SL.obj,Qbar.1W=Qbar.1W,Qbar.0W=Qbar.0W,Qbar.AW=Qbar.AW,Z=Z))
		} else {
			return(list(est=blip,SL=SL.obj,Qbar.1W=Qbar.1W,Z=Z))
		}
	}
}