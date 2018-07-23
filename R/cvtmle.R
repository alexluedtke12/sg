# Standard deviation truncated at sig.trunc
# g.est is vector giving estimated probability of A=1 given W

.sgtmle = function(W,A,Y,Q.est,blip.est,g.est,family=binomial(),kappa=1,sig.trunc=0.1,alpha=0.05,trunc.Q=c(0.005,0.995),obsWeights=NULL,RR=FALSE){
	n = nrow(W)
	if(is.null(obsWeights)) obsWeights = rep(1,length(Y))	
	if(family$family=='binomial') for(i in 1:2){Q.est[,i] = pmin(pmax(Q.est[,i],trunc.Q[1]),trunc.Q[2])}

	offs = family$linkfun(A*Q.est[,2] + (1-A)*Q.est[,1])

	tau = max(quantile(blip.est,1-kappa,type=1),0)

	if(mean(blip.est==tau)==0) {
		g.star = as.numeric(blip.est>tau)
	} else g.star = (blip.est==tau)*(tau>0)*(kappa - mean(blip.est>tau))/mean(blip.est==tau) + (blip.est>tau)

	if(sum(g.star)>0){
		H = (2*A-1)/(A*g.est + (1-A)*(1-g.est)) * g.star
		eps = suppressWarnings(coef(glm(Y ~ -1 + offset(offs) + H,family=family,weights=obsWeights)))
		Q0 = family$linkinv(family$linkfun(Q.est[,1]) - eps*g.star/(1-g.est))
		Q1 = family$linkinv(family$linkfun(Q.est[,2]) + eps*g.star/g.est)

		est.add = mean((Q1-Q0)*g.star)
		ic.add = obsWeights * (H*(Y-(A*Q1 + (1-A)*Q0)) - tau*(g.star-kappa) + g.star*(Q1-Q0) - est.add)
		if(RR){ # Get a targeted estimate Q0.star of EY0, and then add this to the contrast from the previous
			offs.0 = family$linkfun(Q.est[,1])
			H.0 = (1-A)/(1-g.est)
			eps.0 = suppressWarnings(coef(glm(Y ~ -1 + offset(offs.0) + H.0,family=family,weights=obsWeights)))
			Q0.star = family$linkinv(offs.0 + eps.0/(1-g.est))
			EY0.star = mean(Q0.star)
			est = (1 - (EY0.star + est.add))/(1-EY0.star) # recall: relative risk of 1-Y, since Y beneficial
			ic = -(ic.add/(1 - (EY0.star + est.add)) + obsWeights * (H.0 * (Y-Q0.star) + Q0.star - EY0.star) * (1/(1 - (EY0.star + est.add))-1/(1-EY0.star)))
			ci = exp(log(est) + c(-1,1) * qnorm(1-alpha/2)*max(sd(ic),sig.trunc)/sqrt(n))
			print(mean(ic))
		} else {
			est = est.add; ic = ic.add
			ci = c(est-qnorm(1-alpha/2)*max(sd(ic),sig.trunc)/sqrt(n),est+qnorm(1-alpha/2)*max(sd(ic),sig.trunc)/sqrt(n))
		}
	} else {
		est = ifelse(RR,1,0)
		if(RR){
			ci = exp(c(-1,1) * qnorm(1-alpha/2)*sig.trunc/sqrt(n))
		} else {
			ci = c(-1,1) * qnorm(1-alpha/2)*sig.trunc/sqrt(n)
		}
		ic = rep(0,length(Y))
	}
	return(list(est=est,ci=ci,ic=ic))
}



# Helper function used by sg.cvtmle.bothactive and sg.cvtmle
# Inputs the same as for those functions

.sgcvtmle.preprocess = function(W,A,Y,SL.library,g0=NULL,family=binomial(),num.folds=10,num.SL.rep=5,id=NULL,obsWeights=NULL,stratifyCV=FALSE,ipcw=FALSE){
	require(SuperLearner)
	n = nrow(W)
	if(!is.null(id) & stratifyCV==TRUE) stop('Stratified sampling with id not currently implemented.')
	folds = CVFolds(n,id,Y,SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.folds))
	ests = do.call(rbind,lapply(folds,function(val.inds){

		W.trn = W[-val.inds,]
		A.trn = A[-val.inds]
		Y.trn = Y[-val.inds]
		W.val = W[val.inds,]
		A.val = A[val.inds]
		Y.val = Y[val.inds]

		id.trn = id[-val.inds]
		obsWeights.trn = obsWeights[-val.inds]

		if(ipcw){
			Qbar = rep(mean(Y.trn),2*(nrow(W.val)+nrow(W.trn)))
		} else {
			Qbar = suppressWarnings(Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(Y.trn,data.frame(W.trn,A=A.trn),newX=data.frame(rbind(W.val,W.val,W.trn,W.trn),A=c(rep(c(0,1),each=nrow(W.val)),rep(c(0,1),each=nrow(W.trn)))),family=family,SL.library=SL.library,obsWeights=obsWeights.trn,id=id.trn,cvControl=SuperLearner.CV.control(stratifyCV=stratifyCV&(family$family=='binomial'),V=num.folds))$SL.predict[,1]}))/num.SL.rep)
		}
		Q.val = cbind(Qbar[1:nrow(W.val)],Qbar[(nrow(W.val)+1):(2*nrow(W.val))])

		if(length(g0)==0){
			g = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(A.trn,W.trn,newX=rbind(W.val,W.trn),family=binomial(),SL.library=SL.library,obsWeights=obsWeights.trn,id=id.trn,cvControl=SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.folds))$SL.predict[,1]}))/num.SL.rep
			g.val = head(g,nrow(W.val))
			g.trn = tail(g,n=nrow(W.trn))
		} else {
			g.val = g0[val.inds]
			g.trn = g0[-val.inds]
		}

		Z.trn = (2*A.trn-1)/(A.trn*g.trn + (1-A.trn)*(1-g.trn)) * (Y.trn-(1-A.trn)*Qbar[(2*nrow(W.val)+1):(2*nrow(W.val)+nrow(W.trn))]-A.trn*Qbar[(2*nrow(W.val)+nrow(W.trn)+1):(2*nrow(W.val)+2*nrow(W.trn))]) + Qbar[(2*nrow(W.val)+nrow(W.trn)+1):(2*nrow(W.val)+2*nrow(W.trn))]-Qbar[(2*nrow(W.val)+1):(2*nrow(W.val)+nrow(W.trn))]
		tmp = Reduce('+',lapply(1:num.SL.rep,function(i){
			SL = SuperLearner(Z.trn,W.trn,newX=W.val,SL.library=SL.library,obsWeights=obsWeights.trn,id=id.trn,cvControl=SuperLearner.CV.control(validRows=CVFolds(length(Y.trn),id.trn,Y.trn,SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.folds)),V=num.folds))
			cbind(SL$SL.predict[,1],SL$library.predict)}))/num.SL.rep
		blip.est = tmp[,1]
		lib.blip.est = tmp[,-1]
		return(cbind(Q.val,blip.est,g.val,lib.blip.est))
	}))[order(unlist(folds)),]
	Q.est = ests[,1:2]
	blip.est = ests[,3]
	g.est = ests[,4]
	lib.blip.est = ests[,-(1:4)]
	return(list(Q.est=Q.est,blip.est=blip.est,g.est=g.est,lib.blip.est=lib.blip.est))
}

#' CV-TMLE Estimating Impact of Treating Optimal Subgroup
#' 
#' @description
#' This function estimates the impact of treating the optimal subgroup versus giving the opposite treatment using a CV-TMLE.
#' @details
#' Cross-validated targeted minimum loss-based estimator (CV-TMLE) for the impact of treating the optimal subgroup versus either treating no one or randomizing the treatment via a fair coin toss
#'
#' When \code{bothactive} is FALSE, \code{sig.trunc} is useful if the treatment effect is negative (or small) in all strata of covariates, since in these cases we expect that the variance estimate will converge to zero and coverage of the non-truncated confidence interval may suffer. The truncation ensures that the standard deviation estimate is never too small so that coverage of the lower bound is at worst conservative. Coverage of the upper bound also relies on being able to estimate the optimal subgroup well in terms of mean outcome (see the cited papers).
#'
#' We do not have any theoretical justication for the CV-TMLE confidence interval when the treatment effect falls on the decision boundary with positive probability (decision boundary is zero), though we have seen that it performs well in simulations.
#' @keywords CV-TMLE, subgroup
#' @param W data frame with observations in the rows and baseline covariates used to form the subgroup in columns.
#' @param A binary treatment vector, where 1 indicates that an individual was treated.
#' @param Y real-valued outcome for which large values are preferred (if use relative risk contrast, then Y should be an indicator of the absence of an adverse event, and the relative risk returned is the relative risk of the adverse event).
#' @param SL.library SuperLearner library (see documentation for \code{SuperLearner} in the corresponding package) used to estimate outcome regression, treatment mechanism (if unknown), and conditional average treatment effect function.
#' @param bothactive If TRUE, then returns the impact of treating the optimal subgroup versus randomizing the treatment via a fair coin toss. Otherwise, returns an estimate of the impact of treating the optimal subgroup versus treating no one.
#' @param kappa maximum allowable probability of a randomly drawn individual belonging to the optimal subgroup. The default of 1 indicates no constraint. Only relevant if bothactive is FALSE
#' @param g0 if known (as in a randomized controlled trial), a vector of probabilities of treatment A=1 given covariates. If \code{NULL}, \code{SuperLearner} will be used to estimate these probabilities.
#' @param family \code{binomial()} if outcome bounded in [0,1], or \code{gaussian()} otherwise.
#' @param sig.trunc value at which the standard deviation estimate is truncated (only used if bothactive is TRUE)
#' @param alpha confidence level for returned confidence interval set to (1-alpha)*100\%.
#' @param num.folds number of folds to use in cross-validation step of the CV-TMLE.
#' @param num.SL.rep number of super-learner repetitions (increasing this number should make the algorithm more stable across seeds).
#' @param num.est.rep number of repetitions of estimator, minimizing variation over cross-validation fold assignment (increasing this number should make the algorithm more stable across seeds)
#' @param id optional cluster identification variable. Will ensure rows with same id remain in same validation fold each time cross-validation used
#' @param obsWeights observation weights
#' @param stratifyCV stratify validation folds by event counts (does this for estimation of outcome regression, treatment mechanism, and conditional average treatment effect function). Useful for rare outcomes
#' @param RR estimates relative risk (TRUE) or additive contrast (FALSE) between the mean outcome under optimal versus randomizing treatment via a fair coin toss. For relative risk, estimates the additive outcome of Y not occurring (since throughut we assume Y is beneficial)
#' @param ipcw if TRUE, then does not estimate outcome regression (just sets it to zero for all covariate-treatment combinations)
#' @param lib.ests Also return estimates based on candidate optimal rule estimates in the super-learner library
#' @param verbose give status updates
#' @return a list containing
#' \item{est}{Vector containing estimates of the impact of treating the optimal subgroup. Items in the vector correspond to different choices of algorithms for estimating the optimal treatment rule (if \code{lib.ests} is FALSE, only returns SuperLearner estimate).}
#' \item{ci}{Matrix containing confidence intervals for the impact of treating the optimal subgroup. Left column contains lower bounds, right column contains upper bounds. Rows correspond to different choices of algorithms for estimating the optimal treatment rule (if \code{lib.ests} is FALSE, only returns SuperLearner estimate).}
#' @references
#' ``Evaluating the Impact of Treating the Optimal Subgroup,'' technical report to be released soon.
#' 
#' M. J. van der Laan and A. R. Luedtke, ``Targeted learning of the mean outcome under an optimal dynamic treatment rule,'' \emph{Journal of Causal Inference}, vol. 3, no. 1, pp. 61-95, 2015.
#' @examples
#' SL.library = c('SL.mean','SL.glm')
#' Qbar = function(a,w){plogis(a*(w$W1>0))}
#' n = 500
#' W = data.frame(W1=rnorm(n),W2=rnorm(n),W3=rnorm(n),W4=rnorm(n))
#' A = rbinom(n,1,1/2)
#' Y = rbinom(n,1,Qbar(A,W))
#' 
#' sg.cvtmle(W,A,Y,bothactive=TRUE,SL.library=SL.library,num.SL.rep=2,num.folds=5,family=binomial())
#' sg.cvtmle(W,A,Y,bothactive=TRUE,SL.library=SL.library,num.SL.rep=2,num.folds=5,family=binomial(),stratifyCV=TRUE)
#' sg.cvtmle(W,A,Y,bothactive=TRUE,SL.library=SL.library,num.SL.rep=2,num.folds=5,family=binomial(),id=rep(1:(n/2),2),obsWeights=1+3*runif(n))
#' 
#' sg.cvtmle(W,A,Y,bothactive=FALSE,SL.library=SL.library,num.SL.rep=2,num.folds=5,family=binomial(),sig.trunc=0.001)
#' sg.cvtmle(W,A,Y,bothactive=FALSE,SL.library=SL.library,num.SL.rep=2,num.folds=5,kappa=0.25,family=binomial())
#' sg.cvtmle(W,A,Y,bothactive=FALSE,SL.library=SL.library,num.SL.rep=2,num.folds=5,family=binomial(),stratifyCV=TRUE)
#' sg.cvtmle(W,A,Y,bothactive=FALSE,SL.library=SL.library,num.SL.rep=2,num.folds=5,family=binomial(),id=rep(1:(n/2),2),obsWeights=rbinom(n,1,2/3)*runif(n))
#' @export

sg.cvtmle = function(W,A,Y,SL.library,bothactive=FALSE,kappa=1,g0=NULL,family=binomial(),sig.trunc=1e-10,alpha=0.05,num.folds=10,num.SL.rep=5,num.est.rep=5,id=NULL,obsWeights=NULL,stratifyCV=FALSE,RR=FALSE,ipcw=FALSE,lib.ests=FALSE,verbose=TRUE){
	require(SuperLearner)
	if(!(sum(obsWeights) %in% c(0,length(Y)))){
		warning('Renormalizing obsWeights so they sum to sample size.')
		obsWeights = obsWeights/mean(obsWeights)
	}

	fits = lapply(1:num.est.rep,function(i){
		message(paste0('Running estimator iteration ',i,'.'))
		init.ests = .sgcvtmle.preprocess(W,A,Y,SL.library,g0=g0,family=family,num.folds=num.folds,num.SL.rep=num.SL.rep,id=id,obsWeights=obsWeights,stratifyCV=stratifyCV,ipcw=ipcw)
		Q.est = init.ests$Q.est
		blip.est = init.ests$blip.est
		g.est = init.ests$g.est
		lib.blip.est = init.ests$lib.blip.est

		if(bothactive){
			out = .sgtmle.bothactive(W,A,Y,Q.est,blip.est,g.est,alpha=alpha,obsWeights=obsWeights,RR=RR)
		} else {
			out = .sgtmle(W,A,Y,Q.est,blip.est,g.est,kappa=kappa,sig.trunc=sig.trunc,alpha=alpha,obsWeights=obsWeights,RR=RR)
		}
		if(lib.ests){
			out.lib = lapply(1:ncol(lib.blip.est),function(i){
				if(bothactive){
					return(.sgtmle.bothactive(W,A,Y,Q.est,lib.blip.est[,i],g.est,alpha=alpha,obsWeights=obsWeights,RR=RR))
				} else {
					return(.sgtmle(W,A,Y,Q.est,lib.blip.est[,i],g.est,kappa=kappa,sig.trunc=sig.trunc,alpha=alpha,obsWeights=obsWeights,RR=RR))
				}
			})
			return(list(ests=c(out$est,sapply(out.lib,function(x){x$est})),vars=c(var(out$ic),sapply(out.lib,function(x){var(x$ic)})),algs=c('SuperLearner',colnames(lib.blip.est))))
		} else {
			return(list(ests=out$est,vars=var(out$ic),algs='SuperLearner'))
		}
	})

	nms = fits[[1]]$algs

	est.mat = do.call(rbind,lapply(fits,function(x){x$ests}))
	rownames(est.mat) = paste0('Repetition ',1:nrow(est.mat))
	est = colMeans(est.mat)
	se = sqrt(colMeans(do.call(rbind,lapply(fits,function(x){x$vars})))/nrow(W))
	ci = cbind(est - 1.96*se,est + 1.96*se)
	colnames(ci) = c('lb','ub')

	colnames(est.mat) <- rownames(ci) <- names(est) <- nms

	return(list(est=est,ci=ci,est.mat=est.mat))
}
