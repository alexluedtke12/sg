# Standard deviation truncated at sig.trunc
# g.est is vector giving estimated probability of A=1 given W

.sgtmle = function(W,A,Y,Q.est,blip.est,g.est,baseline.probs,family=binomial(),kappa=1,sig.trunc=0.1,alpha=0.05,trunc.Q=c(0.005,0.995),obsWeights=NULL,RR=FALSE){
	n = nrow(W)
	if(is.null(obsWeights)) obsWeights = rep(1,length(Y))	
	if(family$family=='binomial') for(i in 1:2){Q.est[,i] = pmin(pmax(Q.est[,i],trunc.Q[1]),trunc.Q[2])}

	offs = family$linkfun((A==1)*Q.est[,2] + (A==0)*Q.est[,1] + (A!=0 & A!=1)/2)

	tau = max(quantile(blip.est,1-kappa,type=1),0)

	if(mean(blip.est==tau)==0) {
		g.star = as.numeric(blip.est>tau)
	} else g.star = (blip.est==tau)*(tau>0)*(kappa - mean(blip.est>tau))/mean(blip.est==tau) + (blip.est>tau)

	if((sum(g.star)>0) | ((sum(g.star)==0) & !(baseline.probs[2]==0 & baseline.probs[1]==1))){
		H = (A==1)*(g.star-baseline.probs[2])/g.est[,2] + (A==0)*(1-g.star - baseline.probs[1])/g.est[,1]
		eps = suppressWarnings(coef(glm(Y ~ -1 + offset(offs) + H,family=family,weights=obsWeights)))

		Q0 = family$linkinv(family$linkfun(Q.est[,1]) + eps*I(1-g.star - baseline.probs[1])/g.est[,1])
		Q1 = family$linkinv(family$linkfun(Q.est[,2]) + eps*I(g.star-baseline.probs[2])/g.est[,2])

		est.add = mean(Q1*(g.star-baseline.probs[2]) + Q0*(1-g.star - baseline.probs[1]))
		ic.add = obsWeights * (H*(Y-((A==1)*Q1 + (A==0)*Q0)) - tau*(g.star-kappa) + (g.star-baseline.probs[2])*Q1+(1-g.star-baseline.probs[1])*Q0 - est.add)

		if(RR){ # Get a targeted estimate of mean outcome under baseline.probs, and then add this to the contrast from the previous
			offs.bp = family$linkfun(baseline.probs[1]*Q.est[,1]+baseline.probs[2]*Q.est[,2])
			H.bp = baseline.probs[1]*(A==0)/g.est[,1] + baseline.probs[2]*(A==1)/g.est[,2]
			eps.bp = suppressWarnings(coef(glm(Y ~ -1 + offset(offs.bp) + H.bp,family=family,weights=obsWeights)))
			Qbp.star = family$linkinv(offs.bp + eps.bp * (baseline.probs[1]/g.est[,1] + baseline.probs[2]/g.est[,2]))
			EYbp.star = mean(Qbp.star)
			est = (1 - (EYbp.star + est.add))/(1-EYbp.star) # recall: relative risk of 1-Y, since Y beneficial
			ic = -(ic.add/(1 - (EYbp.star + est.add)) + obsWeights * (H.bp * (Y-Qbp.star) + Qbp.star - EYbp.star) * (1/(1 - (EYbp.star + est.add))-1/(1-EYbp.star)))
			ci = exp(log(est) + c(-1,1) * qnorm(1-alpha/2)*max(sd(ic),sig.trunc)/sqrt(n))
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



# Helper function used by sg.cvtmle
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
			Qbar = do.call(cbind,lapply(0:1,function(a){suppressWarnings(
				Reduce('+',lapply(1:num.SL.rep,function(i){
						if(sum(A.trn==a)>0){
							return(SuperLearner(
								Y.trn[A.trn==a],
								data.frame(W.trn,A=A.trn)[A.trn==a,],
								newX=data.frame(rbind(W.val,W.trn),A=a),
								family=family,
								SL.library=SL.library,
								obsWeights=obsWeights.trn[A.trn==a],
								id=id.trn[A.trn==a],
								cvControl=SuperLearner.CV.control(stratifyCV=stratifyCV&(family$family=='binomial'),V=num.folds)
								)$SL.predict[,1])
						} else {
							# warning("No individuals in a training fold received one of the two treatments when estimating the outcome regression.")
							return(rep(mean(Y.trn),nrow(W.val)+nrow(W.trn)))
						}
					}))/num.SL.rep)}
						))
		}
		Q.val = Qbar[1:nrow(W.val),]

		if(length(g0)==0){
			g1 = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(as.numeric(A.trn==1),W.trn,newX=rbind(W.val,W.trn),family=binomial(),SL.library=SL.library,obsWeights=obsWeights.trn,id=id.trn,cvControl=SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.folds))$SL.predict[,1]}))/num.SL.rep
			if(length(unique(A))>2){
				g0 = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(as.numeric(A.trn==0),W.trn,newX=rbind(W.val,W.trn),family=binomial(),SL.library=SL.library,obsWeights=obsWeights.trn,id=id.trn,cvControl=SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.folds))$SL.predict[,1]}))/num.SL.rep
				if(sum(g0+g1>1)>0){
					sm = g0 + g1
					g0 = g0/sm
					g1 = g1/sm
				}
			} else {
				g0 = 1-g1
			}
			g = cbind(g0,g1)
			g.val = g[1:nrow(W.val),]
			g.trn = g[(nrow(W.val)+1):nrow(W),]
		} else {
			g.val = g0[val.inds,]
			g.trn = g0[-val.inds,]
		}

		Z.trn = ((A.trn==1)-(A.trn==0))/((A.trn==1)*g.trn[,2] + (A.trn==0)*g.trn[,1] + (A.trn!=0 & A.trn!=1)) * (Y.trn-(A.trn==0)*Qbar[(nrow(W.val)+1):nrow(W),1]-(A.trn==1)*Qbar[(nrow(W.val)+1):nrow(W),2]) + Qbar[(nrow(W.val)+1):nrow(W),2]-Qbar[(nrow(W.val)+1):nrow(W),1]
		tmp = Reduce('+',lapply(1:num.SL.rep,function(i){
			SL = SuperLearner(Z.trn,W.trn,newX=W.val,SL.library=SL.library,obsWeights=obsWeights.trn,id=id.trn,cvControl=SuperLearner.CV.control(validRows=CVFolds(length(Y.trn),id.trn,Y.trn,SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.folds)),V=num.folds),method="method.NNLS2")
			cbind(SL$SL.predict[,1],SL$library.predict)}))/num.SL.rep
		blip.est = tmp[,1]
		lib.blip.est = tmp[,-1]
		return(cbind(Q.val,blip.est,g.val,lib.blip.est))
	}))[order(unlist(folds)),]
	Q.est = ests[,1:2]
	blip.est = ests[,3]
	g.est = ests[,4:5]
	lib.blip.est = ests[,-(1:5)]
	return(list(Q.est=Q.est,blip.est=blip.est,g.est=g.est,lib.blip.est=lib.blip.est))
}

#' CV-TMLE Estimating Impact of Treating Optimal Subgroup
#' 
#' @description
#' This function estimates the impact of treating the optimal subgroup versus giving the opposite treatment using a CV-TMLE.
#' @details
#' Cross-validated targeted minimum loss-based estimator (CV-TMLE) for the impact of treating the optimal subgroup versus either a baseline treatment strategy
#'
#' \code{sig.trunc} is useful if the treatment effect is negative (or small) in all strata of covariates and \code{baseline.probs}=c(1,0), since in these cases we expect that the variance estimate will converge to zero and coverage of the non-truncated confidence interval may suffer. The truncation ensures that the standard deviation estimate is never too small so that coverage of the lower bound is at worst conservative. Coverage of the upper bound also relies on being able to estimate the optimal subgroup well in terms of mean outcome (see the cited papers).
#'
#' We do not have any theoretical justication for the CV-TMLE confidence interval when the treatment effect falls on the decision boundary with positive probability (decision boundary is zero), though we have seen that it performs well in simulations.
#' @keywords CV-TMLE, subgroup
#' @param W data frame with observations in the rows and baseline covariates used to form the subgroup in columns.
#' @param A a vector of treatments. This vector may have more than two levels, but the two levels of interest in the treatment assignment problem must be specified in \code{txs}.
#' @param Y real-valued outcome for which large values are preferred (if use relative risk contrast, then Y should be an indicator of the absence of an adverse event, and the relative risk returned is the relative risk of the adverse event).
#' @param SL.library SuperLearner library (see documentation for \code{SuperLearner} in the corresponding package) used to estimate outcome regression, treatment mechanism (if unknown), and conditional average treatment effect function.
#' @param txs A vector of length two indicating the two treatments of interest in A that will be used for the treatment assignment problem. Individuals with treatment \code{A} equal to a value not in \code{txs} will be treated as censored at treatment.
#' @param baseline.probs A vector of length two indicating the (stochastic) treatment rule to use as a baseline when evaluating performance of the estimated optimal treatment rule. In this treatment rule, the first treatment in \code{txs} is assigned with probability \code{baseline.probs[1]} and the second is assigned with probability \code{baseline.probs[2]}. To obtain the marginal mean outcome under the optimal treatment strategy, i.e. not contrasting against any baseline, set \code{baseline.probs=c(0,0)}.
#' @param kappa maximum allowable probability of a randomly drawn individual belonging to the optimal subgroup. The default of 1 indicates no constraint.
#' @param g0 if known (as in a randomized controlled trial), a matrix of probabilities of receiving the first tx in \code{txs} given covariates (first column) and the second treatment in \code{txs} given covariates (second column). Rows correspond to individuals with (\code{W},\code{A},\code{Y}) observed. If \code{NULL}, \code{SuperLearner} will be used to estimate these probabilities.
#' @param family \code{binomial()} if outcome bounded in [0,1], or \code{gaussian()} otherwise.
#' @param sig.trunc value at which the standard deviation estimate is truncated. 
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
#' @param init.ests.out Set this option to TRUE to return the initial SuperLearner estimates. Can be fed to a new call of this function using init.ests.in to speed up that call. E.g., useful if want to call this function at many values of \code{kappa.}
#' @param init.ests.in Can be used to feed the function the initial SuperLearner estimates from a previous call of this function (see \code{init.ests.out}). Dramatically reduces runtime. See Example below.
#' @param verbose give status updates
#' @return a list containing
#' \item{est}{Vector containing estimates of the impact of treating the optimal subgroup. Items in the vector correspond to different choices of algorithms for estimating the optimal treatment rule (if \code{lib.ests} is FALSE, only returns SuperLearner estimate).}
#' \item{ci}{Matrix containing confidence intervals for the impact of treating the optimal subgroup. Left column contains lower bounds, right column contains upper bounds. Rows correspond to different choices of algorithms for estimating the optimal treatment rule (if \code{lib.ests} is FALSE, only returns SuperLearner estimate).}
#' \item{est.mat}{Estimates across repetitions.}
#' @references
#' ``Evaluating the Impact of Treating the Optimal Subgroup,'' technical report to be released soon.
#' 
#' M. J. van der Laan and A. R. Luedtke, ``Targeted learning of the mean outcome under an optimal dynamic treatment rule,'' \emph{Journal of Causal Inference}, vol. 3, no. 1, pp. 61-95, 2015.
#' @examples
#' SL.library = c('SL.mean','SL.glm')
#' Qbar = function(a,w){plogis(a*w$W1)}
#' n = 500
#' W = data.frame(W1=rnorm(n),W2=rnorm(n),W3=rnorm(n),W4=rnorm(n))
#' A = rbinom(n,1,1/2)
#' Y = rbinom(n,1,Qbar(A,W))
#' 
#' # comparing the mean outcome under the optimal rule to the mean outcome
#' # when treating half of the population at random
#' sg.cvtmle(W,A,Y,baseline.probs=c(0.5,0.5),SL.library=SL.library,num.SL.rep=2,num.folds=5,family=binomial())
#' # same as above, but adding ids (used in CV splits) and in observation weights
#' sg.cvtmle(W,A,Y,baseline.probs=c(0.5,0.5),SL.library=SL.library,num.SL.rep=2,num.folds=5,family=binomial(),id=rep(1:(n/2),2),obsWeights=1+3*runif(n))
#' 
#' # comparing the mean outcome under the optimal rule against the mean outcome under treating no one
#' sg.cvtmle(W,A,Y,baseline.probs=c(1,0),SL.library=SL.library,num.SL.rep=2,num.folds=5,family=binomial(),sig.trunc=0.001)
#' # comparing the mean outcome under an optimal rule that treats at most 25 percent of people
#' # to the mean outcome under treating 25 percent of people at random
#' sg.cvtmle(W,A,Y,baseline.probs=c(0.75,0.25),SL.library=SL.library,num.SL.rep=2,num.folds=5,kappa=0.25,family=binomial())
#' # same as above, but adding ids (used in CV splits) and in observation weights
#' sg.cvtmle(W,A,Y,baseline.probs=c(1,0),SL.library=SL.library,num.SL.rep=2,num.folds=5,family=binomial(),id=rep(1:(n/2),2),obsWeights=rbinom(n,1,2/3)*runif(n))
#'
#' # comparing the mean outcomes under optimal rules that treats at most prop percent of people
#' out_10 = sg.cvtmle(W,A,Y,baseline.probs=c(0,0),SL.library=SL.library,num.SL.rep=2,num.folds=5,kappa=0.10,family=binomial(),init.ests.out=TRUE)
#' init.ests = out_10$init.ests
#' for(prop in seq(0.10,0.25,by=0.05)){
#'   print(paste0("Can treat ",prop," proportion of population."))
#'   out = sg.cvtmle(W,A,Y,baseline.probs=c(0,0),SL.library=SL.library,num.SL.rep=2,num.folds=5,kappa=prop,family=binomial(),init.ests.out=FALSE,init.ests.in=init.ests,verbose=FALSE)
#' 	 print(out$est)
#' }
#' @export

sg.cvtmle = function(W,A,Y,SL.library,txs=c(0,1),baseline.probs=c(0.5,0.5),kappa=1,g0=NULL,family=binomial(),sig.trunc=1e-10,alpha=0.05,num.folds=10,num.SL.rep=5,num.est.rep=5,id=NULL,obsWeights=NULL,stratifyCV=FALSE,RR=FALSE,ipcw=FALSE,lib.ests=FALSE,init.ests.out=FALSE,init.ests.in=NULL,verbose=TRUE,...){
	require(SuperLearner)

	if(any(names(list(...))=="bothactive")) {
		warning("The bothactive option is deprecated. Its functionality can now be imitated using the baseline.probs argument. To imitate bothactive=TRUE, set baseline.probs=c(0.5,0.5). To imitate bothactive=FALSE, set txs=c(0,1) and baseline.probs=c(1,0).")
	}

	if(!(sum(obsWeights) %in% c(0,length(Y)))){
		warning('Renormalizing obsWeights so they sum to sample size.')
		obsWeights = obsWeights/mean(obsWeights)
	}

	# Recode A so that first value in txs is coded as zero, second is coded as 1
	A.new = rep(NA,length(A))
	A.new[A==txs[1]] = 0
	A.new[A==txs[2]] = 1
	A.new[!(A%in%txs)] = 1e10
	A = A.new

	# if init.ests.in was used, make surehas the same length as num.SL.rep
	# (since otherwise there may be a user input error)
	if(length(init.ests.in)>0 & (length(init.ests.in)!=num.est.rep)){
		stop("Length of init.ests.in should be the same as num.est.rep")
	}

	fits = lapply(1:num.est.rep,function(i){
		if(verbose) message(paste0('Running estimator iteration ',i,'.'))
		if(length(init.ests.in)>0){
			init.ests = init.ests.in[[i]]
		} else {
			init.ests = .sgcvtmle.preprocess(W,A,Y,SL.library,g0=g0,family=family,num.folds=num.folds,num.SL.rep=num.SL.rep,id=id,obsWeights=obsWeights,stratifyCV=stratifyCV,ipcw=ipcw)
		}
		Q.est = init.ests$Q.est
		blip.est = init.ests$blip.est
		g.est = init.ests$g.est
		lib.blip.est = init.ests$lib.blip.est

		out = .sgtmle(W,A,Y,Q.est,blip.est,g.est,baseline.probs,kappa=kappa,sig.trunc=sig.trunc,alpha=alpha,obsWeights=obsWeights,RR=RR)
		if(lib.ests){
			out.lib = lapply(1:ncol(lib.blip.est),function(i){
				return(.sgtmle(W,A,Y,Q.est,lib.blip.est[,i],g.est,baseline.probs,kappa=kappa,sig.trunc=sig.trunc,alpha=alpha,obsWeights=obsWeights,RR=RR))
			})
			out.list = list(ests=c(out$est,sapply(out.lib,function(x){x$est})),vars=c(var(out$ic),sapply(out.lib,function(x){var(x$ic)})),algs=c('SuperLearner',colnames(lib.blip.est)))
		} else {
			out.list = list(ests=out$est,vars=var(out$ic),algs='SuperLearner')
		}
		if(init.ests.out){
			out.list = append(out.list,list(init.ests=init.ests))
		}
		return(out.list)
	})

	nms = fits[[1]]$algs

	est.mat = do.call(rbind,lapply(fits,function(x){x$ests}))
	rownames(est.mat) = paste0('Repetition ',1:nrow(est.mat))
	est = colMeans(est.mat)
	se = sqrt(colMeans(do.call(rbind,lapply(fits,function(x){x$vars})))/nrow(W))
	ci = cbind(est - qnorm(1-alpha/2)*se,est + qnorm(1-alpha/2)*se)
	colnames(ci) = c('lb','ub')

	colnames(est.mat) <- rownames(ci) <- names(est) <- nms

	out.list = list(est=est,ci=ci,est.mat=est.mat)
	if(init.ests.out){
		out.list = append(out.list,list(init.ests=lapply(fits,function(fit){fit$init.ests})))
	}

	return(out.list)
}
