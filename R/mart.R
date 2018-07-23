
#' Stabilized One-Step Estimating Impact of Treating Optimal Subgroup
#' 
#' @description
#' This function estimates the impact of treating the optimal subgroup using a stabilized one-step estimator.
#' @details
#' Stabilized one-step estimator for the impact of treating the optimal subgroup. Users will generally find this estimator more conservative (biased down) than the CV-TMLE. While asymptotic results and prior simulations show that this estimator performs as well as the CV-TMLE in large samples, the CV-TMLE outperforms this estimator provided the necessary regularity conditions hold for the CV-TMLE. Unlike for the CV-TMLE confidence interval, we have theoretical justication for the stabilized one-step confidence interval when the treatment effect is null with positive probability.
#' 
#' \code{sig.truncs} is useful if the treatment effect is negative (or small) in all strata of covariates, since in these cases we expect that the variance estimate will converge to zero and coverage of the non-truncated confidence interval may suffer. The truncation ensures that the standard deviation estimate is never too small so that coverage of the lower bound is at worst conservative. Coverage of the upper bound also relies on being able to estimate the optimal subgroup well in terms of mean outcome (see the cited papers).
#' 
#' Note: Have not implemented resource-constrained stabilized one-step estimator.
#' @keywords stabilized one-step, subgroup
#' @param W data frame with observations in the rows and baseline covariates used to form the subgroup in columns.
#' @param A binary treatment vector, where 1 indicates that an individual was treated.
#' @param Y real-valued outcome.
#' @param SL.library SuperLearner library (see documentation for \code{SuperLearner} in the corresponding package) used to estimate outcome regression, treatment mechanism (if unknown), and conditional average treatment effect function.
#' @param g0 if known (as in a randomized controlled trial), a vector of probabilities of treatment A=1 given covariates. If \code{NULL}, \code{SuperLearner} will be used to estimate these probabilities.
#' @param family \code{binomial()} if outcome bounded in [0,1], or \code{gaussian()} otherwise.
#' @param sig.truncs vector of values at which the standard deviation estimate is truncated. See \code{Details}.
#' @param alpha confidence level for returned confidence interval set to (1-alpha)*100\%.
#' @param elln number of observations to set aside for the first evaluation of the nuisance functions. Recommend around nrow(W)/10, unless willing to significantly slow down run time by setting num.chunks to (n-elln)/k for a small constant k (e.g. 1 or 2) -- in that case can set elln=20 or so.
#' @param num.chunks number of chunks to use when computing the nuisance functions -- using more chunks should slightly improve statistical properties but slow down the runtime. Allowed values between 1 and \code{nrow(W)-elln}. Recommend at least 10.
#' @return a list containing
#' \item{est}{Vector containing estimates of the impact of treating the optimal subgroup. Items in the vector correspond to different choices of sig.truncs (order maintained).}
#' \item{ci}{Matrix containing confidence intervals for the impact of treating the optimal subgroup. Left column contains lower bounds, right column contains upper bounds. Rows correspond to different choices of sig.truncs (order maintained).}
#' @references 
#' ``Evaluating the Impact of Treating the Optimal Subgroup,'' technical report to be released soon.
#' 
#' A. R. Luedtke and M. J. van der Laan, ``Statistical inference for the mean outcome under a possibly non-unique optimal treatment strategy,'' \emph{Annals of Statistics}, vol. 44, no. 2, pp. 713-742, 2016.
#' @examples
#' SL.library = c('SL.mean','SL.glm')
#' Qbar = function(a,w){plogis(0.2*a*(w$W1>0))}
#' n = 1000
#' W = data.frame(W1=rnorm(n),W2=rnorm(n),W3=rnorm(n),W4=rnorm(n))
#' A = rbinom(n,1,1/2)
#' Y = rnorm(n,Qbar(A,W),sd=0.025)
#' 
#' sg.mart(W,A,Y,SL.library=SL.library,family=binomial(),sig.truncs=c(0.001,0.1))
#' @export

sg.mart = function(W,A,Y,SL.library,g0=NULL,family=binomial(),sig.truncs=0.1,alpha=0.05,elln=ceiling(nrow(W)/10),num.chunks=10,num.SL.rep=5) {
	require(SuperLearner)
	n = nrow(W)

	mt = Reduce('+',lapply(0:(num.chunks-1),function(i){
		old.inds = 1:(elln + i*floor((n-elln)/num.chunks) + max(i-(num.chunks-(n-elln)%%num.chunks),0))
		new.inds = (max(old.inds) + 1):min(max(old.inds) + floor((n-elln)/num.chunks) + (i>=num.chunks-(n-elln)%%num.chunks),n)

		W.old = W[old.inds,]
		A.old = A[old.inds]
		Y.old = Y[old.inds]

		W.new = W[new.inds,]
		A.new = A[new.inds]
		Y.new = Y[new.inds]

		tmp = suppressWarnings(Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(Y.old,data.frame(W.old,A=A.old),newX=data.frame(rbind(W.new,W.new,W.old,W.old),A=c(rep(c(0,1),each=nrow(W.new)),rep(c(0,1),each=nrow(W.old)))),family=family,SL.library=SL.library)$SL.predict[,1]}))/num.SL.rep)
		# if(sum(tmp$errorsInCVLibrary | tmp$errorsInLibrary)>0) tmp = suppressWarnings(SuperLearner(Y.old,data.frame(W.old,A=A.old),newX=data.frame(rbind(W.new,W.new,W.old,W.old),A=c(rep(c(0,1),each=nrow(W.new)),rep(c(0,1),each=nrow(W.old)))),family=family,SL.library=SL.library[-which(tmp$errorsInCVLibrary | tmp$errorsInLibrary)]))
		Qbar.0W.new = tmp[1:nrow(W.new)]
		Qbar.1W.new = tmp[(nrow(W.new)+1):(2*nrow(W.new))]
		Qbar.0W.old = tmp[(2*nrow(W.new) + 1):(2*nrow(W.new) + nrow(W.old))]
		Qbar.1W.old = tmp[(2*nrow(W.new) + nrow(W.old)+1):(2*nrow(W.new) + 2*nrow(W.old))]
		rm(tmp)

		if(length(g0)==0){
			g = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(A.old,W.old,newX=rbind(W.old,W.new),family=binomial(),SL.library=SL.library)$SL.predict[,1]}))/num.SL.rep
			g.old = head(g,n=nrow(W.old))
			g.new = tail(g,n=nrow(W.new))
		} else {
			g.old = g0[old.inds]
			g.new = g0[new.inds]
		}


		Z.old = (2*A.old-1)/(A.old*g.old + (1-A.old)*(1-g.old)) * (Y.old-(1-A.old)*Qbar.0W.old-A.old*Qbar.1W.old) + Qbar.1W.old-Qbar.0W.old
		uz = unique(Z.old)
		if(length(uz)>1){
			tmp = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(Z.old,W.old,newX=rbind(W.new,W.old),family=gaussian(),SL.library=SL.library)$SL.predict[,1]}))/num.SL.rep
			# if(sum(tmp$errorsInCVLibrary | tmp$errorsInLibrary)>0) tmp = SuperLearner(Z.old,W.old,newX=rbind(W.new,W.old),family=gaussian(),SL.library=SL.library[-which(tmp$errorsInCVLibrary | tmp$errorsInLibrary)])
			blip.new = tmp[1:nrow(W.new)]
			blip.old = tmp[(nrow(W.new)+1):(nrow(W.new)+nrow(W.old))]
			rm(tmp)
		} else {
			blip.new = rep(uz,nrow(W.new))
			blip.old = rep(uz,nrow(W.old))
		}

		return(do.call(rbind,lapply(sig.truncs,function(sig.trunc){
			dn.old = as.numeric(blip.old>0)
			dn.new = as.numeric(blip.new>0)

			IC.new = dn.new * ((2*A.new-1)/(A.new*g.new + (1-A.new)*(1-g.new))  * (Y.new-A.new*Qbar.1W.new - (1-A.new)*Qbar.0W.new) + Qbar.1W.new - Qbar.0W.new)
			IC.old = dn.old * ((2*A.old-1)/(A.old*g.old + (1-A.old)*(1-g.old))  * (Y.old-A.old*Qbar.1W.old - (1-A.old)*Qbar.0W.old) + Qbar.1W.old - Qbar.0W.old)

			wt = 1/max(sd(IC.old),sig.trunc)
			IC.sum = sum(IC.new * wt)
			c(IC.sum,nrow(W.new)*wt)
		})))
	}),0)

	m = n-elln
	delta = alpha/2

	est = mt[,1]/mt[,2]
	sig.bar = m/mt[,2]

	ci = cbind(est-qnorm(1-alpha/2)*sig.bar/sqrt(m),est+qnorm(1-alpha/2)*sig.bar/sqrt(m))

	return(list(est=est,ci=ci))
}

