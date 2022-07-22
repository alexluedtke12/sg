# Standard deviation truncated at sig.trunc
# g.est is vector giving estimated probability of A=1 given W

.sgtmle = function(W,A,Y,txs,Q.est,blip.est,g.est,baseline.probs,family,kappa=1,sig.trunc=0.1,alpha=0.05,trunc.Q=c(0.005,0.995),obsWeights=NULL,RR=FALSE){
	require(Hmisc)
	n = nrow(W)
	if(is.null(obsWeights)) obsWeights = rep(1,length(Y))
	obsWeights = obsWeights/mean(obsWeights)

	if(family$family=='binomial') {
		for(i in seq(txs)){Q.est[,i] = pmin(pmax(Q.est[,i],trunc.Q[1]),trunc.Q[2])}
		offs = family$linkfun(Reduce("+",lapply(seq(txs),function(i){
			(A==txs[i])*Q.est[,i]
			})))
		offs[!(A %in% txs)] = 0
	} else if(family$family=='gaussian'){
		offs = Reduce("+",lapply(seq(txs),function(i){(A==txs[i])*Q.est[,i]}))
	} else {
		stop("family must be either binomial or gaussian.")
	}

	# treatment effect of assigning first treatment in txs vs assigning the next best treatment
	# useful when kappa<1 (so that there is a constraint on this first treatment resource)
	one.vs.next.best = blip.est[,1] - apply(blip.est[,-1,drop=FALSE],1,max)
	tau = wtd.quantile(one.vs.next.best, obsWeights, type="i/n", probs=1-kappa)
	if(mean((one.vs.next.best>tau) * obsWeights)/mean(obsWeights)>kappa){
	  tau = min(one.vs.next.best[one.vs.next.best>tau])
	}

	# if one.vs.next.best<tau, then the individual will not be treated with tx 1 in the presence of the resource constraint
	blip.est[one.vs.next.best<tau,1] = -Inf

	# probability of always treating someone with treatment 1 (one.vs.next.best>tau)
	alway.trt1.prob = mean((one.vs.next.best>tau) * obsWeights)
	# probability of sometimes treating someone with treatment 1 (one.vs.next.best==tau) -- prob determined by resource constraint
	sometimes.trt1.prob = mean((one.vs.next.best==tau) * obsWeights)

	# vector of optimal treatment probabilities
	# (allowed to be stochastic to deal with resource constraint. typically is deterministic)
	g.star = t(apply(blip.est,1,function(xx){
		wm = which.max(xx)
		out = rep(0,length(xx))
		if((wm==1) & (xx[1]-max(xx[-1])==tau)){
			# only assign treatment 1 with the probability that is allowed by the resource constraint
			out[1] = (tau>0)*(kappa - alway.trt1.prob)/sometimes.trt1.prob
			# otherwise assign the second best treatment
			out[which.max(xx[-1])+1] = 1-out[1]
		} else {
			out[wm] = 1
		}
		return(out)
		}))

	H = Reduce("+",lapply(seq(txs),function(i){
			a = txs[i]
			(A==a)*(g.star[,i]-baseline.probs[i])/g.est[,i]
		}))

	eps = suppressWarnings(coef(glm(Y ~ -1 + offset(offs) + H,family=family,weights=obsWeights)))

	Q = do.call(cbind,lapply(seq(txs),function(i){
		family$linkinv(family$linkfun(Q.est[,i]) + eps*I(g.star[,i]-baseline.probs[i])/g.est[,i])
		}))

	QA = Reduce("+",lapply(seq(txs),function(i){(A==txs[i])*Q[,i]}))

	Q.contrast = Reduce("+",lapply(seq(txs),function(i){
		Q[,i]*(g.star[,i]-baseline.probs[i])
		}))

	est.add = mean(obsWeights * Q.contrast)

	ic.add = obsWeights * (H*(Y-QA) - max(tau,0)*(g.star[,1]-kappa) + Q.contrast - est.add)

	if(RR){ # Get a targeted estimate of mean outcome under baseline.probs, and then add this to the contrast from the previous
		offs.bp = family$linkfun(Reduce("+",lapply(seq(txs),function(i){baseline.probs[i]*Q.est[,i]})))
		H.bp = Reduce("+",lapply(seq(txs),function(i){baseline.probs[i]*(A==txs[i])/g.est[,i]}))
		eps.bp = suppressWarnings(coef(glm(Y ~ -1 + offset(offs.bp) + H.bp,family=family,weights=obsWeights)))
		Qbp.star = family$linkinv(offs.bp + eps.bp * Reduce("+",lapply(seq(txs),function(i){baseline.probs[i]/g.est[,i]})))
		EYbp.star = mean(obsWeights * Qbp.star)
		est = (1 - (EYbp.star + est.add))/(1-EYbp.star) # recall: relative risk of 1-Y, since Y beneficial
		ic = -(ic.add/(1 - (EYbp.star + est.add)) + obsWeights * (H.bp * (Y-Qbp.star) + Qbp.star - EYbp.star) * (1/(1 - (EYbp.star + est.add))-1/(1-EYbp.star)))
		ci = exp(log(est) + c(-1,1) * qnorm(1-alpha/2)*max(sd(ic),sig.trunc)/sqrt(n))
	} else {
		est = est.add
		ic = ic.add
		ci = c(est-qnorm(1-alpha/2)*max(sd(ic),sig.trunc)/sqrt(n),est+qnorm(1-alpha/2)*max(sd(ic),sig.trunc)/sqrt(n))
	}
	return(list(est=est,ci=ci,ic=ic))
}



# Helper function used by sg.cvtmle
# Inputs the same as for those functions

.sgcvtmle.preprocess = function(W,A,Y,Delta,SL.library,OR.SL.library,prop.SL.library,missingness.SL.library,txs,family,g0=NULL,Q0=NULL,num.folds=10,num.SL.rep=5,id=NULL,folds=NULL,obsWeights=NULL,stratifyCV=FALSE,SL.method="method.NNLS2"){
	require(SuperLearner)

	# Recode missing Y values to 0
	Y[Delta==0] = 0

	n = nrow(W)
	if(!is.null(id) & stratifyCV==TRUE) stop('Stratified sampling with id not currently implemented.')
	if(is.null(folds)){
		folds = CVFolds(n,id,Y,SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.folds))
	}

	ests = lapply(seq(folds),function(k){
		val.inds = folds[[k]]

		W.trn = W[-val.inds,,drop=FALSE]
		A.trn = A[-val.inds]
		Y.trn = Y[-val.inds]
		Delta.trn = Delta[-val.inds]
		W.val = W[val.inds,,drop=FALSE]
		A.val = A[val.inds]
		Y.val = Y[val.inds]
		Delta.val = Delta[val.inds]

		id.trn = id[-val.inds]
		obsWeights.trn = obsWeights[-val.inds]

		if(length(g0)==0){
			if((length(txs)==2) & all(Delta.trn==1)){ # no need to fit propensity for both treatments if there are only two treatments and no outcomes are missing (probabilities add to 1)
				g1 = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(as.numeric(A.trn==txs[1]),W.trn,newX=rbind(W.val,W.trn),family=binomial(),SL.library=prop.SL.library,obsWeights=obsWeights.trn,id=id.trn,cvControl=SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.folds))$SL.predict[,1]}))/num.SL.rep
				g = cbind(g1,1-g1)
			} else {
				g = do.call(cbind,lapply(txs,function(a){
						Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(as.numeric(A.trn==a),W.trn,newX=rbind(W.val,W.trn),family=binomial(),SL.library=prop.SL.library,obsWeights=obsWeights.trn,id=id.trn,cvControl=SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.folds))$SL.predict[,1]}))/num.SL.rep
					}))
				if(any(rowSums(g)>1)){ # if the estimated treatment probs sum to more than 1, then normalize
					sm = rowSums(g)
					inds = which(sm>1)
					for(i in 1:ncol(g)){
						g[inds,i] = g[inds,i]/sm[inds]
					}
				}
			}
			g.val = g[1:nrow(W.val),]
			g.trn = g[(nrow(W.val)+1):nrow(W),]
		} else {
			g.val = g0[val.inds,]
			g.trn = g0[-val.inds,]
		}

		if(all(Delta.trn==1)){
			g.delta = rep(1,length(Delta))
		} else {
			g.delta = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(Delta.trn,data.frame(W.trn,A=A.trn),newX=rbind(data.frame(W.val,A=A.val),data.frame(W.trn,A=A.trn)),family=binomial(),SL.library=missingness.SL.library,obsWeights=obsWeights.trn,id=id.trn,cvControl=SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.folds))$SL.predict[,1]}))/num.SL.rep
		}
		g.delta.val = g.delta[1:nrow(W.val)]
		g.delta.trn = g.delta[(nrow(W.val)+1):nrow(W)]

		# recode A so that it's set to missing (some value not in txs -- max(txs)+1 will work)
		# whenever Delta is set to missing
		A.val[Delta.val==0] = max(txs) + 1
		A.trn[Delta.trn==0] = max(txs) + 1
		# modify treatment probability to also account for missingness probability
		g.val = g.val * do.call(cbind,lapply(1:length(txs),function(i){g.delta.val}))
		g.trn = g.trn * do.call(cbind,lapply(1:length(txs),function(i){g.delta.trn}))

		if(length(Q0)==0){
			# fit outcome regression
			Qbar = do.call(cbind,lapply(txs,function(a){suppressWarnings(
				Reduce('+',lapply(1:num.SL.rep,function(i){
						if(sum(A.trn==a)>0){
							return(SuperLearner(
								Y.trn[A.trn==a],
								data.frame(W.trn,A=A.trn)[A.trn==a,],
								newX=data.frame(rbind(W.val,W.trn),A=a),
								family=family,
								SL.library=OR.SL.library,
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
		} else {
			Qbar = rbind(Q0[[k]][val.inds,],Q0[[k]][-val.inds,])
		}
		Q.val = Qbar[1:nrow(W.val),]

		tmp = lapply(seq(txs),function(i){
			a = txs[i]
			Z.trn = (A.trn==a)/((A.trn==a)*g.trn[,i] + (A.trn!=a)) * (Y.trn-Qbar[(nrow(W.val)+1):nrow(W),i]) + Qbar[(nrow(W.val)+1):nrow(W),i] - Delta.trn*Y.trn/g.delta.trn

			Reduce('+',lapply(1:num.SL.rep,function(j){
				SL = SuperLearner(Z.trn,W.trn,newX=W.val,SL.library=SL.library,obsWeights=obsWeights.trn,id=id.trn,cvControl=SuperLearner.CV.control(validRows=CVFolds(length(Y.trn),id.trn,Y.trn,SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.folds)),V=num.folds),method=SL.method)
				cbind(SL$SL.predict[,1],SL$library.predict)}))/num.SL.rep
		})

		blip.est = do.call(cbind,lapply(tmp,function(tx_fit){tx_fit[,1]}))
		# list of blip.est-like objects, one for each learning in the SL library
		lib.blip.est = lapply(seq(SL.library),function(i){
			do.call(cbind,lapply(tmp,function(tx_fit){tx_fit[,i+1]}))
			})

		return(list(Q.val=Q.val,blip.est=blip.est,g.val=g.val,lib.blip.est=lib.blip.est))
	})

	Q.est = do.call(rbind,lapply(ests,function(fold_ests){fold_ests$Q.val}))[order(unlist(folds)),]
	blip.est = do.call(rbind,lapply(ests,function(fold_ests){fold_ests$blip.est}))[order(unlist(folds)),]
	g.est = do.call(rbind,lapply(ests,function(fold_ests){fold_ests$g.val}))[order(unlist(folds)),]
	lib.blip.est = lapply(seq(SL.library),function(i){do.call(rbind,lapply(ests,function(fold_ests){fold_ests$lib.blip.est[[i]]}))[order(unlist(folds)),]})

	return(list(Q.est=Q.est,blip.est=blip.est,g.est=g.est,lib.blip.est=lib.blip.est))
}

#' CV-TMLE Estimating Impact of Treating Optimal Subgroup

sg.cvtmle = function(W,A,Y,SL.library,Delta=rep(1,length(A)),OR.SL.library=SL.library,prop.SL.library=SL.library,missingness.SL.library=SL.library,txs=c(0,1),baseline.probs=c(0.5,0.5),kappa=1,g0=NULL,Q0=NULL,family=binomial(),sig.trunc=1e-10,alpha=0.05,num.folds=10,num.SL.rep=5,SL.method="method.NNLS2",num.est.rep=5,id=NULL,folds=NULL,obsWeights=NULL,stratifyCV=FALSE,RR=FALSE,lib.ests=FALSE,init.ests.out=FALSE,init.ests.in=NULL,verbose=TRUE,...){
	require(SuperLearner)

	if(any(names(list(...))=="bothactive")) {
		warning("The bothactive option is deprecated. Its functionality can now be imitated using the baseline.probs argument. To imitate bothactive=TRUE, set baseline.probs=c(0.5,0.5). To imitate bothactive=FALSE, set txs=c(0,1) and baseline.probs=c(1,0).")
	}

	if(any(names(list(...))=="ipcw")) {
		warning("The ipcw argument is deprecated.")
	}

	if(!(sum(obsWeights) %in% c(0,length(Y)))){
		warning('Renormalizing obsWeights so they sum to sample size.')
		obsWeights = obsWeights/mean(obsWeights)
	}

	if(length(txs)!=length(unique(txs))){
		stop("txs should not contain any duplicate entries.")
	}

	if(length(txs)==1){
		stop("txs should contain at least two entries.")
	}

	# if init.ests.in was used, make sure has the same length as num.SL.rep
	# (since otherwise there may be a user input error)
	if(length(init.ests.in)>0 & (length(init.ests.in)!=num.est.rep)){
		stop("Length of init.ests.in should be the same as num.est.rep")
	}

	# if folds is specified, then ignore the value of num.est.rep and only use the
	# supplied folds (corresponds to a particular realization of the folds that could
	# have been used if num.est.rep had been 1)
	if(num.est.rep>1 & !is.null(folds)){
		warning("Because the folds input was specified, automatically setting num.est.rep to 1.")
		num.est.rep = 1
	}

	# if folds is specified, then stratifyCV is not respected for the purpose of choosing folds in
	# the outer optimization. Instead, the user should ensure that the supplied folds are
	# stratified by event counts if desired
	if(stratifyCV & !is.null(folds)){
		warning("Because the folds input was specified, stratifyCV will only be used to stratify inner cross-validation used to estimate nuisance functions. If stratification based on event count is desired for choosing the folds in the outer layer of cross-validation, then the user should ensure that the supplied folds respect these event counts.")
	}

	# if folds is specified, then num.folds must match the number of folds in folds. If it doesn't,
	# it will just be set to this number.
	if(num.folds!=length(folds) & !is.null(folds)){
		warning("The num.folds input must match the length of num.folds. Since it doesn't, num.folds has been reset to be equal to length(num.folds).")
		num.folds = length(folds)
	}

	# To ensure proper cross-fitting, folds must be specified along with Q0. 
	if(!is.null(Q0) & is.null(folds)){
		stop("folds must specified if Q0 is specified.")
	}

	# Print a message clarifying to the user what the resource constraint is imposing.
	if(kappa<1){
		message(paste0("The resource constraint imposes that at most a kappa=",kappa," proportion of the population can receive treatment ",txs[1],"."))
	}

	# If baseline.probs is set to NULL, then the mean outcome under the optimal rule is reported (i.e., this value is not contrasted against anything).
	# To achieve this, we set baseline.probs to be a vector of zeros
	if(length(baseline.probs)==0){
		baseline.probs = rep(0,length(txs))
	}

	# reformat SL.library so that it is a list of length-1 or length-2 vectors
	# (where the first entry in a length-2 vector is the learning algorithm,
	#  the second is the screening algorithm)
	SL.library = do.call(c,lapply(SL.library,function(z){
	  if(length(z)>2){
	    lapply(z[2:length(z)],function(zz){c(z[1],zz)})
	  } else {
	    return(list(z))
	  }
	}))

	fits = lapply(1:num.est.rep,function(i){
		if(verbose) message(paste0('Running estimator iteration ',i,'.'))
		if(length(init.ests.in)>0){
			init.ests = init.ests.in[[i]]
		} else {
			init.ests = .sgcvtmle.preprocess(W,A,Y,Delta,SL.library,OR.SL.library,prop.SL.library,missingness.SL.library,txs,family=family,g0=g0,Q0=Q0,num.folds=num.folds,num.SL.rep=num.SL.rep,id=id,folds=folds,obsWeights=obsWeights,stratifyCV=stratifyCV,SL.method=SL.method)
		}
		Q.est = init.ests$Q.est
		blip.est = init.ests$blip.est
		g.est = init.ests$g.est
		lib.blip.est = init.ests$lib.blip.est

		out = .sgtmle(W,A,Y,txs,Q.est,blip.est,g.est,baseline.probs,family=family,kappa=kappa,sig.trunc=sig.trunc,alpha=alpha,obsWeights=obsWeights,RR=RR)
		if(lib.ests){
			out.lib = lapply(lib.blip.est,function(lbe.curr){
				return(.sgtmle(W,A,Y,txs,Q.est,lbe.curr,g.est,baseline.probs,family=family,kappa=kappa,sig.trunc=sig.trunc,alpha=alpha,obsWeights=obsWeights,RR=RR))
			})
			out.list = list(ests=c(out$est,sapply(out.lib,function(x){x$est})),vars=c(var(out$ic),sapply(out.lib,function(x){var(x$ic)})),algs=c('SuperLearner',SL.library))
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
