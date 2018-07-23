# g.est is vector giving estimated probability of A=1 given W

.sgtmle.bothactive = function(W,A,Y,Q.est,blip.est,g.est,family=binomial(),alpha=0.05,trunc.Q=c(0.005,0.995),obsWeights=NULL,RR=FALSE,ipcw=FALSE){
	n = nrow(W)
	if(is.null(obsWeights)) obsWeights = rep(1,length(Y))	
# **** Decide if want to keep this
# Q.est.tmp = Q.est
# Q.est[,1] = (Q.est.tmp[,1] + Q.est.tmp[,2] - blip.est)/2
# Q.est[,2] = (Q.est.tmp[,1] + Q.est.tmp[,2] + blip.est)/2
	if(family$family=='binomial') for(i in 1:2){Q.est[,i] = pmin(pmax(Q.est[,i],trunc.Q[1]),trunc.Q[2])}

	if(RR){
		# remember: for relative risk, contrast probabilities of event Y /not/ occurring (since Y assumed beneficial)
		blip.est.fl = 1-blip.est
		Q.est.fl = 1-Q.est
		Y.fl = 1-Y

		d = as.numeric(blip.est.fl>0) 

		offs.d = family$linkfun(Q.est.fl[,1] + (Q.est.fl[,2]-Q.est.fl[,1]) * d)
		offs.rct = family$linkfun(1/2 * Q.est.fl[,1] + 1/2 * Q.est.fl[,2])

		H.d = (A==d)/(A*g.est + (1-A)*(1-g.est)) 
		H.rct = (1/2) * (1/g.est + 1/(1-g.est))

		eps.d = suppressWarnings(coef(glm(Y.fl ~ -1 + offset(offs.d) + H.d,family=family,weights=obsWeights)))
		eps.rct = suppressWarnings(coef(glm(Y.fl ~ -1 + offset(offs.rct) + H.rct,family=family,weights=obsWeights)))

		Q.d =  family$linkinv(offs.d + eps.d/(d*g.est + (1-d)*(1-g.est)))
		Q.rct =  family$linkinv(offs.rct + eps.rct * H.rct)

		EY.fl.d = mean(obsWeights*Q.d)
		EY.fl.rct = mean(obsWeights*Q.rct)

		est = EY.fl.d/EY.fl.rct
		ic = obsWeights * (1/EY.fl.d * (H.d * (Y.fl-Q.d) + Q.d-EY.fl.d) - 1/EY.fl.rct * (H.rct * (Y.fl-Q.rct) + Q.rct - EY.fl.rct))
		ci = exp(log(est) + c(-1,1) * qnorm(1-alpha/2)*sd(ic)/sqrt(n))
	} else {
		g.star = 2*as.numeric(blip.est>0)-1

		offs = family$linkfun(A*Q.est[,2] + (1-A)*Q.est[,1])
		H = (2*A-1)/(A*g.est + (1-A)*(1-g.est)) * g.star
		eps = suppressWarnings(coef(glm(Y ~ -1 + offset(offs) + H,family=family,weights=obsWeights)))

		Q0 = family$linkinv(family$linkfun(Q.est[,1]) - eps*g.star/(1-g.est))
		# Q1 = family$linkinv(family$linkfun(Q.est[,2]) + eps*g.star/g.est)
		Q1 = family$linkinv(family$linkfun(Q.est[,2]) + eps*g.star/g.est)
		# Divide by two below so get mean under optimal treatment strategy minus mean under 50-50 randomization
		est = mean(obsWeights*(Q1-Q0)*g.star) / 2
		ic = obsWeights * (H*(Y-(A*Q1 + (1-A)*Q0)) + g.star*(Q1-Q0) - 2*est) / 2
		ci = est + c(-1,1) * qnorm(1-alpha/2)*sd(ic)/sqrt(n)
	}

	# ic is the influence curve of the ATE if RR=FALSE, otherwise is the influence curve of the log relative risk
	return(list(est=est,ci=ci,ic=ic))
}