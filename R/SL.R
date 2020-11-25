
#' SuperLearner for Estimating the Conditional Average Treatment Effect

sg.SL = function(W,A,Y,Delta=rep(1,length(A)),SL.library,OR.SL.library=SL.library,prop.SL.library=SL.library,missingness.SL.library=SL.library,txs=c(0,1),g0=NULL,family=binomial(),num.SL.folds=10,num.SL.rep=5,SL.method="method.NNLS2",id=NULL,obsWeights=NULL,stratifyCV=FALSE,lib.ests=FALSE,...){
	require(SuperLearner)

	if(any(names(list(...))=="separate.reg")) warning("The separate.reg option is deprecated. All outcome regressions are now stratified based on treatment status, i.e. the deprecated separate.reg always equals TRUE.")

	if(any(names(list(...))=="ipcw")) {
		warning("The ipcw argument is deprecated.")
	}

	if(any(names(list(...))=="project")) {
		warning("The project argument deprecated.")
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

	# Recode missing Y values to 0
	Y[Delta==0] = 0

	if((family$family=='binomial') & stratifyCV) {
		blip.cvControl = list(V=num.SL.folds,validRows=CVFolds(length(Y),id,Y,SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.SL.folds)))
	} else {
		if(sum(!(Y%in%c(0,1)))>0 & stratifyCV){
			warning('Can only stratify on binary outcome. Will stratify on estimation of treatment mechanism only.')
		}
		blip.cvControl = list(V=num.SL.folds)
	}
	n = nrow(W)

	if(length(g0)==0){
		if((length(txs)==2) & all(Delta==1)){ # no need to fit propensity for both treatments if there are only two treatments and no outcomes are missing (probabilities add to 1)
			g1 = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(as.numeric(A==1),W,family=binomial(),SL.library=prop.SL.library,cvControl=list(V=num.SL.folds,stratifyCV=stratifyCV),id=id,obsWeights=obsWeights)$SL.predict[,1]}))/num.SL.rep
			g = cbind(g1,1-g1)
		} else {
			g = do.call(cbind,lapply(txs,function(a){
				Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(as.numeric(A==a),W,family=binomial(),SL.library=prop.SL.library,cvControl=list(V=num.SL.folds,stratifyCV=stratifyCV),id=id,obsWeights=obsWeights)$SL.predict[,1]}))/num.SL.rep
				}))
			if(any(rowSums(g)>1)){ # if the estimated treatment probs sum to more than 1, then normalize
				sm = rowSums(g)
				inds = which(sm>1)
				for(i in 1:ncol(g)){
					g[inds,i] = g[inds,i]/sm[inds]
				}
			}
		}
	} else {
		g=g0
	}

	if(all(Delta==1)){
		g.delta = Delta
	} else {
		g.delta = Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(Delta,data.frame(W,A=A),family=binomial(),SL.library=missingness.SL.library,obsWeights=obsWeights,id=id,cvControl=SuperLearner.CV.control(stratifyCV=stratifyCV,V=num.SL.folds))$SL.predict[,1]}))/num.SL.rep
	}

	# recode A so that it's set to missing (some value not in txs -- max(txs)+1 will work)
	# whenever Delta is set to missing
	A[Delta==0] = max(txs) + 1
	# modify treatment probability to also account for missingness probability
	g = g * do.call(cbind,lapply(1:length(txs),function(i){g.delta}))

	Qbar = do.call(cbind,lapply(txs,function(a){
		Reduce('+',lapply(1:num.SL.rep,function(i){SuperLearner(Y[A==a],W[A==a,],newX=W,family=family,SL.library=OR.SL.library,cvControl=list(V=num.SL.folds,stratifyCV=(family$family=='binomial') & stratifyCV),id=id[A==a],obsWeights=obsWeights[A==a])$SL.predict[,1]}))/num.SL.rep
		}))

	# for each treatment, SL.objs contains a list of SL fits estimating the average treatment
	# effect of assigning that treatment versus assigning the given treatment at random according
	# to the propensity observed in the population
	SL.objs = lapply(seq(txs),function(i){
		a = txs[i]
		Z = (A==a)/((A==a)*g[,i] + (A!=a)) * (Y-Qbar[,i]) + Qbar[,i] - Delta*Y/g.delta

		lapply(1:num.SL.rep,function(i){
			SuperLearner(Z,W,SL.library=SL.library,cvControl=blip.cvControl,id=id,obsWeights=obsWeights,method=SL.method)})
		})
	SL.cate.fun = function(w){do.call(cbind,lapply(seq(txs),function(i){
		out = Reduce('+',lapply(1:num.SL.rep,function(j){predict(SL.objs[[i]][[j]],newdata=w)$pred}))/num.SL.rep
		if(family$family=='binomial') out = pmin(pmax(out,-1),1)
		return(out)}))}
	blip = SL.cate.fun(W)

	if(lib.ests){
		lib.cate.fun = function(w){
			list.out = lapply(seq(SL.library),function(k){
				do.call(cbind,lapply(seq(txs),function(i){
					out = Reduce('+',lapply(1:num.SL.rep,function(j){predict(SL.objs[[i]][[j]],newdata=w)$library.predict[,k]}))/num.SL.rep
					if(family$family=='binomial') out = pmin(pmax(out,-1),1)
					return(out)}))
				})
			names(list.out) = SL.library
			return(list.out) }
		lib.ests = lib.cate.fun(W)
		return(list(est=blip,SL.cate.fun=SL.cate.fun,SL=SL.objs,Qbar=Qbar,lib.ests=lib.ests,lib.cate.fun=lib.cate.fun))
	} else {
		return(list(est=blip,SL.cate.fun=SL.cate.fun,SL=SL.objs,Qbar=Qbar))
	}
}
