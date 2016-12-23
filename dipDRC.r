getPackageIfNeeded <- function(pkg) {
  if (!require(pkg, character.only=TRUE))
    install.packages(pkgs=pkg, dependencies=TRUE)
}

pkgs	<-	c("drc", "car")

sapply(pkgs,getPackageIfNeeded)


require(drc)			# for dose-response curves; version 2.5-12
require(car)			# for linear detailed regression analysis

# fifth order polynomial
p5	<-	function(x,int,b1,b2,b3,b4,b5) int+b1*x+b2*x^2+b3*x^3+b4*x^4+b5*x^5

deriv_poly_coef<-function(co) {
    stopifnot(
    	all(sapply(names(co), FUN=function(x) x %in% c('int','b1','b2','b3','b4','b5')))
    )
    new.co	<-	c(co['b1'],2*co['b2'],3*co['b3'],4*co['b4'],5*co['b5'],0)
    new.co[is.na(new.co)] <- 0
    names(new.co)	<-	c('int','b1','b2','b3','b4','b5')
    new.co
}

# obtaining best fit parameters of a 5th order polynomial
# linear model to 5th order polynomial using p5 function defined above
fit_p5	<-	function(dtf)
{
	x		<-	colnames(dtf)[1]
	y		<-	colnames(dtf)[2]
	
	form	<-	formula(paste(y, '~ p5(', x, ', int, b1, b2, b3, b4, b5)'))
	
	m.p5	<-	nls(form, data=dtf, start=list(int=1,b1=1,b2=1,b3=1,b4=1,b5=1))
	m.p5.coef.1stderiv	<-	deriv_poly_coef(coef(m.p5))
	m.p5.coef.2ndderiv	<-	deriv_poly_coef(m.p5.coef.1stderiv)
	
	list(m=m.p5, coef.1stderiv=m.p5.coef.1stderiv, coef.2ndderiv=m.p5.coef.2ndderiv)
}

log2norm <- function(count, ids, norm_type=c('idx','ref')[1], 
	norm_id="0", norm_vals, norm_idx=1, zero=log2(0.999))
{
# this function will normalize in one of two distinct ways:
# 1) using the position of each vector specified by 'idx'
# 2) using the position identified from matching 'norm_id'
# to its position of the vector from 'norm_val'
# e.g. norm_id = 3 could refer to a vector 'day' containing a time variable
# made this modification on 2014-09-23, but should still work with older code
    
    # l2
    l2 <- log2(count)
		# finds time points with no cells (0), and replace it
		# with zero (log2(0.999)), so that the data can be displayed 
		# in log scale, yet easily found.
    l2[is.infinite(l2)] <- zero
	norm <- numeric()
	group <- as.character(unique(ids))

    if(norm_type=='idx')
    {
		for(i in group)
		{
			d <- l2[ids == i]
			norm <- append(norm, d - d[norm_idx])
		}
	
		norm
		
	}	else	{
	
		for(i in group)
		{
			norm_pos		<-	match(norm_id, as.character(norm_vals)[ids == i])
			d <- l2[ids == i]
			norm <- append(norm, d - d[norm_pos])
		}
		norm
	}

}

# root-mean-square deviation (RMSD) of an estimator == root-mean-square error (RMSE) 
# when the estimator is unbiased, RMSD is the square root of variance == standard error
# RMSD represents the sample standard deviation of the differences between predicted values and observed values.
# resid == vector or matrix of residuals (difference between the observed and predicted values)
rmsd	<-	function(resid)
{
	sqrt(sum((resid)^2)/length(resid))
}



makeUCond	<-	function(dat,var)
{
	if(class(dat) != 'data.frame' & class(var) != 'character')
	{
		message('data or variables sent to makeUCond invalid')
		return(NA)
	}
	if(length(dat[,var])==0)
	{
		message(paste('Could not find',var,'in data frame'))
		return(NA)
	}
	apply(dat[,var],1, function(x) paste(x, collapse='_'))
}



plotGC_DIPfit	<-	function(dtp, tit='unknown', toFile=FALSE, newDev=TRUE, add.line.met='none',...)
{
	stuff	<-	list(...)
	if('o' %in% names(stuff) & tit=='rmse')	tit <- paste(tit,'with o =',stuff[['o']])
	dip	<-	findDIP(dtp,...)
	add.line	<-	FALSE
	if(add.line.met != 'none')	
	{	
		add.line	<-	TRUE
		dip2 <- findDIP(dtp,met=add.line.met)
	}

	fn	<-	paste(tit,nrow(dtp),'points.pdf')
	if('metric' %in% names(stuff))	tit <- paste(tit,stuff[['metric']])
	
	if(newDev & !toFile)	dev.new(width=7.5, height=3)
	if(newDev &toFile)	pdf(file=fn, width=7.5, height=3)
	if(newDev) par(mfrow=c(1,3), oma=c(0,0,1,0))

	plot(dtp, main=NA, xlab=NA, ylab=NA)
	mtext(side=1, 'Time (h)', font=2, line=2)
	mtext(side=2, 'log2(cell number)', font=2, line=2)
	dip.val	<-	round(dip$dip,4)
	dip.95conf	<-	round(abs(dip.val-confint(dip$best.model)[2,1]),5)
	legend("bottomright", c(paste('DIP =',dip.val),paste0('  ±',dip.95conf),paste0('start =',round(dip$start.time,1))), bty='n', pch="")
	curve(coef(dip$best.model)[1]+coef(dip$best.model)[2]*x,from=0,to=150,add=TRUE, col='red', lwd=3)
	
	try(polygon(	x=c(dip$start.time, 150, 150),
		y=c(coef(dip$best.model)[1]+coef(dip$best.model)[2]*dip$start.time,
		coef(dip$best.model)[1]+confint(dip$best.model)[2,1]*150,
		coef(dip$best.model)[1]+confint(dip$best.model)[2,2]*150), col=adjustcolor("gray",alpha.f=0.4), density=NA))
	
	abline(v=dip$start.time, lty=2)
	if(add.line)	abline(v=dip2$start.time, lty=3, col='red')

	# plotting graph of adjusted R2
	plot(dip$eval.times,dip$ar2, ylim=c(0,1), xlab=NA, ylab=NA, main=NA)
	mtext(side=1, 'Time (h)', font=2, line=2)
	mtext(side=2, 'adj R2', font=2, line=2)
	abline(v=dip$start.time, lty=2)
	if(add.line)	abline(v=dip2$start.time, lty=3, col='red')
	
	# plotting graph of RMSE
	plot(dip$eval.times,dip$rmse, ylim=c(0,0.5), xlab=NA, ylab=NA, main=NA)
	mtext(side=1, 'Time (h)', font=2, line=2)
	mtext(side=2, 'RMSE', font=2, line=2)
	abline(v=dip$start.time, lty=2)
	if(add.line)	abline(v=dip2$start.time, lty=3, col='red')
	r1d	<-	coef(fit_p5(data.frame(dip$eval.times,dip$rmse))$m)
	curve(p5(x,r1d[1],r1d[2],r1d[3],r1d[4],r1d[5],r1d[6]), from=0, to=120, col='blue', lwd=3, add=TRUE)

	if(newDev)	mtext(tit, outer=TRUE, side=3, font=2, line=-1.5)
	if(newDev & toFile)	dev.off()
	invisible(dip)
}

findDIP	<-	function(dtf,name='unknown',all.models=FALSE, metric=c('opt','ar2','rmse')[1], o=0, dat.type='cell.counts', print.dip=FALSE)
{
	if(nrow(dtf)<4 | ncol(dtf)<2)
	{
		message('findDIP requires data frame with nrows >= 4 and ncol >= 2')
		return(NA)
	}

	x		<-	dtf[,1]
	if(dat.type=='cell.counts')	{y <- log2(dtf[,2])} else {y <- dtf[,2]}
	n		<-	nrow(dtf)-2
	
	m		<-	list()
	rmse	<-	numeric()
	ar2		<-	numeric()
	p		<-	numeric()
	
	b		<-	y
	a		<-	x
	i		<-	1

	for(i in seq(n))
	{	
	# fit a linear model
		m[[i]]	<-	lm(b ~ a)
	# root-mean-squared-error of the residuals of the data to the best-fit model
		rmse	<-	append(rmse,rmsd(m[[i]]$residuals))
	# adjusted R^2 value of the best-fit model
		ar2		<-	append(ar2,summary(m[[i]])$adj.r.squared)
	# p-value
		p		<-	append(p,summary(m[[i]])$coefficients[2,4])
	# remove the first element from each vector
		a		<-	a[-1]
		b		<-	b[-1]
	}
	
	eval.times	<-	x[seq(n)]
	# fit a 5th order polynomial to the values of rmse for linear models starting at each time point
	# simply used as an way to describe how the values of rmse change as starting time points are dropped 
	rmse.p5		<-	tryCatch({fit_p5(data.frame(x=eval.times,rmse=rmse))},error=function(cond){NA})
	rmse.p5.coef	<-	tryCatch({coef(rmse.p5$m)},error=function(cond){NA})
	# first derivative of the best-fit 5th order polynomial is used to estimate when
	# the change of rmse over time (i.e. the first derivative value at a given time point) approaches zero
	rmse.p5.1std.coef	<-	tryCatch({rmse.p5$coef.1stderiv},error=function(cond){NA})
	
	f	<-	function(...,offset) p5(...) + offset	


	# opt is a combined metric for choosing the best linear model using
	# adjusted R2 value, rmse of the residuals and fraction of the total data points used in model
	# this weightings of each component of this metric were optimized empirically!
	opt	<-	ar2 * (1-rmse)^2 * (length(eval.times)-seq(eval.times))^0.25
	
	# idx is the index of the time after which all data points are included
	# for the best fit linear model
	# three different metrics are provided for identifying this time, adjusted R2 (ar2),
	# root mean square error of the residuals (rmse) or the opt metric described above
	idx	<-	switch(metric,
		# index using best adjusted R2 value within the first 90% of data
		ar2 = ifelse(n<=5,match(max(ar2),ar2),match(max(ar2[seq(floor(n*.9))]),ar2)),			
		
		# index using time at which 1st derivative of rmse becomes 0
		# offset value used if curve does not cross 0
		rmse = tryCatch(												
			{floor(uniroot(f, interval=c(0,72), tol=1e-6, extendInt='upX',
				offset= o,						
				int=rmse.p5.1std.coef[1],
				b1=rmse.p5.1std.coef[2],
				b2=rmse.p5.1std.coef[3],
				b3=rmse.p5.1std.coef[4],
				b4=rmse.p5.1std.coef[5],
				b5=rmse.p5.1std.coef[6]
			)$root)},
			error=function(cond) {
				message("Could not find root of 1st deriv of RMSE")
				message(cond)
            	1
            }
        ),
		opt = match(max(opt),opt)
	)
		

	dip	<-	coef(m[[idx]])[2]
	ci	<-	diff(confint(m[[idx]])[2,])/2
	names(ci)	<-	'±'
	if(print.dip) print(paste('DIP =',round(dip,4),'starting at',round(x[idx],2),'with',n+2,'data points'))
	out=list(data=data.frame(Time_h=x,l2=y), model=m, metric.used=metric, n=n, idx=idx, best.model=m[[idx]], opt=opt, 
		eval.times=x[seq(n)], rmse=rmse,ar2=ar2,p=p,start.time=x[idx],dip=dip, dip.95ci=ci, rmse.p5.1std.coef=rmse.p5.1std.coef)
	if(!all.models)
	{
		out[["model"]] <- NULL
	}
	invisible(out)
}

sumRep	<-	function(count,ids)
{
	uid		<-	unique(ids)
	sums	<-	sapply(unique(ids), FUN=function(x) sum(count[ids==x]))
	cbind(uid,sums)
}


dipDRC <- function(dtf, xName='time', yName='cell.count', var=c('cell.line','drug1','drug1.conc','expt.date'), 
	print.dip=FALSE, norm=FALSE, plotIt=TRUE, toFile=FALSE, fct=LL.4(), uidName='uid', ...)
{	
	# Function to extract DIP rate from single cell line and single drug + control 
	# and calculates a 4-param log-logistic fit by default
	if(plotIt & toFile)	pdf('dipDRC_graph.pdf')
	concName	<-	var[grep('[Cc]onc',var)]
	exptID		<-	var[grepl('[Dd]ate',var) | grepl('[Ii][Dd]',var)][1]
	Uconc		<-	sort(unique(dtf[,concName]),decreasing=TRUE)
	dip.rates	<-	dtf[,c(var,uidName)]
	dip.rates	<-	dip.rates[!duplicated(dip.rates[,uidName]),]
	rownames(dip.rates)	<-	NULL
	dip.rates	<-	cbind(dip=NA,dip.rates)
	for(id in unique(dtf[,uidName]))
	{
		dip.rates[dip.rates[,uidName]==id,'dip']	<-	
			tryCatch({findDIP(	dtf[dtf[,uidName]==id,c(xName,yName)], print.dip=print.dip)$dip
				},
				error=function(cond) {
					message(paste("Error in dipDRC/findDIP:",uidName,'=',id))
					message(cond)
					return(NA)
				}				
			)   
	}
	
	dip.rates$norm.dip	<-	dip.rates$dip/mean(dip.rates[dip.rates[,concName]==0,'dip'])
	if(norm)
	{	
		f <- formula(paste0('norm.dip ~ ',concName))
		out <- tryCatch({drm(f,data=dip.rates,fct=fct)},error=function(cond) {return(dip.rates)})
	} else
	{
		f <- formula(paste0('dip ~ ',concName))
		out <- tryCatch({drm(f,data=dip.rates,fct=fct)	},error=function(cond) {return(dip.rates)})
	}
	if(plotIt & class(out)=='drc')
	{
		plot.dipDRC(out, ...)
		abline(h=0, col=grey(0.5), lty=2)
	}
	if(plotIt & class(out)!='drc')
	{
		temp <- out
		temp[temp$drug1.conc==0,'drug1.conc'] <- min(temp[temp$drug1.conc!=0,'drug1.conc'])/10
		myargs <- list(...)
		myargs['type'] <- NULL
		myargs <- c(list(formula("dip ~ log10(drug1.conc)"), data=temp),myargs)
		do.call(plot,myargs)
		abline(h=0, col=grey(0.5), lty=2)
		legend("bottomleft",'No DRC fit',bty='n')
	}

	if(plotIt & toFile) dev.off()
	invisible(out)
}



plot.dipDRC	<-	function(drmo, plot.type='confidence', showEC50=TRUE, ...)
{
	if(is.na(drmo[1])) 
	{
		message('DRM object not available to plot')
		return(NA)
	}
	param	<-	list(...)
	
	cell.line <- unique(drmo$origData$cell.line)[1]
	drug <- unique(drmo$origData$drug1)[1]
	
	if(!('main' %in% names(param)))
	{
		tit <- paste(cell.line,drug,sep='_')
		param[['main']] <- tit
	}
	if(!('type' %in% names(param)))		param[['type']] <- plot.type

	myargs <- append(list(x=drmo),param)
	do.call(plot,args=myargs)
	if(showEC50) abline(v=ED(drmo,50,interval='delta',display=FALSE)[1],col='red')
	abline(h=0, col=grey(0.5), lty=2)
}

myll4	<-	function(x,b,c,d,e)
{
#                           d - c           
# f(x)  =  c  +  ---------------------------
#                1 + exp(b(log(x) - log(e)))
#	b: Hill coefficient
#	c: lower limit (Emax)
#	d: upper limit (Emin)
#	e: EC50 (not log scaled)
	c + ( (d - c) / (1 + exp(b*(log(x) - log(e)))))
}

addLL4curve <- function(drmodel, fromval=1e-12, toval=1e-5, ...)
{
	param <- coef(drmodel)
	names(param) <- letters[2:5]

	curve(do.call(myll4,args=append(list(x),as.list(param))), from=fromval, to=toval,add=TRUE,...)
}

getParam <- function(drmod)
{
	if(class(drmod) != 'drc') {message('getParam() requires drm object'); return(invisible(NA))}
	# estimates and confidence intervals of log-logistic model fit parameters
	ci <- as.data.frame(confint(drmod))
	ci <- cbind(est=coef(drmod),ci)
	rownames(ci) <- c('slope','Emax','E0','EC50')
	colnames(ci) <- c('est','lower','upper')
	
	# e_halfmax = DIP rate at EC50
	e_halfmax <- predict(drmod, data.frame(drug1.conc=ci['EC50','est']), interval='confidence')
	names(e_halfmax) <- c('est','lower','upper')
	ic10 <- ED(drmod, ci['E0','est']*0.9, type='absolute', interval='delta', display=FALSE)
	rownames(ic10) <- 'IC10'
	ic50 <- ED(drmod, ci['E0','est']/2, type='absolute', interval='delta', display=FALSE)
	rownames(ic50) <- 'IC50'
	ic100 <- ED(drmod, 0, type='absolute', interval='delta', display=FALSE)
	rownames(ic100) <- 'IC100'
	out <- rbind(ci,e_halfmax,ic10[c(1,3,4)],ic50[c(1,3,4)],ic100[c(1,3,4)])
	rownames(out) <- c(rownames(out)[1:4],'E50','IC10','IC50','IC100')
	out
}

