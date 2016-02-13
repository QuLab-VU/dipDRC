# DIP rate paper functions and all required libraries
# 

require(deSolve)		# ODE solver usedfor generating outputs from 2-state model
						# without partial equilibrium assumption
require(gplots)			# needed for colors
require(grid)
require(drc)			# for dose-response curves; version 2.5-12
library(car)			# for linear detailed regression analysis


twoStateModel	<-	function(times, states, parameters)
{
	with(as.list(c(states, parameters)),{
		# rates
		dCell		<-	(DIP0 - k_on * Drug) * Cell + (k_off * CellPrime)
		dCellPrime	<-	(DIPmax - k_off) * CellPrime + (k_on * Drug * Cell)
		
		# return list of rates
		list(c(dCell,dCellPrime))
	})
}

dipPEA		<-	function(k_on, k_off,  DIP0, DIPmax, drug)
# analytical solution for value of DIP rate under PEA
{
	( 1 / log(2) ) * ( k_off * DIP0 + (k_on * drug * DIPmax) ) / (k_on * drug + k_off )
}

biasF	<- function(kx,ky,t)	exp((kx-ky)*t)


#
# 4-parameter log-logistic function
#
#                           d - c           
# f(x)  =  c  +  ---------------------------
#                1 + exp(b(log(x) - log(e)))
#
#	b: Hill coefficient
#	c: lower limit (Emax)
#	d: upper limit (Emin)
#	e: EC50
#

getAA	<-	function(drmod, conc=drug.conc, max=1)
{
	# drmod is dose–response model (drm)
	# normalize to max of 0 by subtracting 1 from all values
	# and divide by number of measurements
	vals	<-	PR(drmod, conc)
	-sum(vals-1) / length(vals)
}

getGI	<-	function(t0val,tnval,cnval)
# t0val = starting cell #
# tnval = cell # at time n (in specific condition)
# cnval = control cell # at time n
{
	(tnval-t0val) / (cnval-t0val)
}

getGInci	<-	function(t0val,tnval,c0val=t0val,cnval)
# t0val = starting cell #
# tnval = cell # at time n (in specific condition)
# c0val = starting cell # of control
# cnval = control cell # at time n
{
	ifelse(	tnval > t0val, 
			(tnval-t0val) / (cnval-c0val),
			(tnval-t0val) / c0val
	)
}


getCellCount	<-	function(parameters)
{
	for(dc in drug.conc)
	{	
		p	<-	append(parameters,dc)
		names(p)[5]	<-	'Drug'
		temp		<- as.data.frame(ode(y = states, times = times, func = twoStateModel, parms = p))
		Cell.tot	<-	as.numeric(temp$Cell+temp$CellPrime)
		ifelse(dc==drug.conc[1], out <- data.frame(times,Cell.tot), out <- cbind(out,Cell.tot))
	}
	colnames(out)	<-	c('time',as.character(drug.conc))
	out
}


getDRM	<-	function(d=PC9)
{
	allDRM	<-	list()
	for(cl in unique(d$cellLine))	
	{
		m				<-	drm(resp.ratio ~ conc, data=d[d$cellLine==cl,], fct=LL.4())		# to see curves from each expt add: , ECID
		name			<-	ifelse(cl=='par','PC9',cl)
		allDRM[[name]]	<-	m
	}
	allDRM
}

assembleDIP	<-	function(dat)
{
	mean.rates	<-	data.frame(subline=sort(unique(dat$Subline)),DIP=0,DIP.dev=0)
	all.rates	<-	data.frame()
	for(sl in sort(unique(dat$Subline)))
	{
		dtf	<-	dat[dat$Subline==sl,]
# uncomment if you want to see the raw data to which the model is being fit
#		plot(log2(Cell.count) ~ Time_h, data=dtf[dtf$Time_h>72,], main=sl)
		m	<-	lm(log2(Cell.count) ~ Time_h * factor(UID), data=dtf[dtf$Time_h>72,])
		dip	<-	coef(m)[grep('Time_h',names(coef(m)))]
		dip	<-	c(dip[1],dip[-1]+dip[1])
		names(dip)	<-	gsub('Time_h:factor\\(UID\\)','',names(dip))
		names(dip)[1]	<-	setdiff(unique(dtf$UID),names(dip))
		all.rates	<-	rbind(all.rates,data.frame(Subline=sl,DIP=dip, UID=names(dip)))
		mean.val	<-	t.test(dip)[['estimate']]
		dev.val		<-	( t.test(dip)['conf.int'][[1]][2]-t.test(dip)['conf.int'][[1]][1] ) / 2
		mean.rates[mean.rates$subline==sl,2:3]	<-	c(mean.val,dev.val)
	}
	list(all.rates=all.rates, mean.rates=mean.rates)
}

generateStats	<-	function(dat=DS.corr)
{
	stats	<-	aggregate(val ~ subl + var, data=dat, mean)
	colnames(stats)[3]	<-	'mean'
	stats	<-	cbind(stats,aggregate(val ~ subl + var, data=dat, length)['val'])
	colnames(stats)[4]	<-	'n'

	temp	<-	as.data.frame(aggregate(val ~ subl + var, data=dat, FUN=function(x) as.numeric(t.test(x)$conf.int))[,3])
	dev		<-	(temp[,2]-temp[,1])/2
	stats	<-	cbind(stats,dev)

	# add EC50 values
	ec50.mean	<-	sapply(allDRM,FUN=function(x) ED(x,50,interval = "delta", display=FALSE)[1])
	ec50.dev	<-	(sapply(allDRM,FUN=function(x) ED(x,50,interval = "delta", display=FALSE)[3])-
		sapply(allDRM,FUN=function(x) ED(x,50,interval = "delta", display=FALSE)[2]))/2
	ec50.n		<-	sapply(allDRM,FUN=function(x) length(unique(x$origData$exptDate)))
	ec50		<-	cbind(mean=ec50.mean,dev=ec50.dev,n=ec50.n)
	ec50[,1:2]	<-	apply(ec50[,1:2], c(1,2), FUN=function(x) signif(x*1e9,3))		# convert to nM

	stats		<-	rbind(stats,data.frame(subl=rownames(ec50),var='EC50',ec50))
	stats		<-	stats[order(stats$var,stats$subl),]
	stats[,3:4]	<-	apply(stats[,3:4], c(1,2), FUN=function(x) signif(x, digits=3))
	stats		<-	stats[!stats$subl %in% c('DS8','DS9','PC9'),]
	stats$dev	<-	sapply(stats$dev,FUN=function(x) signif(x,3))
	rm(list=c('ec50.mean','ec50.dev','ec50.n'))

	type	<-	as.character(unique(stats$var))
	stats
}

scatterPlotErr	<-	function(dtp,x.var,y.var,add.line=FALSE,xlab.txt=x.var,ylab.txt=y.var,...)
{
	x		<-	dtp[dtp$var==x.var,]
	y		<-	dtp[dtp$var==y.var,]
			
	# plots with y-axis errors
	plotCI(
		x$mean,
		y$mean,
		uiw=y$dev,
		err='y',
		pch=NA,
		...
	)

	# plots with x-axis errors
	plotCI(
		x$mean,
		y$mean,
		uiw=x$dev,
		err='x',
		pch=NA,
		add=TRUE,
		...
	)
	points(x$mean,y$mean,pch=19)
	mtext(side=1, xlab.txt, font=2, line=2)
	mtext(side=2, ylab.txt, font=2, line=2)	
	my.fit<- lm(y$mean ~ x$mean)
	ar2	<-	round(summary(my.fit)$adj.r.squared,2)
	mtext(as.expression(bquote(Adj~R^2~'='~.(ar2))), side=1, line=-1.25, cex=0.75,adj=.95)
	if(add.line)	abline(coef(my.fit),lty=2,lwd=2)	
}

getDIP		<-	function(dtf)	{
	myMod	<-	lm(nl2 ~ time, data=dtf)
	coef(myMod)['time']
}

getDRM	<-	function(dat=PC9)
{
	allDRM	<-	list()
	for(cl in unique(dat$cellLine))	
	{
		m				<-	drm(resp.ratio ~ conc, data=dat[dat$cellLine==cl,], fct=LL.4())		# to see curves from each expt add: , ECID
		name			<-	ifelse(cl=='par','PC9',cl)
		allDRM[[name]]	<-	m
	}
	allDRM
}

# Function to extract DIP rate across multiple condititons 
# and calculate a 4-param logistic fit

drcDIP	<-	function(dtf,cl=cl,drg=drg,PIP.range=PIP.range,norm=T)
{	
	rates	<-	data.frame(cellLine=character(), conc=integer(),
		DIPrate=numeric(),drug=character())
	dtp		<-	dtf[dtf$cellLine==cl & dtf$drug==drg,]
	Uconc	<-	unique(dtp$conc)
	for(co in 2:length(Uconc))
	{
		rates	<-	rbind(rates, 
			data.frame(cellLine=cl,conc=Uconc[co],
			DIPrate=getDIP(dtp[dtp$conc==Uconc[co] & dtp$time>=PIP.range[1] & dtp$time<=PIP.range[2],]),drug=drg))
	}
	rates$normDIP	<-	rates$DIPrate/rates[rates$conc==min(rates$conc),'DIPrate']
	if(norm)
	{	
		drm(normDIP~conc,data=rates,fct=LL.4())
	} else
	{
		drm(DIPrate~conc,data=rates,fct=LL.4())	
	}
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

rmsd	<-	function(resid)
{
# root-mean-square deviation (RMSD) of an estimator == root-mean-square error (RMSE) 
# when the estimator is unbiased, RMSD is the square root of variance == standard error
# RMSD represents the sample standard deviation of the differences between predicted values and observed values.
# resid == vector or matrix of residuals (difference between the observed and predicted values)
	sqrt(sum((resid)^2)/length(resid))
}


findDIP	<-	function(dtf,name='unknown',all.models=FALSE, metric=c('opt','ar2','rmse')[1], o=0)
{
	if(nrow(dtf)<5 | ncol(dtf)<2)
	{
		message('findDIP requires data frame with nrows >= 5 and ncol >= 2')
		return(NA)
	}

	x		<-	dtf[,1]
	y		<-	dtf[,2]
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
	print(paste('DIP =',dip,'starting at',x[idx],'with',n+2,'data points'))
	out=list(data=data.frame(Time_h=x,l2=y), model=m, metric.used=metric, n=n, idx=idx, best.model=m[[idx]], opt=opt, 
		eval.times=x[seq(n)], rmse=rmse,ar2=ar2,p=p,start.time=x[idx],dip=dip, rmse.p5.1std.coef=rmse.p5.1std.coef)
	if(!all.models)
	{
		out[["model"]] <- NULL
	}
	return(out)
}



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

f	<-	function(...,offset) p5(...) + offset	

# plotting growth curve and predicted curve of 5th order polynomial fit
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
	legend("bottomright", c(paste('DIP =',round(dip$dip,4)),paste0('start =',round(dip$start.time,1))), bty='n', pch="")
	curve(coef(dip$best.model)[1]+coef(dip$best.model)[2]*x,from=0,to=150,add=TRUE, col='red', lwd=3)
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


f	<-	function(...,offset=0) p5(...) + offset	

getWellRates	<-	function(raw, time.range=c(70,120))
{
	timeName	<-	colnames(raw)[grep('[Tt]ime', colnames(raw))]
	wellName	<-	colnames(raw)[grep('[Ww]ell',colnames(raw))]
	dateName	<-	colnames(raw)[grep('[Dd]ate',colnames(raw))]
	if(length(wellName)>1)	wellName	<-	wellName[nchar(wellName)==4]
	f	<-	formula(paste('nl2 ~ ',timeName,' * ',wellName))
	m	<-	lm(f, data=raw[raw[,timeName] > time.range[1] & raw[,timeName] < time.range[2],])
	wells	<-	unique(raw[,wellName])
	rates	<-	coef(m)[grep(timeName,names(coef(m)))]
	rates	<-	c(rates[1],rates[-1]+rates[1])
	cl		<-	unique(raw$cellLine)
	expt	<-	ifelse(is.null(unique(raw[,dateName])), 'unknown date',unique(raw[,dateName]))
	out		<-	data.frame(Well=wells, DIP=rates, cellLine=cl, Date=expt)
	rownames(out)	<-	NULL
	out
}
