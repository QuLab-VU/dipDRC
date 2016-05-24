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

plotGC			<-	function(dtp, new.win=TRUE, add.lines=FALSE, get.DIP=FALSE,...)
{
	if(new.win)	dev.new(width=3, height=3)
	par(mar=c(4,3,2,0.5))
	time.name		<-	colnames(dtp)[grep('[Tt]ime',colnames(dtp))]
	date.name		<-	colnames(dtp)[grep('[Dd]ate',colnames(dtp))]
	drug.name		<-	colnames(dtp)[grep('[Dd]rug',colnames(dtp))]
	if(length(drug.name)==0)	drug.name		<-	colnames(dtp)[grep('[Tt]reatment',colnames(dtp))]
	
	cl	<-	unique(dtp[,grep('[Ll]ine',colnames(dtp))])
	u.dates		<-	unique(dtp[,date.name])
	date.col	<-	topo.colors(length(u.dates))

	plot(dtp[,c(time.name,'nl2')], main=paste(cl,unique(dtp[,drug.name])), xlab="", ylab="", xlim=c(0,140), ylim=c(-1,5.5), type='n',...)
	if(get.DIP) DIP <- numeric()
	for(co in unique(dtp$conc)) 
	{
		if(get.DIP) DIP <- append(DIP,findDIP(dtp[dtp$conc==co,c(time.name,'nl2')])$DIP,dat.type='nl2')

		for(uid in unique(dtp[dtp$conc==co,]$uid))
		{
			dfl	<-	dtp[dtp$uid==uid & dtp$conc==co,]
			i	<-	match(unique(dfl$conc), unique(dtp$conc))
			lines(dfl[,time.name], dfl[,'nl2'], col=date.col[match(dfl[,date.name],u.dates)])
			if(add.lines)	curve(x*DIP[i]+int[i], from=72, to=tail(dtp[dtp$conc==co,time.name],1), 
				col='green', lwd=2, add=TRUE)
			text(tail(dfl[,time.name],1), tail(dfl$nl2,1), paste(co,'µM'), cex=0.75, pos=4)
		}
	}
	mtext('Time (h)', side=1, font=2, line=2)
	mtext('Population doublings', side=2, font=2, line=2)
	legend('topleft',legend=u.dates, col=date.col, lwd=1,bty='n',cex=0.75)
	if(get.DIP) 
	{
		out	<-	data.frame(DIP=DIP, cellLine=cl, conc=unique(dtp$conc))
		rownames(out)	<-	NULL
	} else { out <- NA }
	invisible(out)
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
	if(nrow(dtf)<5 | ncol(dtf)<2)
	{
		message('findDIP requires data frame with nrows >= 5 and ncol >= 2')
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


# Function to extract DIP rate across multiple condititons 
# and calculate a 4-param logistic fit
dipDRC	<-	function(dtf, xName='time', yName='cell.count', var=c('cell.line','drug','conc','expt.date'), 
	print.dip=FALSE, norm=FALSE, plotIt=TRUE, toFile=FALSE, showEC50=TRUE, ...)
{	
	if(plotIt & toFile)	pdf('dipDRC_graphs.pdf')
	concName	<-	var[grep('[Cc]onc',var)]
	exptID		<-	var[grepl('[Dd]ate',var) | grepl('[Ii][Dd]',var)][1]
	dtf$u.cond		<-	makeUCond(dtf,var)
	dtf$cell.drug	<-	makeUCond(dtf,var[1:2])
	out	<-	list()
	for(ucd	in unique(dtf$cell.drug))
	{
		temp		<-	dtf[dtf$cell.drug==ucd,]
		Uconc		<-	unique(temp[,concName])
		dip.rates	<-	temp[match(unique(temp$u.cond),temp$u.cond),var]
		rownames(dip.rates)	<-	NULL
		dip.rates	<-	cbind(dip=NA,dip.rates)
		for(r in unique(temp[,exptID]))
		{
			for(co in unique(Uconc)) 
			{
				dip.rates[dip.rates[,concName]==co & dip.rates[,exptID]==r,'dip']	<-	
					findDIP(sumRep(	temp[temp[,concName]==co & temp[,exptID]==r,yName],
									temp[temp[,concName]==co & temp[,exptID]==r,xName]), print.dip=print.dip)$dip
			}
		}
		dip.rates$norm.dip	<-	dip.rates$dip/dip.rates[dip.rates[,concName]==min(dip.rates[,concName]),'dip']
		if(norm)
		{	
			# need to make formula using correct names
			out[[ucd]] <- tryCatch({drm(norm.dip~conc,data=dip.rates,fct=LL.4())},error=function(cond) {return(NA)})
		} else
		{
			out[[ucd]] <- tryCatch({drm(dip~conc,data=dip.rates,fct=LL.4())	},error=function(cond) {return(NA)})
		}
		if(plotIt & !is.na(out[[ucd]][1]))
		{
			plot(out[[ucd]],main=ucd, ...)
			if(showEC50) abline(v=ED(out[[ucd]],50,interval='delta',display=FALSE)[1],col='red')
			abline(h=0, col=grey(0.5), lty=2)
		}
	}
	if(plotIt & toFile) dev.off()
	invisible(out)
}


