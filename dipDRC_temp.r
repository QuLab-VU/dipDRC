dipDRC	<-	function(dtf, xName='time', yName='cell.count', var=c('conc','expt.date'), 
	print.dip=FALSE, norm=FALSE, plotIt=TRUE, toFile=FALSE, showEC50=TRUE, ...)
{	
	if(plotIt & toFile)	pdf('dipDRC_graph.pdf')
	concName	<-	var[grep('[Cc]onc',var)]
	exptID		<-	var[grepl('[Dd]ate',var) | grepl('[Ii][Dd]',var)][1]
	out	<-	list()

	Uconc		<-	unique(dtf[,concName])
	dip.rates	<-	dtf[,var]
	rownames(dip.rates)	<-	NULL
	dip.rates	<-	cbind(dip=NA,dip.rates)
	for(r in unique(dtf[,exptID]))
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
		f <- formula(paste0('norm.dip ~ ',concName))
		out[[ucd]] <- tryCatch({drm(f,data=dip.rates,fct=LL.4())},error=function(cond) {return(NA)})
	} else
	{
		f <- formula(paste0('dip ~ ',concName))
		out[[ucd]] <- tryCatch({drm(f,data=dip.rates,fct=LL.4())	},error=function(cond) {return(NA)})
	}
	if(plotIt & !is.na(out[[ucd]][1]))
	{
		plot(out[[ucd]], ...)
		if(showEC50) abline(v=ED(out[[ucd]],50,interval='delta',display=FALSE)[1],col='red')
		abline(h=0, col=grey(0.5), lty=2)
	}

	if(plotIt & toFile) dev.off()
	invisible(out)
}

plot.dipDRC	<-	function(drm, cell.line="Cell line", drug="drug", type='confidence', showEC50=TRUE, ...)
{
	if(is.na(drm[1])) 
	{
		message('DRM object not available to plot')
		return(NA)
	}
	if(!exists(title)
	dtp	<-	drm[[names(drm) == paste(cell.line,drug,sep='_')]]
	plot(drm, main=paste()
	abline(h=0, col=grey(0.5), lty=2)
}

# root-mean-square deviation (RMSD) of an estimator == root-mean-square error (RMSE) 
# when the estimator is unbiased, RMSD is the square root of variance == standard error
# RMSD represents the sample standard deviation of the differences between predicted values and observed values.
# resid == vector or matrix of residuals (difference between the observed and predicted values)
rmsd	<-	function(resid)
{
	sqrt(sum((resid)^2)/length(resid))
}


densMax <- function(times, counts, ids, min.ar2=0.99)
{
# To determine whether cell counts are exponentially increasing throughout the entire time span
# if not, the time at which linearity fails is returned
# assumes all data are from the same condition (i.e. same cell line)
# ids = unique identifier for each sample (usually a well from an experiment)
# assumes counts in linear scale (i.e. direct cell counts)

	m		<-	numeric()
	
	for(i in unique(ids))
	{	
		ttf <- times[ids==i]
		ctf <- log2(counts[ids==i])
		
		if(length(ttf) < 5)
		{
			message(paste('Need at least 5 data points for ',i))
			next
		}

		# m.ad = model with all time points
		m.ad	<-	lm(ctf ~ ttf)
		ar2	<-	summary(m.ad)$adj.r.squared

		for(l in (length(ttf)):5)
		{
		# include data through the 'l'th time point
			x <- ttf[seq(l)]
			y <- ctf[seq(l)]
		# capture position of max counts when adjusted R^2 value >= min.ar2
			if(summary(lm(y ~ x))$adj.r.squared >= min.ar2) break
		}
#		message(paste(i,l,2^ctf[l]))
		max.dens <- 2^ctf[l]

		m <- append(m,max.dens)
	}
	return(m)
}

findCtrlDIP <- function(times, counts, ids, min.ar2=0.99)
{
# To determine whether cell counts are exponentially increasing throughout the entire time span
# times and counts for which adjusted R2 value is >= min.ar2 argument are returned
# ids = unique identifier for each sample (usually a well from an experiment)
# assumes counts in linear scale (i.e. direct cell counts)
	m			<-	numeric()
	dip			<-	numeric()
	
	for(i in unique(ids))
	{	
		ttf <- times[ids==i]
		ctf <- log2(counts[ids==i])
		
		if(length(ttf) < 5)
		{
			message(paste('Need at least 5 data points for ',i))
			next
		}

		# m.ad = model with all time points
		m.ad	<-	lm(ctf ~ ttf)
		ar2	<-	summary(m.ad)$adj.r.squared

		for(l in (length(ttf)):5)
		{
		# include data through the 'l'th time point
			x <- ttf[seq(l)]
			y <- ctf[seq(l)]
		# capture last position where adjusted R^2 value >= min.ar2
			if(summary(lm(y ~ x))$adj.r.squared >= min.ar2) break
		}
		
		if(summary(lm(y ~ x))$adj.r.squared < min.ar2)
		{ next } else { dip <- append(dip,coef(lm(y ~ x))[2]) }
	}
	dip
}

filterCtrlData <- function(times, counts, ids, min.ar2=0.99)
{
# To determine whether cell counts are exponentially increasing throughout the entire time span
# times and counts for which adjusted R2 value is >= min.ar2 argument are returned
# ids = unique identifier for each sample (usually a well from an experiment)
# assumes counts in linear scale (i.e. direct cell counts)
	dip			<-	numeric()
	dtk.times	<-	integer()
	dtk.counts	<-	integer()
	dtk.ids		<-	character()
	for(i in unique(ids)) x <- append(x,length(ids[ids==i]))
	min.pts		<-	max(x)
	
	for(i in unique(ids))
	{	
		ttf <- times[ids==i]
		ctf <- log2(counts[ids==i])
		
		if(length(ttf) < 5)
		{
			message(paste('Need at least 5 data points for ',i))
			next
		}

		# m.ad = model with all time points
		m.ad	<-	lm(ctf ~ ttf)
		ar2	<-	summary(m.ad)$adj.r.squared

		for(l in (length(ttf)):5)
		{
		# include data through the 'l'th time point
			x <- ttf[seq(l)]
			y <- ctf[seq(l)]
		# capture last position where adjusted R^2 value >= min.ar2
			if(summary(lm(y ~ x))$adj.r.squared >= min.ar2) break
		}
		
		if(summary(lm(y ~ x))$adj.r.squared < min.ar2)
		{
			message(paste('ids =',i,'is not linear within ar2 =',min.ar2))
			next
		} else { 
			
			dtk.times <- append(dtk.times,x)
			dtk.counts <- append(dtk.counts,floor(2^y))
			dtk.ids <- append(dtk.ids, rep(i,length(x)))
		}
	}
	
	data.frame(times=dtk.times,counts=dtk.counts,ids=dtk.ids)
}
