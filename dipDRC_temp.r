
nEach	<-	function(ids)
{
	out	<-	integer()
	for(i in unique(ids)) out <- append(out,length(ids[ids==i]))
	names(out) <- unique(ids)
	out
}

firstInstPos	<-	function(vec)	match(unique(vec), vec)

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
			message(paste('The first 5 pts of ids =',i,'is not linear within ar2 =',min.ar2))
			next
		} else { 
			dtk.times <- append(dtk.times,x)
			dtk.counts <- append(dtk.counts,floor(2^y))
			dtk.ids <- append(dtk.ids, rep(i,length(x)))
		}
	}
	
	data.frame(times=dtk.times,counts=dtk.counts,ids=dtk.ids)
}

findCtrlDIP <- function(times, counts, ids, 
	col.names=c('Time','Cell.counts','Well'),
	type='mean')
{
# data should be filtered first using filterCtrlData()
# ids = unique identifier for each sample (usually a well from an experiment)
# assumes counts in linear scale (i.e. direct cell counts)

	#	
	# types
	#
	# max: identify the sample(s) with the most data points, sum them and fit a lm
	# sum: use all samples trimmed to the fewest data points, sum them and fit a lm
	# mean: use all samples, fitting lm separately to each and return the average

	dip			<-	numeric()
	dip.temp	<-	numeric()

	if(type=='mean')
	{
		for(i in unique(ids))
		{	
			x <- times[ids==i]
			y <- log2(counts[ids==i])

			dip.temp <- append(dip.temp,coef(lm(y ~ x))[2])
			names(dip.temp)[length(names(dip.temp))]	<-	i
		}
		dip <- mean(dip.temp)
		dip	<-	append(dip, sd(dip.temp)/sqrt(length(dip.temp)))
		names(dip)	<-	c('mean.dip','std.err')
		
		count.t0	<-	floor(mean(sapply(unique(ids), FUN=function(x) head(counts[ids==x],1))))
		dts.times	<-	unique(times)
		dts.counts	<-	c(count.t0,floor(count.t0*2^(dip[1]*dts.times[-1])))
	}
	
	if(type=='max')
	{
		n			<-	nEach(ids)
		n.max		<-	max(n)
		max.ids		<-	names(n)[n==n.max]
		dts.times	<-	times[ids %in% max.ids]
		dts.counts	<-	counts[ids %in% max.ids]
		a			<-	aggregate(dts.counts ~ dts.times, FUN=sum)
		m			<-	lm(log2(dts.counts) ~ dts.times, data=a)
		dip			<-	coef(m)[2]
		dip			<-	append(dip, summary(m)$coefficients[2,2])
		names(dip)	<-	c('max.dip','std.err')
		dts.times	<-	a$dts.times
		dts.counts	<-	a$dts.counts
	}
	
	if(type=='sum')
	{
		n			<-	nEach(ids)
		n.min		<-	min(n)
		dts.times	<-	as.numeric(sapply(unique(ids), FUN=function(x) head(times[ids==x],n.min)))
		dts.counts	<-	as.numeric(sapply(unique(ids), FUN=function(x) head(counts[ids==x],n.min)))
		a			<-	aggregate(dts.counts ~ dts.times, FUN=sum)
		m			<-	lm(log2(dts.counts) ~ dts.times, data=a)
		dip			<-	coef(m)[2]
		dip			<-	append(dip, summary(m)$coefficients[2,2])
		names(dip)	<-	c('sum.dip','std.err')
		dts.times	<-	a$dts.times
		dts.counts	<-	a$dts.counts
	}
	
	out.data	<-	data.frame(dts.times,dts.counts)
	colnames(out.data)	<-	col.names[1:2]
	
	out <- list(out.data,dip)
	names(out)	<- c('data','dip')
	
	out
}

controlQC	<-	function(times, counts, ids, col.names=c('Time','Cell.counts','Well'), 
	plotIt=TRUE, ctrl.type='mean', cell.line.name="", ...)
{
	# expecting data to be from a single cell line
	# will filter data for consistency with exponential growth and return
	# a data.frame with a single set of time points and the estimated cell counts
	# undefined arguments will be passed to plot function
	
	fd	<-	filterCtrlData(times, counts, ids)
	if(plotIt)	plotGC(fd[,1],fd[,2],fd[,3],fd[,3], main=cell.line.name, ...)
	out	<-	findCtrlDIP(fd[,1],fd[,2],fd[,3], col.names=col.names, type=ctrl.type)$data
	if(plotIt)	lines(out[,1],log2norm(out[,2],ids=1), lwd=3)
	invisible(out)
}

findMaxDens <- function(times, counts, ids, min.ar2=0.99)
{
# To determine at what density cell counts are no longer increasing exponentially
# times and counts for which adjusted R2 value is >= min.ar2 argument are returned
# ids = unique identifier for each sample (usually a well from an experiment)
# assumes counts in linear scale (i.e. direct cell counts)
	max.exp.counts	<-	integer()
	max.ids		<-	character()
	
	for(i in unique(ids))
	{	
		ttf <- times[ids==i]
		ctf <- log2(counts[ids==i])
		
		if(length(ttf) < 5)
		{
			message(paste('Need at least 5 data points for ',i))
			next
		}

		for(l in 5:(length(ttf)))
		{
		# include data through the 'l'th time point
			x <- ttf[seq(l)]
			y <- ctf[seq(l)]
		# capture last position where adjusted R^2 value <= min.ar2
			if(summary(lm(y ~ x))$adj.r.squared <= min.ar2) break
		}
		
		if(summary(lm(y ~ x))$adj.r.squared < min.ar2 & l ==5)
		{
			message(paste('The first 5 pts of ids =',i,'is not linear within ar2 =',min.ar2))
			next
		} else { 
			max.exp.counts <- append(max.exp.counts,floor(2^y[l]))
			max.ids <- append(max.ids, i)
		}
	}
	
	data.frame(max.exp.counts=max.exp.counts,ids=max.ids)
}

prepCountData <- function(dat, time.count.id.names=c('Time','Cell.count','Well'), 
	max.cell.count=3000, cell.line.name=unique(dat$Cell.line))
{
	ctrl <- dat[dat$Drug1=='control',]
	dat <- dat[dat$Drug1!='control',]
	cn	<- colnames(dat)
	if(length(cell.line.name)>1)
	{
		message('WARNING: Length of prepCountData cell.line.name > 1, using first value only')
		cell.line.name <- cell.line.name[1]
	}
	dat[dat[time.count.id.names[2]] > max.cell.count,time.count.id.names[2]]	<-	NA
	dat$Drug1.conc <- signif(dat$Drug1.conc,2)

	# add control data
	ctrl <- controlQC(ctrl[,time.count.id.names[1]],
		ctrl[,time.count.id.names[2]],ctrl[,time.count.id.names[3]], col.names=time.count.id.names, plotIt=FALSE)
	
	temp <- tail(dat,nrow(ctrl))
	temp[,time.count.id.names[1:2]]	<-	ctrl
	temp[,!(colnames(temp) %in% c(time.count.id.names,'Cell.line','Date'))] <- NA
	temp$uid <- paste0(cell.line.name,'_control')
	temp$Drug1 <- 'control'
	temp$Drug1.conc <- 0
	temp[,time.count.id.names[3]] <- 'CTRL'
	rbind(dat,temp)
}

plotGC_DIPfit2	<-	function(dtp, tit='unknown', toFile=FALSE, newDev=TRUE, add.line.met='none',...)
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
