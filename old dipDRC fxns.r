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