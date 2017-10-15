if(!(exists('findDIP'))) message('WARNING: findDIP not found. plotGC_DIPfit requires dipDRC to be loaded.')

plotGC <- function(x, y, uid, rep=uid, y.type='count', color=TRUE, leg=TRUE, ...)
{
    # x = time
    # y = counts (linear or log scale)
    # uid = unique condition
    # rep = replicate id
    ifelse(y.type=='count',
        plot.y <- log2norm(y,ids=uid),
        plot.y <- y)
    urep <- unique(rep)
    ifelse(color, rep.col <- gplots::colorpanel(n=length(urep),low='blue',mid='orange',high='red'), rep.col <- rep('black',length(urep)))
    plot(x,plot.y,type='n',...)
    for(myrep in urep)
        for(id in unique(uid))
            lines(x[uid==id & rep==myrep],plot.y[uid==id & rep==myrep], col=rep.col[match(myrep,urep)])
    if(leg) legend('topleft',legend=urep,col=rep.col,lwd=1)
}

plotGC_DIPfit <- function(dtp, tit='unknown', toFile=FALSE, newDev=TRUE, add.line.met='none',...)
{
    stuff <- list(...)
    if('o' %in% names(stuff) & tit=='rmse')    tit <- paste(tit,'with o =',stuff[['o']])
    dip <- findDIP(dtp,...)
    add.line <- FALSE
    if(add.line.met != 'none')    
    {    
        add.line <- TRUE
        dip2 <- findDIP(dtp,met=add.line.met)
    }

    fn <- paste(tit,nrow(dtp),'points.pdf')
    if('metric' %in% names(stuff))    tit <- paste(tit,stuff[['metric']])
    
    if(newDev & !toFile)    dev.new(width=7.5, height=3)
    if(newDev &toFile)    pdf(file=fn, width=7.5, height=3)
    if(newDev) par(mfrow=c(1,3), oma=c(0,0,1,0))

    plot(dip$data, main=NA, xlab=NA, ylab=NA)
    mtext(side=1, 'Time (h)', font=2, line=2)
    mtext(side=2, 'log2(cell number)', font=2, line=2)
    dip.val <- round(dip$dip,4)
    dip.95conf <- round(abs(dip.val-confint(dip$best.model)[2,1]),5)
    legend("bottomright", c(paste('DIP =',dip.val),paste0('  ±',dip.95conf),paste0('start =',round(dip$start.time,1))), bty='n', pch="")
    curve(coef(dip$best.model)[1]+coef(dip$best.model)[2]*x,from=0,to=150,add=TRUE, col='red', lwd=3)
    
    try(polygon(    x=c(dip$start.time, 150, 150),
        y=c(coef(dip$best.model)[1]+coef(dip$best.model)[2]*dip$start.time,
        coef(dip$best.model)[1]+confint(dip$best.model)[2,1]*150,
        coef(dip$best.model)[1]+confint(dip$best.model)[2,2]*150), col=adjustcolor("gray",alpha.f=0.4), density=NA))
    
    abline(v=dip$start.time, lty=2)
    if(add.line)    abline(v=dip2$start.time, lty=3, col='red')

    # plotting graph of adjusted R2
    plot(dip$eval.times,dip$ar2, ylim=c(0,1), xlab=NA, ylab=NA, main=NA)
    mtext(side=1, 'Time (h)', font=2, line=2)
    mtext(side=2, 'adj R2', font=2, line=2)
    abline(v=dip$start.time, lty=2)
    if(add.line)    abline(v=dip2$start.time, lty=3, col='red')
    
    # plotting graph of RMSE
    plot(dip$eval.times,dip$rmse, ylim=c(0,0.5), xlab=NA, ylab=NA, main=NA)
    mtext(side=1, 'Time (h)', font=2, line=2)
    mtext(side=2, 'RMSE', font=2, line=2)
    abline(v=dip$start.time, lty=2)
    if(add.line)    abline(v=dip2$start.time, lty=3, col='red')
    r1d <- coef(fit_p5(data.frame(as.numeric(dip$eval.times),dip$rmse))$m)
    curve(p5(x,r1d[1],r1d[2],r1d[3],r1d[4],r1d[5],r1d[6]), from=0, to=120, col='blue', lwd=3, add=TRUE)

    if(newDev)    mtext(tit, outer=TRUE, side=3, font=2, line=-1.5)
    if(newDev & toFile)    dev.off()
    invisible(dip)
}

plotCtrlGC <- function(dat=d, type=c('raw','lin')[1], uniq='plate.name', ...)
{
    prepMultiGCplot()
    for(u in unique(dat[,uniq]))
    {
        dtp <- d[d[,uniq]==u & d$drug1=='control',]
        dtp <- dtp[order(dtp$well,dtp$time),]
        switch(type,
            lin = controlQC(dtp$time, dtp$cell.count, dtp$well,
                c('time','cell.count','well'), cell.line.name=pn, ...),
            raw = plotGC(dtp$time, dtp$cell.count, dtp$uid, dtp$uid, ...)
        )
    }
}

plotRaw <- function(dat=ad, id.name='cell.line', toFile=FALSE, fn='RawGC', ...) 
{
    # id.name = parameter to subset data by
    for(id in unique(dat[,id.name]))
    {
        dtp <- dat[dat[,id.name]==id,]
        dtp <- dtp[order(dtp$drug1.conc,dtp$time),]
        prepMultiGCplot(toFile=toFile,fn=paste(fn,'_',id,'.pdf',sep=''))
        plotAllDrugs(dtp, ylim=c(-2,7), ...)
        if(toFile) dev.off()
    }
}

plotCellLineDrug <- function(dat=ad, expt.date='20160513',
    toFile=FALSE,max.count=NA, max.count.type='by.cell.line',plotFit=FALSE, verbose=FALSE)
{
    out <- data.frame()
    if(is.na(max.count) & max.count.type!='by.cell.line')
        max.count <- min(findMaxCounts(dat$time,dat$cell.count,paste(dat$cell.line,dat$well,sep='_'),verbose=verbose))

    for(cl in unique(dat$cell.line))
    {
        if(is.na(max.count) & max.count.type == 'by.cell.line')
        {
            dfm <- dat[dat$Cell.line==cl & dat$drug1=='control',]
            max.count <- min(findMaxCounts(dfm$time,dfm$cell.count,dfm$well,verbose=verbose))
        }
        d <- prepCountData(dat[dat$cell.line==cl,], max.cell.count=max.count)
        
        fn <- paste(expt.date,cl,'growth curves.pdf')
        if(!toFile)    dev.new(width=11,height=8.5)
        if(toFile)    pdf(file=fn, width=11,height=8.5)
        par(mfrow=c(3,5), mar=c(2,3,1.5,.5))

        for(dr in unique(d$drug1))
        {
            if(dr=='control') next
            dtp <- d[d$drug1==dr | d$drug1=='control',]
            if(plotFit)
            {
	# should modify plotGC_DIPfit function to allow summing replicate wells and add
	# lines of estimated DIP rate for each concentration of drug
	# currently no line fit can be plotted (no plotGC_DIPfit2 function)
                try(plotGC_DIPfit2(
                    data.frame(x=dtp$time,
                    y=dtp$cell.count,
                    uid=dtp$well,
                    rep=dtp$drug1.conc),
                    main=paste(cl,dr),
                    ylim=c(-2,6), xlab=NA,ylab=NA))
            } else {
                plotGC(
                    x=dtp$time,
                    y=dtp$cell.count,
                    uid=dtp$well,
                    rep=dtp$drug1.conc,
                    main=paste(cl,dr),
                    ylim=c(-2,6), xlab=NA,ylab=NA)
            }
        out <- rbind(out,d)
        }
        if(toFile)    dev.off()
    }
    invisible(out)
}

plotAllDrugs <- function(dat, drug.col='drug1', uid='well', output=FALSE, ...)
{
    # dat expected output from prepCountData 
    cl <- unique(dat$cell.line)
    if(length(cl)>1)
    {
        message('prepDRCplot needs data.frame containing single cell line')
        return(invisible(NA))
    }

    out <- list()
    drugs <- sort(unique(dat[,drug.col]))
    drugs <- setdiff(drugs,"control")
    for(dr in drugs)
    {
        df.name <- paste(cl,dr,sep='_')
        dtp <- dat[dat[,drug.col]==dr | dat[,drug.col]=='control',]
        dtp <- dtp[order(dtp$drug1.conc,dtp[,uid],dtp$time),]
        plotGC(dtp$time,dtp$cell.count, dtp[,uid], signif(dtp$drug1.conc,2), main=paste(cl,dr), ...)
        out[[df.name]] <- dtp
    }
    if(output) invisible(out)
}

plotFrac <- function(dat,frac.cn='ch2posfrac',...)
{
    x <- dat$time
    y <- dat[,frac.cn]
    uid <- dat$uid
    rep <- signif(dat$drug1.conc,1)
    y.type='log'
    others <- as.list(substitute(list(...)))[-1L]
    arglist <- list(x,y,uid,rep,y.type)
    arglist <- append(arglist,others)
    do.call(plotGC, arglist)
}

plotAllCh2 <- function(datlist=pd, toFile=FALSE)
{
    lapply(names(datlist), function(cl) 
        {
            prepMultiGCplot(toFile=toFile,fn=paste('Ch2PosPlots_',cl,'.pdf',sep=''))
            drugs <- setdiff(unique(datlist[[cl]]$drug1),'control')
            lapply(drugs, function(dr) 
            {
                dat <- ssCLD(datlist[[cl]],cl=cl,drug=dr)
                dat <- dat[order(dat$drug1.conc),]
                if(is.null(dat$ch2posfrac)) dat$ch2posfrac <- dat$ch2.pos/dat$cell.count
                ylim <- c(0,1)
                main <- paste(cl,substr(dr,1,15),sep='+')
                ylab <- expression(Fraction~of~Ch2^'+'~cells)
                xlab <- "Time (h)"
                myargs <- list(dat=dat,ylim=ylim,main=main, xlab=xlab, ylab=ylab)
                do.call(plotFrac,myargs)
                return()
            })
            if(toFile) dev.off()
            return()
        }
    )
}
