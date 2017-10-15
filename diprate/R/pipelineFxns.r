assembleData <- function(data.loc=getwd())
{
    file.names <- dir(data.loc, full.names=TRUE)
    map.names <- file.names[grepl('[Mm]ap',file.names)]
    file.names <- setdiff(file.names[grepl('csv',file.names)],map.names)
    for(fn in file.names)
    {
        map.pre <- unlist(strsplit(fn,'.csv'))
        if(is.na(map.pre))
        {
            message(paste('Could not find a matching plate map for: ',fn))
            return(NA)
        } else {         
            map.name <- map.names[grepl(map.pre,map.names,fixed=TRUE)]
        }
        temp <- cellCountCV(read.csv(fn,as.is=TRUE))
        temp <- addMapInfo3(temp,map.name)
        ifelse(fn==file.names[1], out <- temp,out <- rbind(out,temp))
    }
    out$uid <- paste(out$cell.line,out$drug1,sep='_')
    out$nl2 <- NULL
    colnames(out) <- tolower(colnames(out))
    out
}


checkControlData <- function(ctrl=ad[ad$drug1=='control',],plotIt=TRUE)
{
	# examine control data for each Cell.line
    ctrl.fd <- data.frame()                    # object for filtered data
    for(cl in unique(ctrl$cell.line))
    {
        dfc <- ctrl[ctrl$cell.line==cl,]
        temp <- controlQC(dfc$time, dfc$cell.count, dfc$well, xlim=c(0,150), ylim=c(0,7), cell.line.name=cl, plotIt=plotIt)
        temp$cell.line <- cl
        ctrl.fd <- rbind(ctrl.fd,temp)
    }
    ctrl.fd
}

prepCountData <- function(dat, time.count.id.names=c('time','cell.count','well'), 
    max.cell.count=3000, cell.line.name=unique(dat$cell.line),ctrl.type='mean',...)
{
    tryCatch(
    {
        ctrl <- dat[dat$drug1=='control',]
        dat <- dat[dat$drug1!='control',]
        cn    <- colnames(dat)
        if(length(cell.line.name)>1)
        {
            message('WARNING: Length of prepCountData cell.line.name > 1, using first value only')
            cell.line.name <- cell.line.name[1]
        }

        # replace cell.count values greater than max.cell.count with NA and remove
        dat[dat[time.count.id.names[2]] > max.cell.count,time.count.id.names[2]] <- NA
        dat <- dat[!is.na(dat[,time.count.id.names[2]]),]
    
        # add control data
        ctrl <- tryCatch(
            {    controlQC(ctrl[,time.count.id.names[1]],
                    ctrl[,time.count.id.names[2]],
                    ctrl[,time.count.id.names[3]], 
                    col.names=time.count.id.names, plotIt=FALSE,ctrl.type=ctrl.type,...)
            }, error=function(cond)
            {
                message(paste('error in controlQC parsing data from:',cell.line.name))
                message(cond)
                return(NA)    
            }
        )
        if(length(ctrl)!=1)
        {
            temp <- tail(dat,nrow(ctrl))
            temp[,time.count.id.names[1:2]] <- ctrl
            temp[,!(colnames(temp) %in% c(time.count.id.names,'cell.line','expt.date'))] <- NA
            temp$uid <- paste0(cell.line.name,'_control')
            if("ucd" %in% cn) temp$ucd <- paste0(cell.line.name,'_control')
            if("plate.name" %in% cn & all(dat$plate.name==dat$plate.name[1])) 
                temp$plate.name <- dat$plate.name[1]
            temp$drug1 <- 'control'
            temp$drug1.conc <- 0
            temp[,time.count.id.names[3]] <- 'CTRL'
            return(rbind(dat,temp))
        }
    }, error=function(cond) 
    {
        message(paste('error in prepCountData processing:',cell.line.name))
        message(cond)
        cat('\n')
        return(NA)
    })
}

prepCountByCL <- function(dat)
{
	# determine max cell counts in controls maintaining exponential growth
    out <- list()
    for(cl in unique(dat$cell.line))
    {
        d <- dat[dat$cell.line==cl,]
        max.count <- 
            min(
                findMaxCounts(
                    d[d$drug1=='control',]$time,
                    d[d$drug1=='control',]$cell.count,
                    paste(d[d$drug1=='control',]$cell.line,d[d$drug1=='control',]$well,sep='_'),
                    verbose=FALSE
                )$max.exp.counts
            )
        out[[cl]] <- prepCountData(d, max.cell.count=max.count)
    }
    out
}


prepMultiGCplot <- function(w=11,h=8.5, toFile=FALSE, fn='multiGCplot.pdf', ...)
{
    if(toFile)    pdf(width=w,height=h, file=fn)
    if(!toFile)    dev.new(width=w,height=h)
    par(mfrow=c(3,5), mar=c(2,3,1.5,.5))

    std.par <- list(ylim=c(-2,6), xlab=NA, ylab=NA)
    invisible(std.par)
}

prepMultiDRCplot <- function(w=11,h=8.5, toFile=FALSE, fn='multiDRM.pdf', ...)
{
    if(!toFile) dev.new(width=w,height=h)
    if(toFile) pdf(width=w,height=h,file=fn)

    par(mfrow=c(3,5), mar=c(2,3,1.5,.5))
}

getAllGC <- function(dat, plotIt=TRUE, ...)
{
    # dat is output of prepCountData as list(); name is cell line (may include seeding density)
    # ... are passed to plat function
    ad2 <- list()
    for(n in names(dat)) 
    {
        if(plotIt) prepMultiGCplot()
        ad2 <- append(ad2,plotAllDrugs(dat[[n]], output=plotIt, ...))
    }
    ad2
}

getAllDipRates <- function(dat, uidName='uid', ucdName='ucd')
{
    out <- sapply(unique(dat[,uidName]), simplify=FALSE, FUN=function(x)
        if(grepl('_control',dat[match(x,dat[,uidName])[1],ucdName]))
        {
            temp <- filterCtrlData(dat[dat[,uidName]==x,'time'],dat[dat[,uidName]==x,'cell.count'],dat[dat[,uidName]==x,uidName])
            findDIP(temp[,c('times','counts')], name=x)
        } else
            findDIP(dat[dat[,uidName]==x,c('time','cell.count')], name=x)
    )
    out
}


assembleDIP <- function(fdo, od, var=c('cell.line','drug1','drug1.conc'), name='uid')
{
    # fdo == findDIP object (list)
    # od == original data sent to findDIP
    # names(fdo) should be found in one of the columns in od; default matching colname in od is "uid"
    id.name <- names(fdo)
    dip.rates <- od[firstInstPos(od[,name]),var]
    rownames(dip.rates) <- NULL
    dip.rates[,name] <- names(fdo)
    dip.rates$dip <- unlist(sapply(fdo,FUN=function(x) x['dip']))
    dip.rates$dip.95ci <- unlist(sapply(dip,FUN=function(x) x['dip.95ci']))
    dip.rates$ucd <- makeUCond(dip.rates, c('cell.line','drug1'))
    dip.rates
}

getAllDRM <- function(dat, plotIt=TRUE, toFile=FALSE, ...)
{
    all.drm <- lapply(names(dat), FUN=function(x)
    {
        cl <- unlist(lapply(strsplit(x,'-',fixed=TRUE),function(z) z[1]))
        drugs <- unique(dat[[x]]$drug1)[!grepl('control',unique(dat[[x]]$drug1))]
    
        # prepare plotting area
        prepMultiDRCplot(toFile=toFile, fn=paste('multiDRM_',x,'.pdf',sep=''))

        temp <- sapply(drugs, function(dr) 
        {
            # get data to fit (subset for a specific cell line and drug)
            dtf <- ssCLD(dat[[x]],cl,dr)
            # make a unique name (cell line - drug)            
            n <- unique(paste(cl,dr,sep='_'))
            as.list(dipDRC(dtf, plotIt=plotIt, main=n, ...))
        }, simplify=FALSE)
        if(toFile) dev.off()
        temp
    })
    names(all.drm) <- names(dat)
    invisible(all.drm)
}



ssCLD <- function(dat,cl=unique(dat$cell.line),drug='erlotinib')
{
    # return dat selected by cell line and drug (including controls)
    # controls labeled as "control" in drug (character object)
        dat[dat$cell.line==cl & (dat$drug1==drug | dat$drug1=='control'),]
}
