assembleData <- function(data.loc=getwd())
{
    #' Assemble DIP rates data
    #'
    #' Deprecated - to be removed
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
	#' Check control data
	#'
	#' Helper function calling controlQC; deprecated - to be removed
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
    #' Prepare cell count data by quality-contrl filtering control wells
    #'
    #' Function to process cell count data to reduce control data to a single well
    #'
    #' @param dat data.frame of data to be prepared
    #' @param time.count.id.names character vector of column names for times, cell counts and ids
    #' @param max.cell.count numeric of maximum acceptable value of cell counta
    #' @param cell.line.name character of cell line name
    #' @param ctrl.type character of type of control data aggregation
    #'

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
    #' Prepare cell counts by cell line
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
    #' Prepare device window or file for multiple plots of growth curves
    if(toFile)    pdf(width=w,height=h, file=fn)
    if(!toFile)    dev.new(width=w,height=h)
    par(mfrow=c(3,5), mar=c(2,3,1.5,.5))

    std.par <- list(ylim=c(-2,6), xlab=NA, ylab=NA)
    invisible(std.par)
}

prepMultiDRCplot <- function(w=11,h=8.5, toFile=FALSE, fn='multiDRM.pdf', ...)
{
    #' Prepare device window or file for multiple plots of dose-response curves
    if(!toFile) dev.new(width=w,height=h)
    if(toFile) pdf(width=w,height=h,file=fn)

    par(mfrow=c(3,5), mar=c(2,3,1.5,.5))
}


assembleDIP <- function(fdo, od, var=c('cell.line','drug1','drug1.conc'), name='uid')
{
    #' Function to extract DIP rates from findDIP object (list)
    #'
    #' Function takes a list of findDIP objects (list), extracts the original data used to estimate
    #'  the DIP rate, and assembles them into a data.frame of same structure as original data
    #' @param fdo list of findDIP objects (list)
    #' @param od data.frame of original data sent to findDIP
    #' @param var character vector of variables
    #' @param name character of column name of unique identifier for each DIP rate (default is 'uid')
    #'
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
    #' Get all dose-response models
    #'
    #' Function will identify unique conditions based on \emph{drug1} and \emph{cell.line}
    #' First, the name of each object in the list \emph{dat} will be split on the character '-'
    #'  assuming the first character string corresponds to \emph{cell.line} and
    #' @param dat named list of data.frames containing cell count data
    #' @param plotIt logical indicating whether to plot graphs of dose-response models
    #' @param toFile logicsal indicating whether to save the plotted graphs to a file
    #' @param ... Remaining arguments will get passed to dipDRC function
    #'
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
    #' Subset data by cell line and drug
    # return dat selected by cell line and drug (including controls)
    # controls labeled as "control" in drug (character object)
        dat[dat$cell.line==cl & (dat$drug1==drug | dat$drug1=='control'),]
}
