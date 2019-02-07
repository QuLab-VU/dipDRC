
nEach <- function(ids)
{
    #' Count the number of each unique item
    #'
    #' Given a vector of items, identify and count how many of each unique item is present
    #'
    #' @param ids vector
    #' @return named integer of counts for each unique item in \code{ids}
    #'
    out <- integer()
    for(i in unique(ids)) out <- append(out,length(ids[ids==i]))
    names(out) <- unique(ids)
    out
}

firstInstPos <- function(vec)
{
    #' Find the position of the first occurrence of each unique item
    #'
    #' Given a vector of items, identify the first position of each unique item
    #'
    #' @param vec a character, numeric or integer vector
    #' @return named integer of indexed positions for each unique item in \code{vec}
    #'
    out <- match(unique(vec), vec)
    names(out) <- unique(vec)
    out
}

filterCtrlData <- function(times=NULL, counts=NULL, ids=NULL, dat=NULL, min.ar2=0.99, verbose=FALSE)
{
	#' Filter control data
	#'
	#' Determine whether cell counts are exponentially increasing throughout the entire time span.
	#'  Assumes \code{counts} are in linear scale (i.e. direct cell counts)
    #'
    #' @param times vector of times
    #' @param counts vector of cell counts
    #' @param ids vector of unique identifiers used to separate groups (usually a well from an experiment)
    #' @param dat data.frame of times, cell counts and unique identifiers in columns 1:3 
    #' @param min.ar2 numeric of minimum value for adjusted R-squared value of linear model fit
    #' @param verbose logical whether to show progress
    #'
    #' @return data.frame of times, counts, and ids for control data passing filter (i.e.,
    #'  linear (in log scale) with adj R-squared value less than \code{min.ar2})
    #' 
    #' Input parameters can be either \code{times}, \code{counts}, and \code{ids}, or a \emph{data.frame}
    #'  of these in the first three columns
    #'  
    
    dip   <- numeric()
    dtk.times <- integer()
    dtk.counts <- integer()
    dtk.ids  <- character()

    if(!is.null(dat))
    {
        times <- dat[,1]
        counts <- dat[,2]
        ids <- dat[,3]
    }
    for(i in unique(ids))
    {
     if(verbose) message(paste('Processing',i))

     ttf <- times[ids==i]
     ctf <- log2(counts[ids==i])

     if(length(ttf) < 5)
     {
      message(paste('Need at least 5 data points for ',i))
      next
     }
     mintimept  <- ifelse(length(ttf)/2 > 5,floor(length(ttf)/2),5)

     # m.ad = model with all time points
     m.ad <- lm(ctf ~ ttf)
     ar2 <- summary(m.ad)$adj.r.squared

     for(l in length(ttf):5)
     {
     # include data through the 'l'th time point
      x <- ttf[seq(l)]
      y <- ctf[seq(l)]
     # capture last position where adjusted R^2 value >= min.ar2
     # added in case of perfect line (e.g. no change in cell number for 5 consecutive time points)
      if(!is.na(summary(lm(y ~ x))$adj.r.squared))
      {
       pass <- summary(lm(y ~ x))$adj.r.squared >= min.ar2
      } else {
       pass <- FALSE
      }

      if(pass)
      {
     #remove points that caused failure
       x <- head(x,-1)
       y <- head(y,-1)
       break
      }
     }

     if(!pass)
     {
      if(verbose) message(paste('The first 5 pts of ids =',i,'is not linear within ar2 =',min.ar2))
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
    col.names=c('time','cell.counts','well'),
    type='mean')
{
    #' Find control DIP rate
    #'
	#' @param times numeric or difftime vector
    #' @param counts vector of cell counts (assumed in linear scale, i.e. direct cell counts)
	#' @param ids character vector of unique identifier for each sample (usually a well from an experiment)
    #' @param col.names character vector of column names for output
    #' @param type character of function to apply to replicates; default is \emph{mean}
    #'
	#' data should be filtered first using filterCtrlData()
    #'
    #' Other acceptable type values: \itemize{
    #' \item{\emph{max}: identify the sample(s) with the most data points, sum them and fit a lm},
    #' \item{\emph{sum}: use all samples trimmed to the fewest data points, sum them and fit a lm},
    #' \item{\emph{mean}: use all samples, fitting lm separately to each and return the average}}
    #'

    dip   <- numeric()
    dip.temp <- numeric()

    if(type=='mean')
    {
     m <- list()
     for(i in unique(ids))
     {
      x <- times[ids==i]
      y <- log2(counts[ids==i])
      m[[i]] <- lm(y ~ x)
      dip.temp <- append(dip.temp,coef(m[[i]])[2])
      names(dip.temp)[length(names(dip.temp))] <- i
     }
     dip <- mean(dip.temp)
     dip <- append(dip, sd(dip.temp)/sqrt(length(dip.temp)))
     names(dip) <- c('mean.dip','std.err')

     count.t0 <- floor(mean(sapply(unique(ids), FUN=function(x) head(counts[ids==x],1))))
     dts.times <- unique(times)
     dts.counts <- c(count.t0,floor(count.t0*2^(dip[1]*dts.times[-1])))
    }

    if(type=='max')
    {
     n   <- nEach(ids)
     n.max  <- max(n)
     max.ids  <- names(n)[n==n.max]
     dts.times <- times[ids %in% max.ids]
     dts.counts <- counts[ids %in% max.ids]
     a   <- aggregate(dts.counts ~ dts.times, FUN=sum)
     m   <- lm(log2(dts.counts) ~ dts.times, data=a)
     dip   <- coef(m)[2]
     dip   <- append(dip, summary(m)$coefficients[2,2])
     names(dip) <- c('max.dip','std.err')
     dts.times <- a$dts.times
     dts.counts <- a$dts.counts
    }

    if(type=='sum')
    {
     n   <- nEach(ids)
     n.min  <- min(n)
     dts.times <- as.numeric(sapply(unique(ids), FUN=function(x) head(times[ids==x],n.min)))
     dts.counts <- as.numeric(sapply(unique(ids), FUN=function(x) head(counts[ids==x],n.min)))
     a   <- aggregate(dts.counts ~ dts.times, FUN=sum)
     m   <- lm(log2(dts.counts) ~ dts.times, data=a)
     dip   <- coef(m)[2]
     dip   <- append(dip, summary(m)$coefficients[2,2])
     names(dip) <- c('sum.dip','std.err')
     dts.times <- a$dts.times
     dts.counts <- a$dts.counts
    }

    out.data <- data.frame(dts.times,dts.counts)
    colnames(out.data) <- col.names[1:2]

    list(data=out.data,dip=dip,model=m)
}

timeAtTH <- function(times,counts,threshold=1000)
{
    #' Time at threshold
    #'
    #' Function to determine first time a threshold value is achieved
    #' @param times numeric of times of measurements (usually relative to time of treatment)
    #' @param counts integer of number of objects (nuclei) at each measurement time
    #' @param threshold integer of threshold value of \emph{counts}; default = 1000
    #' @return numeric of first value of \emph{times} when \emph{counts} exceeds \emph{threshold}
    
    # ensure ordered by time
    counts <- counts[order(times)]
    times <- times[order(times)]
    ifelse(!any(counts > threshold),max(times),min(times[counts > threshold]))
}

controlQC <- function(times, counts, ids, col.names=c('time','cell.counts','uid'),
    plotIt=TRUE, ctrl.type='mean', cell.line.name="", ret.type='counts', ...)
{
    #' Quality control of control wells
    #' 
    #' Function to identify \emph{ids} for which \emph{counts} are exponential over range of \emph{times}.
    #'  Wrapper for \code{filterCtrlData} and \code{findCtrlDIP} functions and other filtering.
    #' First, data are passed to \code{filterCtrlData} which finds linear model fits to data for each
    #'  unique \emph{ids} with at least 5 data points and with minimum adjusted R-squared value 
    #'  greater than 0.99 (default) or the value of the argument \emph{min.ar2}.
    #' Data passing QC can also be plotted (default). Data are expected to be from a single cell line.
    #'  Undefined arguments will be passed to plot function.
    #' 
    #' @param times numeric of times when \emph{counts} were obtained
    #' @param counts integer of cell (nuclei) counts
    #' @param ids character of unique identifiers for each set of \emph{times} and \emph{counts}
    #' @param col.names character of names sent as arguments to \code{filterCtrlData}; default is 
    #'  \code{c('time','cell.counts','uid')}
    #' @param plotIt logical for whether to plot data; default is \code{TRUE}
    #' @param ctrl.type character of function to perform on control data; acceptable values include
    #'  \emph{sum}, \emph{max}, and default is \code{mean}.
    #' @param cell.line.name character of name of cell line from which control data were obtained;
    #'  default is an empty character
    #' @param ret.type character of type of values to return; default is \emph{counts}, which returns
    #'  a \code{data.frame} of \emph{times}, \emph{counts}, and \emph{ids} for control data that passes QC; 
    #'  if value is \emph{all}, returns \code{list} of \emph{control.counts}, a \code{data.frame}
    #'  in same format as for \emph{counts}, \emph{control.type} a character of \emph{ctrl.type}
    #'  argument, \emph{passed.qc} list of arguments passed to \code{findCtrlDIP}, \emph{dip}
    #'  numeric of estimated proliferation rates, and \emph{model} \code{lm} model fit to data. 
    #' 
    #' @return data.frame with a single set of time points and the estimated cell counts (if \code{ret.type} != 'all')
    #' @return list of integer \code{control.counts}, character \code{control.type}, logical \code{passed.qc},
    #'  numeric \code{dip rate}, linear model \code{model} (if \code{ret.type} == 'all') 
    #'  

    arglist <- list(...)
    if('min.ar2' %in% names(arglist))
    {
        fd.args <- list(times,counts,ids,arglist[['min.ar2']])
        names(fd.args) <- c('times','counts','ids','min.ar2')
        arglist['min.ar2'] <- NULL
    } else {
        fd.args <- list(times,counts,ids)
        names(fd.args) <- c('times','counts','ids')
    }
    fd <- do.call(filterCtrlData,fd.args)

    times <- fd$times
    counts <- fd$counts
    ids <- as.character(fd$ids)

    # remove outlier ids ( > 1 sd from the mean)
    #
    # ids = unique identifier for each sample (usually a well from an experiment)
    # assumes counts in linear scale (i.e. direct cell counts)
    m <- lm(log2(counts) ~ times * ids)
    rates <- coef(m)[grepl('times',names(coef(m)))]
    rates <- c(rates[1],rates[-1]+rates[1])

    ids.ok <-  names(rates[rates < mean(rates)+sd(rates) & rates > mean(rates)-sd(rates)])
    ids.ok <- gsub('times:ids','',ids.ok)
    # need to replace first ID since used as comparator in lm fits
    ids.ok <- sub('times',ids[1],ids.ok)

    times <- times[ids %in% ids.ok]
    counts <- counts[ids %in% ids.ok]
    ids <- ids[ids %in% ids.ok]

    fd <- data.frame(times,counts,ids)
    fd$ids <- as.character(fd$ids)

    if(plotIt)
    {
        plot.args <- append(list(fd[,1],fd[,2],fd[,3],fd[,3], cell.line.name),arglist)
        names(plot.args)[1:5] <- c('time','cell.count','ids','rep','main')
        do.call(plotGC, plot.args)
    }
    all.out <- findCtrlDIP(fd[,1],fd[,2],fd[,3], col.names=col.names, type=ctrl.type)
    out <- fd
    if(ret.type=='all')
    out <- list(
        control.counts=out,
        control.type=ctrl.type,
        passed.qc=fd,
        dip=all.out$dip,
        model=all.out$model)
    invisible(out)
}




findMaxCounts <- function(times, counts, ids, min.ar2=0.99, verbose=TRUE)
{
    #' Function to determine cell density above which cells do not proliferate exponentially
    #'
    #' To determine at what density cell counts are no longer increasing exponentially
	#'  times and counts for which adjusted R2 value is >= min.ar2 argument are returned
	#' @param times numeric or difftime vector
    #' @param counts vector of cell counts (assumed in linear scale, i.e. direct cell counts)
	#' @param ids character vector of unique identifier for each sample (usually a well from an experiment)
    #' @param min.ar2 numeric of minimum value of adjusted R2 of linear model fit
    #'  to consider proliferation as exponential
    #' @param verbose logical of amount of information sent to stdout during processing
    #'
	#' Variable \emph{counts} is assumed to be in linear scale (i.e. direct cell counts)
    #'
    max.exp.counts <- integer()
    max.ids  <- character()
    max.count.time <- numeric()

    for(i in unique(ids))
    {
     ttf <- times[ids==i]
     ctf <- log2(counts[ids==i])

     if(length(ttf) < 5)
     {
      message(paste('Need at least 5 data points for ',i))
      next
     }
     mintimept  <- ifelse(length(ttf)/2 > 5,floor(length(ttf)/2),5)

     for(l in mintimept:(length(ttf)))
     {
     # include data through the 'l'th time point
      x <- ttf[seq(l)]
      y <- ctf[seq(l)]
     # capture last position where adjusted R^2 value <= min.ar2
      if(summary(lm(y ~ x))$adj.r.squared <= min.ar2) break
     }

     if(summary(lm(y ~ x))$adj.r.squared < min.ar2 & l == 5)
     {
      if(verbose) message(paste('The first 5 pts of ids =',i,'is not linear within ar2 =',min.ar2))
      next
     } else {
      max.exp.counts <- append(max.exp.counts,floor(2^y[l]))
      max.ids <- append(max.ids, i)
      max.count.time <- append(max.count.time,x[l])
     }
    }

    out <- data.frame(max.exp.counts=max.exp.counts,ids=max.ids,max.count.time=max.count.time)
    out$ids <- as.character(out$ids)
    out
}

prepDataCLD <- function(dat, drug, drug.col="drug1")
{
    #' Prepare data by cell line and drug
    #'
	# examine well data for each Cell.line
    if(length(drug)>1)
    {
     message('Only data from first drug being returned from prepDataCLD')
    } else {
     return(dat[dat[,drug.col]==drug[1],])
    }
}


getGCargs <- function(dat, arg.name=c('time','cell.count','ids'),dat.col=c('time','cell.count','uid'))
{
    #' Get growth curve arguments
    #'
	#' Converts data from \code{data.frame} into a \code{list} of arguments that can be
	#'  passed to other functions.
	#' @param dat data.frame containing colnames in \emph{dat.col}
	#' @param arg.name character of names for output \code{list} of arguments; default is
	#'  c('time','cell.count','ids')
	#' @param dat.col character of colnames present in \emph{dat}; default is
	#'  c('time','cell.count','uid')
	#' @return list of arguments with names from \emph{arg.name}
	# this is a temporary function that will require more work to ensure robustness
    out <- list()
    for(a in arg.name) out[[a]] <- dat[,dat.col[match(a,arg.name)]]
    out
}

aboveMax <- function(times=NULL, counts=NULL, ids=NULL, dat=NULL, max.count=3000)
{
    #' Function to identify vector positions where cell.count > max.count
    #' @param times numeric of times at which cell.count was obtained
    #' @param counts int of number of cells at each \emph{times}
    #' @param ids character of unique identifiers for \emph{times} and \emph{counts}
    #' @param dat data.frame of \emph{times}, \emph{counts} and \emph{ids} in first three columns
    #' @param max.count maximum value of \emph{counts} over which the output value is \code{FALSE}
    #' @return logical of length == \code{length(times)} whether \emph{counts} exceeds \emph{max.count}
    #' 
    #' NOTE: if \emph{dat} is provided, \emph{times}, \emph{counts} and emph{ids} will not be used
    #' @examples
    #' aboveMax( times=rep(c(1:5),2), counts=c((1:5)*1000,(1:5)*700), ids=c(rep('A',5),rep('B',5)) )
    #' 
    if(all(!exists(c('times','counts','ids','dat'))))
    {
        message('aboveMax() needs some data')
        return(NA)
    }
    if(exists('dat') & class(dat)=='data.frame')
    {
        message('aboveMax(): Assuming first three columns in <dat> are: time, cell.count, and uid')
        o <- order(dat[,3],dat[,1])
        if(!all(o == seq(nrow(dat))))
        {
            message('dat must be ordered by ids and time')
            return(NA)
        } else {
            times <- dat[,1]
            counts <- dat[,2]
            ids <- dat[,3]
        }
    }
    out <- rep(FALSE,length(times))
    for(id in unique(ids))
    {
        # find time after which cell.count > max.count
        # if no cell.count > max.count keep all time points
        if(all(counts[ids == id] < max.count)) next
        max.time <- min(times[ids == id & counts >= max.count])
        out[ids == id] <- times[ids == id] >= max.time
    }
    return(out)
}


okControls <- function(times=NULL, counts=NULL, ids=NULL, dat=NULL, minr2=0.95)
{
    #' Function to determine control data demonstrated exponential growth
    #'  over entire duration, based on minimum R2 value of linear model fit
    #' @param times numeric of times at which cell.count was obtained
    #' @param counts integer of number of cells at each \emph{times}
    #' @param ids character of unique identifiers for \emph{times} and \emph{counts}
    #' @param dat data.frame of \emph{times}, \emph{counts} and \emph{ids} in first three columns
    #' @param minr2 numeric of minimum value of R2 needed for all linear models of data
    #'  must be above; default is 0.95
    #' @return logical of whether all \code{lm} fits to data exceeds \emph{minr2}
    #' 
    #' NOTE: if \emph{dat} is provided, \emph{times}, \emph{counts} and emph{ids} will not be used
    #' @examples
    #' dat <- data.frame( times=1:5, counts=500*2^(1:5), ids=rep('A',5) )
    #' okControls(  dat=dat )
    #' dat <- data.frame( times=1:5, counts=c(500*2^(1:4), 5000), ids=rep('A',5) )
    #' okControls(  dat=dat )
    #' 
    if(all(sapply(c('times','counts','ids','dat'),is.null)))
    {
        message('okControls() needs some data')
        return(NA)
    }
    if(!is.null('dat'))
    {
        if(class(dat)!='data.frame' | suppressWarnings(ncol(dat) < 3))
        {
            message('dat arg to okControls needs a data.frame with at least 3 columns')
            return(NA)
        }
        cn <- colnames(dat)
        if(sum(length(grep('time',cn)),length(grep('count',cn)),length(grep('id',cn)))==3)
        {
            # matching colnames(dat) to expected names
            cn[grep('id',cn)] <- 'ids'
            cn[grep('time',cn)] <- 'times'
            cn[grep('count',cn)] <- 'counts'
            colnames(dat) <- cn
        } else {
            message('okControls(): Assuming first three columns in <dat> are: times, counts, and ids')
            o <- order(dat[,3],dat[,1])
            if(!all(o == seq(nrow(dat))))
            {
                message('dat must be ordered by ids and time')
                return(NA)
            }
            colnames(dat) <- c('times','counts','ids')
        }
    } else {
        dat <- data.frame(times=times,counts=counts,ids=ids)
    }
    m <- lme4::lmList(log2(counts) ~ times | ids, data=dat)
    all(unlist(summary(m)$adj.r.squared) >= minr2)
}
