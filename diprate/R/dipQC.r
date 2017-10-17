
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

filterCtrlData <- function(times, counts, ids, min.ar2=0.99, verbose=FALSE)
{
	#' Filter control data
	#' 
	#' Determine whether cell counts are exponentially increasing throughout the entire time span
	#'  times and counts for which adjusted R2 value is >= min.ar2 argument are returned
	#'  ids = unique identifier for each sample (usually a well from an experiment)
	#'  assumes counts in linear scale (i.e. direct cell counts)
    #'
    #' @param times vector of times
    #' @param counts vector of cell counts
    #' @param ids vector of unique identifiers used to separate groups
    #' @param min.ar2 numeric of minimum value for adjusted R-squared value of linear model fit
    #' @param verbose logical whether to show progress
    #' 
    #' @return data.frame of times, counts, and ids for control data passing filter (i.e.,
    #'  linear (in log scale) with adj R-squared value less than \code{min.ar2})
    #'
    dip   <- numeric()
    dtk.times <- integer()
    dtk.counts <- integer()
    dtk.ids  <- character()
    
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
	# data should be filtered first using filterCtrlData()
	# ids = unique identifier for each sample (usually a well from an experiment)
	# assumes counts in linear scale (i.e. direct cell counts)

    # 
    # types
    #
    # max: identify the sample(s) with the most data points, sum them and fit a lm
    # sum: use all samples trimmed to the fewest data points, sum them and fit a lm
    # mean: use all samples, fitting lm separately to each and return the average

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
    # ensure ordered by time
    counts <- counts[order(times)]
    times <- times[order(times)]
    ifelse(!any(counts > threshold),max(times),min(times[counts > threshold]))
}

controlQC <- function(times, counts, ids, col.names=c('time','cell.counts','uid'), 
    plotIt=TRUE, ctrl.type='mean', cell.line.name="", ret.type='counts', ...)
{
    # expecting data to be from a single cell line
    # will filter data for consistency with exponential growth and return
    # a data.frame with a single set of time points and the estimated cell counts
    # undefined arguments will be passed to plot function
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
    
    times <- times[ids %in% ids.ok]
    counts <- counts[ids %in% ids.ok]
    ids <- ids[ids %in% ids.ok]

    fd <- data.frame(times,counts,ids)
    fd$ids <- as.character(fd$ids)
    
    if(plotIt) 
    {
     plot.args <- append(list(fd[,1],fd[,2],fd[,3],fd[,3], cell.line.name),arglist)
     names(plot.args)[1:5] <- c('x','y','uid','rep','main') 
     do.call(plotGC, plot.args)
    }
    all.out <- findCtrlDIP(fd[,1],fd[,2],fd[,3], col.names=col.names, type=ctrl.type)
    out <- all.out$data
    if(plotIt) lines(out[,1],log2norm(out[,2],ids=1), lwd=3)
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
	# To determine at what density cell counts are no longer increasing exponentially
	# times and counts for which adjusted R2 value is >= min.ar2 argument are returned
	# ids = unique identifier for each sample (usually a well from an experiment)
	# assumes counts in linear scale (i.e. direct cell counts)
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
	# examine well data for each Cell.line
    if(length(drug)>1)
    {
     message('Only data from first drug being returned from prepDataCLD')
    } else {
     return(dat[dat[,drug.col]==drug[1],])
    }
}


getGCargs <- function(dat, arg.name=c('time','cell.count','ids'),dat.col=c('time','cell.count','well'))
{
	# converts data from data frame into a list of arguments that can be passed to other functions
	# this is a temporary function that will require more work to ensure robustness
    out <- list()
    for(a in arg.name) out[[a]] <- dat[,dat.col[match(a,arg.name)]]
    out
}