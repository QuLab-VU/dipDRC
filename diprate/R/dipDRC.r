

# fifth order polynomial
.p5 <- function(x,int,b1,b2,b3,b4,b5) int+b1*x+b2*x^2+b3*x^3+b4*x^4+b5*x^5

.deriv_poly_coef <- function(co) {
    stopifnot(
        all(sapply(names(co), FUN=function(x) x %in% c('int','b1','b2','b3','b4','b5')))
    )
    new.co <- c(co['b1'],2*co['b2'],3*co['b3'],4*co['b4'],5*co['b5'],0)
    new.co[is.na(new.co)] <- 0
    names(new.co) <- c('int','b1','b2','b3','b4','b5')
    new.co
}

# obtaining best fit parameters of a 5th order polynomial
# linear model to 5th order polynomial using p5 function defined above
.fit_p5 <- function(dtf)
{
    x     <- colnames(dtf)[1]
    y     <- colnames(dtf)[2]

    form <- formula(paste(y, '~ .p5(', x, ', int, b1, b2, b3, b4, b5)'))

    m.p5 <- nls(form, data=dtf, start=list(int=1,b1=1,b2=1,b3=1,b4=1,b5=1))
    m.p5.coef.1stderiv <- .deriv_poly_coef(coef(m.p5))
    m.p5.coef.2ndderiv <- .deriv_poly_coef(m.p5.coef.1stderiv)

    list(m=m.p5, coef.1stderiv=m.p5.coef.1stderiv, coef.2ndderiv=m.p5.coef.2ndderiv)
}

rmsd <- function(resid)
{
	#' root-mean-square deviation (RMSD) of an estimator == root-mean-square error (RMSE)
	#'
	#' When the estimator is unbiased, RMSD is the square root of variance == standard error
	#'
	#' RMSD represents the sample standard deviation of the differences between predicted values and observed values.
	#' @param resid vector or matrix of residuals (difference between the observed and predicted values)
    sqrt(sum((resid)^2)/length(resid))
}

.makeUCond <- function(dat,var)
{
    #' Make unique condition
    #'
    #' Function to make unique conditions by pasting variables
    #'
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


findDIP <- function(dtf,name='unknown',all.models=FALSE, metric=c('opt','ar2','rmse')[1],
    o=0, dat.type='cell.counts', print.dip=FALSE)
{
    #' Function to determine DIP rate from cell counts over time
    #'
    #' \code{findDIP} uses one of three different metrics to determine the time after which the
    #'  rate of proliferation becomes stable (linear in log scale) and provides an estimate
    #'  of that rate.
    #' @param dtf data.frame of data to fit
    #' @param name name associated with data to fit
    #' @param all.models logical as to whether to return all fit models
    #' @param metric character of \emph{opt}, \emph{ar2}, or \emph{rmse}
    #' @param o numeric value of offset used in polynomial fit for parameter estimation
    #' @param dat.type character of type of data. Default is \emph{cell.counts}
    #' @param print.dip logical as to whether to print DIP rate estimate to stdout
    #'
    if(nrow(dtf)<4 | ncol(dtf)<2)
    {
        message('findDIP requires data frame with nrows >= 4 and ncol >= 2')
        return(NA)
    }

    x     <- dtf[,1]
    if(dat.type=='cell.counts')    {y <- log2(dtf[,2])} else {y <- dtf[,2]}
    n     <- nrow(dtf)-2

    m     <- list()
    rmse  <- numeric()
    ar2   <- numeric()
    p     <- numeric()

    b     <- y
    a     <- x
    i     <- 1

    for(i in seq(n))
    {
    # fit a linear model
        m[[i]] <- lm(b ~ a)
    # root-mean-squared-error of the residuals of the data to the best-fit model
        rmse <- append(rmse,rmsd(m[[i]]$residuals))
    # adjusted R^2 value of the best-fit model
        ar2   <- append(ar2,summary(m[[i]])$adj.r.squared)
    # p-value
        p     <- append(p,summary(m[[i]])$coefficients[2,4])
    # remove the first element from each vector
        a     <- a[-1]
        b     <- b[-1]
    }

    eval.times <- x[seq(n)]

    # fit a 5th order polynomial to the values of rmse for linear models starting at each time point
    # simply used as an way to describe how the values of rmse change as starting time points are dropped
    rmse.p5     <- tryCatch({.fit_p5(data.frame(x=eval.times,rmse=rmse))},error=function(cond){NA})
    rmse.p5.coef <- tryCatch({coef(rmse.p5$m)},error=function(cond){NA})
    # first derivative of the best-fit 5th order polynomial is used to estimate when
    # the change of rmse over time (i.e. the first derivative value at a given time point) approaches zero
    rmse.p5.1std.coef <- tryCatch({rmse.p5$coef.1stderiv},error=function(cond){NA})

    f <- function(...,offset) .p5(...) + offset


    # opt is a combined metric for choosing the best linear model using
    # adjusted R2 value, rmse of the residuals and fraction of the total data points used in model
    # this weightings of each component of this metric were optimized empirically!
    opt <- ar2 * (1-rmse)^2 * (length(eval.times)-seq(eval.times))^0.25

    # idx is the index of the time after which all data points are included
    # for the best fit linear model
    # three different metrics are provided for identifying this time, adjusted R2 (ar2),
    # root mean square error of the residuals (rmse) or the opt metric described above
    idx <- switch(metric,
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


    dip <- coef(m[[idx]])[2]
    ci <- diff(confint(m[[idx]])[2,])/2
    names(ci) <- intToUtf8(177)     # plus-minus symbol (ascii)
    if(print.dip) print(paste('DIP =',round(dip,4),'starting at',round(x[idx],2),'with',n+2,'data points'))
    out=list(data=data.frame(Time_h=x,l2=y), model=m, metric.used=metric, n=n, idx=idx, best.model=m[[idx]], opt=opt,
        eval.times=x[seq(n)], rmse=rmse,ar2=ar2,p=p,start.time=x[idx],dip=dip, dip.95ci=ci, rmse.p5.1std.coef=rmse.p5.1std.coef)
    if(!all.models)
    {
        out[["model"]] <- NULL
    }
    invisible(out)
}

.sumRep <- function(count,ids)
{
    uid     <- unique(ids)
    sums <- sapply(unique(ids), FUN=function(x) sum(count[ids==x]))
    cbind(uid,sums)
}


dipDRC <- function(dtf, xName='time', yName='cell.count', var=c('cell.line','drug1','drug1.conc','expt.date'),
    print.dip=FALSE, norm=FALSE, plotIt=TRUE, toFile=FALSE, fct=LL.4(), uidName='uid', ...)
{
    #' dipDRC: main function for finding DIP rates from structured data.frame
    #'
    #' Function uses data.frame of cell counts over time and associates other variables with
    #'  resultant DIP rate estimates.
    #'
    # Function to extract DIP rate from single cell line and single drug + control
    # and calculates a 4-param log-logistic fit by default
    if(plotIt & toFile)    pdf('dipDRC_graph.pdf')
    concName <- var[grep('[Cc]onc',var)]
    exptID     <- var[grepl('[Dd]ate',var) | grepl('[Ii][Dd]',var)][1]
    Uconc     <- sort(unique(dtf[,concName]),decreasing=TRUE)
    dip.rates <- dtf[,c(var,uidName)]
    dip.rates <- dip.rates[!duplicated(dip.rates[,uidName]),]
    rownames(dip.rates) <- NULL
    dip.rates <- cbind(dip=NA,dip.rates)
    for(id in unique(dtf[,uidName]))
    {
        dip.rates[dip.rates[,uidName]==id,'dip'] <-
            tryCatch({findDIP(    dtf[dtf[,uidName]==id,c(xName,yName)], print.dip=print.dip)$dip
                },
                error=function(cond) {
                    message(paste("Error in dipDRC/findDIP:",uidName,'=',id))
                    message(cond)
                    return(NA)
                }
            )
    }

    dip.rates$norm.dip <- dip.rates$dip/mean(dip.rates[dip.rates[,concName]==0,'dip'])
    if(norm)
    {
        f <- formula(paste0('norm.dip ~ ',concName))
        out <- tryCatch({drm(f,data=dip.rates,fct=fct)},error=function(cond) {return(dip.rates)})
    } else
    {
        f <- formula(paste0('dip ~ ',concName))
        out <- tryCatch({drc::drm(f,data=dip.rates,fct=fct)    },error=function(cond) {return(dip.rates)})
    }
    if(plotIt & class(out)=='drc')
    {
        plot.dipDRC(out, ...)
        abline(h=0, col=grey(0.5), lty=2)
    }
    if(plotIt & class(out)!='drc')
    {
        temp <- out
        temp[temp$drug1.conc==0,'drug1.conc'] <- min(temp[temp$drug1.conc!=0,'drug1.conc'])/10
        myargs <- list(...)
        myargs['type'] <- NULL
        myargs <- c(list(formula("dip ~ log10(drug1.conc)"), data=temp),myargs)
        do.call(plot,myargs)
        abline(h=0, col=grey(0.5), lty=2)
        legend("bottomleft",'No DRC fit',bty='n')
    }

    if(plotIt & toFile) dev.off()
    invisible(out)
}



plot.dipDRC <- function(drmo, plot.type='confidence', showEC50=TRUE, ...)
{
    #' Plot DIP drc (dose-response curve)
    #'
    #'
    if(is.na(drmo[1]))
    {
        message('DRM object not available to plot')
        return(NA)
    }
    param <- list(...)

    cell.line <- unique(drmo$origData$cell.line)[1]
    drug <- unique(drmo$origData$drug1)[1]

    if(!('main' %in% names(param)))
    {
        tit <- paste(cell.line,drug,sep='_')
        param[['main']] <- tit
    }
    if(!('type' %in% names(param)))        param[['type']] <- plot.type

    myargs <- append(list(x=drmo),param)
    do.call(plot,args=myargs)
    if(showEC50) abline(v=ED(drmo,50,interval='delta',display=FALSE)[1],col='red')
    abline(h=0, col=grey(0.5), lty=2)
}

myll4 <- function(x,b,c,d,e)
{
	#' Four-parameter log-logistic function
	#'
	#                           d - c
	# f(x)  =  c  +  ---------------------------
	#                1 + exp(b(log(x) - log(e)))
	#    b: Hill coefficient
	#    c: lower limit (Emax)
	#    d: upper limit (Emin)
	#    e: EC50 (not log scaled)
    c + ( (d - c) / (1 + exp(b*(log(x) - log(e)))))
}

addLL4curve <- function (drmodel, fromval = 1e-12, toval = 1e-05, norm = FALSE, ...)
{
    #' Add 4-parameter log-logistic model curve
    #'
    #' @param drmodel either \code{drc} or \code{numeric} of \code{drc} paramaters
    #' @param fromval numeric of lower limit of concentration range to plot
    #' @param toval numeric of upper limit of concentration range to plot
    #' @param norm logical of whether to scale the curve to max value of 1
    #'
    #' Function to add a curve of a 4-parameter log-logistic function to a preexisting plot
    #'
    if(class(drmodel)=='drc') param <- coef(drmodel)
    if(class(drmodel)=='numeric') param <- drmodel
    names(param) <- letters[2:5]
    if (norm) {
        param[2] <- param[2]/param[3]
        param[3] <- 1
        if (param[2] > param[3]) {
            temp <- param[2]
            param[2] <- param[3]
            param[3] <- temp
        }
    }
    curve(do.call(myll4, args = append(list(x), as.list(param))),
        from = fromval, to = toval, add = TRUE, ...)
}


getAA <- function(p, drugconcrange=c(1e-12,1e-5), minval=-1, norm=TRUE, removeNE=FALSE)
{
    #' Calculate activity area of dose-response model (drm object)
    #'
    #'
    # p is dose-response model (drm) from the drc library
    # response values of less than minval will be replaced with minval
    # removeNE is logical indicating whether to remove AA values < 0 (no inhbitory effect)

    # normalize response values to 0 as maximum level by subtracting E0 from all values
    # total area determined by dividing by the number of measurements
    # This approach is equivalent to that used by Barretina et al. except that higher values
    # correspond to greater effectiveness of the drug
    # See doi:10.1038/nature11735
    if(class(p)=='drc') p <- coef(p)
    if(!(class(p)=='numeric' & length(p)==4))
    {
        message('argument to getAAr must be drm object or its 4 coefficients (assuming LL.4)')
        return(NA)
    }
    names(p) <- letters[2:5]

    # b = slope parameter (Hill coefficient); c = Emax; d = E0; e = EC50
    # E0 == effect in absence of drug
    if(norm)
    # normalize by dividing Emax by E0 and setting E0 to 1
    {
        p[['c']] <- p[['c']]/p[['d']]
        p[['d']] <- 1
    }
    xvals <- 10^seq(log10(drugconcrange[1]),log10(drugconcrange[2]),.1)
    yvals <- do.call(myll4, args=append(list(x=xvals),as.list(p)))-1
    yvals <- sapply(yvals, function(x) ifelse(x < minval , minval, x))
    # to eliminate activity area measurements less than zero (no effect/enhancing proliferation)
    # is slope is less than zero, increasing drug increases proliferation rate --> set output val to 0
    if(removeNE)
    {
        sapply(seq(length(p)),
            function(x) ifelse(p[['b']][x]<0,0,(-sum(yvals) / length(yvals))[x]))
    } else { -sum(yvals) / length(yvals) }
}

getAAr <- function(p,minmax=c(1e-12,1e-5),RespRatio=TRUE)
{
    #' Calculate activity area of dose-response model (drm object) with response ratio as effect
    #'
    #'
    # p is either a dose-response model (drm) from the drc library or its coefficients
    # formula derived by Leonard Harris; https://www.overleaf.com/9362253wtqkzmhwndsz#/33823251/
    # RespRatio = logical determining whether to scale values so minimum effect = 1 (response ratio)
    # if value is less than 0, set to 0
    if(class(p)=='drc') p <- coef(p)
    if(!(class(p)=='numeric' & length(p)==4))
    {
        message('argument to getAAr must be drm object or its 4 coefficients (assuming LL.4)')
        return(NA)
    }
    names(p)     <- c('h','Emax','E0','EC50')
    out <- 1/p['h'] *  log10(
                (p['EC50']^p['h'] + minmax[2]^p['h']) /
                (p['EC50']^p['h'] + minmax[1]^p['h'])
            )
    #
    if(RespRatio) out <- out * (p['E0']-p['Emax'])/p['E0']
    names(out) <- 'AAr'
    out[out < 0] <- 0
    out
}

getParam <- function(drmod)
{
    #' Extract parameters from dose-response models (drm objects)
    #'
    #' Function to extract relevant parameter values and other extracted metrics from
    #'  dose-response model (drm) objects generated by the \code{drc} package used by
    #'  \code{diprate}.
    #' @param drmod \emph{drm} object
    #' @return data.frame of extracted parameters and metrics (with confidence intervals)
    #'
    if(class(drmod) != 'drc') {message('getParam() requires drm object'); return(invisible(NA))}
    # obtain model formula to determine name of dependent variable
    yname <- drmod$dataList$names$orName
    xname <- drmod$dataList$names$dName

    # estimates and confidence intervals of log-logistic model fit parameters
    ci <- as.data.frame(confint(drmod))
    ci <- cbind(est=coef(drmod),ci)
    rownames(ci) <- c('slope','Emax','E0','EC50')
    colnames(ci) <- c('est','lower','upper')

    # rrEmax = max DIP rate scaled to E0
    rrEmax <- ci['Emax',]/ci['E0','est']
    names(rrEmax) <- c('est','lower','upper')
    # e_halfmax = DIP rate at EC50
    e_halfmax <- predict(drmod, data.frame(drug1.conc=ci['EC50','est']), interval='confidence')
    names(e_halfmax) <- c('est','lower','upper')
    # rr_e_halfmax = response ratio at EC50 (see Harris et al., Nature Methods, 2016, Supp Fig 1)
    # ci may not be calculated correctly here
    rr_e_halfmax <- (ci['Emax',]/ci['E0','est']) +
        ((ci['E0','est']-ci['Emax','est']) / (2*ci['E0','est']))
    names(rr_e_halfmax) <- c('est','lower','upper')
    ic10 <- ED(drmod, ci['E0','est']*0.9, type='absolute', interval='delta', display=FALSE)
    rownames(ic10) <- 'IC10'
    ic50 <- ED(drmod, ci['E0','est']/2, type='absolute', interval='delta', display=FALSE)
    rownames(ic50) <- 'IC50'
    ic100 <- ED(drmod, 0, type='absolute', interval='delta', display=FALSE)
    rownames(ic100) <- 'IC100'
    aa <- getAAr(drmod)

    dat <- drmod$origData
    # lowest observed value (Emax(observed) in the highest two concentrations of drug)
    emaxobs <- min(dat[dat[,xname] >= max(dat[,xname])/11,yname])
    # lowest relative observed value (Emax(observed) in the highest two concentrations of drug)
    rremaxobs <- ifelse(yname=='dip', min(dat[dat[,xname] >= max(dat[,xname])/11,'norm.dip']),NA)

    out <- rbind(ci,rrEmax,e_halfmax,rr_e_halfmax,ic10[c(1,3,4)],ic50[c(1,3,4)],ic100[c(1,3,4)],
        emaxobs,rremaxobs,aa)
    rownames(out) <- c(rownames(out)[1:4],'rrEmax','E50','rrE50','IC10','IC50','IC100','Emaxobs','rrEmaxobs','AAr')
    out
}

plotMultiCurve <- function(drm_list, norm=FALSE, leg.scale=0.75, ...)
{
    #' Plot multiple dose-response curves on single graph
    #'
    #'
    if(class(drm_list) != 'list' | !all(lapply(drm_list,class) == 'drc'))
    {
        message('plotMultiCurve needs list of drm objects')
        return(invisible(NA))
    }
    pargs <- as.list(substitute(list(...)))[-1L]
    drm_cl <- unlist(lapply(drm_list, function(x) unique(x$origData$cell.line)))
    drm_dr <- unlist(lapply(drm_list, function(x) unique(x$origData$drug1)[unique(x$origData$drug1) != 'control']))
    if(!'main' %in% names(pargs))
    {
        pargs['main'] <- ifelse(all(drm_dr==drm_dr[1]),drm_dr[1],'')
    }
    drm_names <- paste(drm_cl,drm_dr)
    if(all(drm_dr==drm_dr[1]))    drm_names <- drm_cl

    if(!'ylim' %in% names(pargs))
    {
        if (norm) { pargs <- append(pargs,list(ylim = c(-0.2,1.2))) } else { pargs <- append(pargs,list(ylim = c(-0.02,0.07))) }
    }
    par(mar=c(4,4,2,.5))
    n <- length(drm_list)
    mycolors <- c('red','orange','yellow','green','blue','purple','brown','black')
    if(n <= 8) { rep.col <-    mycolors[seq(n)] } else {
        rep.col <- mycolors[c(1:8,(seq(n)[8:n]%%8)+1)] }
    rep.lty <- rep(seq(ceiling(n/8)),each=8)[1:n]

    allargs <- append(list(x=drm_list[[1]], xlim=c(0,1e-5), log='x', type='none',
        xlab='[drug], M',ylab='Response ratio', lwd=0, col='white'), pargs)
    do.call(plot, allargs)

    lapply(seq(length(drm_list)), function(x) {
        addLL4curve(drm_list[[x]],
            col=rep.col[x], lty=rep.lty[x], norm=norm, ...); invisible(return()) })
    legend('bottomleft',legend=drm_names,col=rep.col,lty=rep.lty,lwd=leg.scale, cex=leg.scale)
}

effectAtMeanIC50 <- function(drmlist=drc,drugname='erlotinib',sig=3)
{
    #' Effect at mean IC50 value
    #'
    #'
    # find mean IC50 across all cell lines for a particular drug
    # then calculate the relative effect induced by that concentration in each cell line
    drms <- drmByDrug(drmlist,drug_name=drugname)
    ic50s <- sapply(drms, function(x) ED(x, coef(x)['d:(Intercept)']/2, type='absolute', interval='delta', display=FALSE)[1])
    conc <- signif(10^mean(log10(ic50s),na.rm=TRUE),3)
    effect <- sapply(drms, function(x) signif(PR(x,conc)/PR(x,0),sig))
    cl <- sapply(names(drms), function(x) strsplit(x,'.',fixed=TRUE)[[1]][1])
    out <- data.frame(cell.line=cl,drug1=drugname,mean.IC50=conc,rEffect.at.meanIC50=effect)
    rownames(out) <- NULL
    out
}

effectAtMedianEC50 <- function(drmlist=drc,drugname='erlotinib',sig=3)
{
    #' Effect at median EC50 value
    #'
    #'
    # find mean IC50 across all cell lines for a particular drug
    # then calculate the relative effect induced by that concentration in each cell line
    drms <- drmByDrug(drmlist,drug_name=drugname)
    med.ec50 <- signif(median(unlist(sapply(drms, function(x) coef(x)['e:(Intercept)'], simplify=FALSE))),sig)
    effect <- sapply(drms, function(x) signif(PR(x,med.ec50)/PR(x,0),sig))
    cl <- sapply(names(drms), function(x) gsub(paste0('.',drugname),'',x))
    out <- data.frame(cell.line=cl,drug1=drugname,median.EC50=med.ec50,rEffect.at.medianEC50=effect)
    out <- out[order(out$cell.line),]
    rownames(out) <- NULL
    out
}
