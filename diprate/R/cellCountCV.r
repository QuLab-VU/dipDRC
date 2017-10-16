cellCountCV <- function(rawCVdata, normPos=1) {
    #' Function to process cell count data in Cellavista multi-measurement format
    #'
    #' @param rawCVdata \emph{data.frame}
    #' @param normPos \emph{integer}
    #'
    #' Function takes \code{data.frame} of raw cell count data exported from Synentec
    #' Cellavista and extracts time of image acquisition, identifies wells from which images
    #' have been obtained, outputs a \code{data.frame} with \code{colnames} of
    #' \emph{Row, Column, Well, Cell.count, nl2}. \emph{nl2} is log2 values of \emph{Cell.count}
    #' normalized to \emph{normPos} argument (default = 1).

    d <- rawCVdata
    d <- d[,c('Column','Row','Cell.Nucleus')]
    # rename Cell.Nucleus to Cell.count
    colnames(d)[3] <- 'Cell.count'

    # ensure that Row and Column are char vectors
    d$Row <- as.character(d$Row)
    d$Column <- as.character(d$Column)

    # make row index of data positions for each time point
    idx <- grep('Timespan',d[,1])

    # extract time of acquisition for each time point
    times <- as.numeric(d[,2][idx])
    time.idx <- seq(length(times))

    # extract the number of wells we have data for
    numWells <- idx[2]-idx[1]-2

    # assemble data into single structure (a)
    for(tp in time.idx)
     {
      temp <- d[(idx[tp]+1):(idx[tp]+numWells),]
      temp$Time  <- times[tp]
      a <- ifelse(tp==time.idx[1], temp, rbind(a,temp))
     }

    rownames(a) <- seq(nrow(a))

    # Generate a column for "Well" in the format RCC
    a[nchar(a$Column)==1,'Column'] <- paste0('0',a[nchar(a$Column)==1,'Column'])
    a$Well  <- paste0(a$Row,a$Column)
    a$Column <- as.integer(a$Column)

    # calculate log2 cell #
    a   <- a[order(a$Well),]
    a   <- a[!is.na(a$Cell.count),]

    a$nl2  <- diprate::log2norm(a$Cell.count, a$Well, norm_idx=normPos)
    rownames(a) <- NULL

    a
}
