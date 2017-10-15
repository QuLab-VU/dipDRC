\name{assemPlateData}
\alias{assemPlateData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
assemPlateData(segdirpath, count_re = "cellcount", pids = "", toFile = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segdirpath}{
%%     ~~Describe \code{segdirpath} here~~
}
  \item{count_re}{
%%     ~~Describe \code{count_re} here~~
}
  \item{pids}{
%%     ~~Describe \code{pids} here~~
}
  \item{toFile}{
%%     ~~Describe \code{toFile} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (segdirpath, count_re = "cellcount", pids = "", toFile = FALSE) 
{
    out <- data.frame()
    if (pids[1] != "") {
        pds <- paste0(segdirpath, paste0("/Plate", pids))
        allpds <- basename(list.dirs(segdirpath, full.names = TRUE, 
            recursive = FALSE))
        if (length(setdiff(basename(pds), allpds)) != 0) {
            message(paste("Not all plate ids passed as <pids> found in", 
                segdirpath))
            return(invisible(NA))
        }
    }
    else {
        pds <- list.dirs(segdirpath, full.names = TRUE, recursive = FALSE)
    }
    pb <- txtProgressBar(1, length(pds), style = 3)
    for (p in pds) {
        files <- list.files(p, full.names = TRUE)
        files <- files[grepl(count_re, files)]
        dat <- sapply(files, function(x) read.csv(x, as.is = TRUE), 
            simplify = FALSE)
        dat <- do.call(rbind, dat)
        rownames(dat) <- NULL
        fn <- basename(p)
        if (ncol(dat) > ncol(out)) {
            newcn <- setdiff(colnames(dat), colnames(out))
            temp <- data.frame(matrix(ncol = length(newcn), nrow = nrow(out)))
            colnames(temp) <- newcn
            out <- cbind(out, temp)
            out <- out[, colnames(dat)]
        }
        if (toFile) 
            write.csv(dat, file = file.path(segdirpath, paste0(fn, 
                "_cellcounts.csv")), row.names = FALSE)
        out <- rbind(out, dat)
        setTxtProgressBar(pb, match(p, pds))
    }
    close(pb)
    invisible(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
