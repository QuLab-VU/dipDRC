\name{controlQC}
\alias{controlQC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
controlQC(times, counts, ids, col.names = c("time", "cell.counts", "uid"), plotIt = TRUE, ctrl.type = "mean", cell.line.name = "", ret.type = "counts", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{counts}{
%%     ~~Describe \code{counts} here~~
}
  \item{ids}{
%%     ~~Describe \code{ids} here~~
}
  \item{col.names}{
%%     ~~Describe \code{col.names} here~~
}
  \item{plotIt}{
%%     ~~Describe \code{plotIt} here~~
}
  \item{ctrl.type}{
%%     ~~Describe \code{ctrl.type} here~~
}
  \item{cell.line.name}{
%%     ~~Describe \code{cell.line.name} here~~
}
  \item{ret.type}{
%%     ~~Describe \code{ret.type} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (times, counts, ids, col.names = c("time", "cell.counts", 
    "uid"), plotIt = TRUE, ctrl.type = "mean", cell.line.name = "", 
    ret.type = "counts", ...) 
{
    arglist <- list(...)
    if ("min.ar2" \%in\% names(arglist)) {
        fd.args <- list(times, counts, ids, arglist[["min.ar2"]])
        names(fd.args) <- c("times", "counts", "ids", "min.ar2")
        arglist["min.ar2"] <- NULL
    }
    else {
        fd.args <- list(times, counts, ids)
        names(fd.args) <- c("times", "counts", "ids")
    }
    fd <- do.call(filterCtrlData, fd.args)
    times <- fd$times
    counts <- fd$counts
    ids <- as.character(fd$ids)
    m <- lm(log2(counts) ~ times * ids)
    rates <- coef(m)[grepl("times", names(coef(m)))]
    rates <- c(rates[1], rates[-1] + rates[1])
    ids.ok <- names(rates[rates < mean(rates) + sd(rates) & rates > 
        mean(rates) - sd(rates)])
    ids.ok <- gsub("times:ids", "", ids.ok)
    times <- times[ids \%in\% ids.ok]
    counts <- counts[ids \%in\% ids.ok]
    ids <- ids[ids \%in\% ids.ok]
    fd <- data.frame(times, counts, ids)
    fd$ids <- as.character(fd$ids)
    if (plotIt) {
        plot.args <- append(list(fd[, 1], fd[, 2], fd[, 3], fd[, 
            3], cell.line.name), arglist)
        names(plot.args)[1:5] <- c("x", "y", "uid", "rep", "main")
        do.call(plotGC, plot.args)
    }
    all.out <- findCtrlDIP(fd[, 1], fd[, 2], fd[, 3], col.names = col.names, 
        type = ctrl.type)
    out <- all.out$data
    if (plotIt) 
        lines(out[, 1], log2norm(out[, 2], ids = 1), lwd = 3)
    if (ret.type == "all") 
        out <- list(control.counts = out, control.type = ctrl.type, 
            passed.qc = fd, dip = all.out$dip, model = all.out$model)
    invisible(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
