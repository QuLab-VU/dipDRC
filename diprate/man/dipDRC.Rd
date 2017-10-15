\name{dipDRC}
\alias{dipDRC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dipDRC(dtf, xName = "time", yName = "cell.count", var = c("cell.line", "drug1", "drug1.conc", "expt.date"), print.dip = FALSE, norm = FALSE, plotIt = TRUE, toFile = FALSE, fct = LL.4(), uidName = "uid", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtf}{
%%     ~~Describe \code{dtf} here~~
}
  \item{xName}{
%%     ~~Describe \code{xName} here~~
}
  \item{yName}{
%%     ~~Describe \code{yName} here~~
}
  \item{var}{
%%     ~~Describe \code{var} here~~
}
  \item{print.dip}{
%%     ~~Describe \code{print.dip} here~~
}
  \item{norm}{
%%     ~~Describe \code{norm} here~~
}
  \item{plotIt}{
%%     ~~Describe \code{plotIt} here~~
}
  \item{toFile}{
%%     ~~Describe \code{toFile} here~~
}
  \item{fct}{
%%     ~~Describe \code{fct} here~~
}
  \item{uidName}{
%%     ~~Describe \code{uidName} here~~
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
function (dtf, xName = "time", yName = "cell.count", var = c("cell.line", 
    "drug1", "drug1.conc", "expt.date"), print.dip = FALSE, norm = FALSE, 
    plotIt = TRUE, toFile = FALSE, fct = LL.4(), uidName = "uid", 
    ...) 
{
    if (plotIt & toFile) 
        pdf("dipDRC_graph.pdf")
    concName <- var[grep("[Cc]onc", var)]
    exptID <- var[grepl("[Dd]ate", var) | grepl("[Ii][Dd]", var)][1]
    Uconc <- sort(unique(dtf[, concName]), decreasing = TRUE)
    dip.rates <- dtf[, c(var, uidName)]
    dip.rates <- dip.rates[!duplicated(dip.rates[, uidName]), 
        ]
    rownames(dip.rates) <- NULL
    dip.rates <- cbind(dip = NA, dip.rates)
    for (id in unique(dtf[, uidName])) {
        dip.rates[dip.rates[, uidName] == id, "dip"] <- tryCatch({
            findDIP(dtf[dtf[, uidName] == id, c(xName, yName)], 
                print.dip = print.dip)$dip
        }, error = function(cond) {
            message(paste("Error in dipDRC/findDIP:", uidName, 
                "=", id))
            message(cond)
            return(NA)
        })
    }
    dip.rates$norm.dip <- dip.rates$dip/mean(dip.rates[dip.rates[, 
        concName] == 0, "dip"])
    if (norm) {
        f <- formula(paste0("norm.dip ~ ", concName))
        out <- tryCatch({
            drm(f, data = dip.rates, fct = fct)
        }, error = function(cond) {
            return(dip.rates)
        })
    }
    else {
        f <- formula(paste0("dip ~ ", concName))
        out <- tryCatch({
            drc::drm(f, data = dip.rates, fct = fct)
        }, error = function(cond) {
            return(dip.rates)
        })
    }
    if (plotIt & class(out) == "drc") {
        plot.dipDRC(out, ...)
        abline(h = 0, col = grey(0.5), lty = 2)
    }
    if (plotIt & class(out) != "drc") {
        temp <- out
        temp[temp$drug1.conc == 0, "drug1.conc"] <- min(temp[temp$drug1.conc != 
            0, "drug1.conc"])/10
        myargs <- list(...)
        myargs["type"] <- NULL
        myargs <- c(list(formula("dip ~ log10(drug1.conc)"), 
            data = temp), myargs)
        do.call(plot, myargs)
        abline(h = 0, col = grey(0.5), lty = 2)
        legend("bottomleft", "No DRC fit", bty = "n")
    }
    if (plotIt & toFile) 
        dev.off()
    invisible(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
