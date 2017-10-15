\name{plotCellLineDrug}
\alias{plotCellLineDrug}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotCellLineDrug(dat = ad, expt.date = "20160513", toFile = FALSE, max.count = NA, max.count.type = "by.cell.line", plotFit = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{expt.date}{
%%     ~~Describe \code{expt.date} here~~
}
  \item{toFile}{
%%     ~~Describe \code{toFile} here~~
}
  \item{max.count}{
%%     ~~Describe \code{max.count} here~~
}
  \item{max.count.type}{
%%     ~~Describe \code{max.count.type} here~~
}
  \item{plotFit}{
%%     ~~Describe \code{plotFit} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (dat = ad, expt.date = "20160513", toFile = FALSE, max.count = NA, 
    max.count.type = "by.cell.line", plotFit = FALSE, verbose = FALSE) 
{
    out <- data.frame()
    if (is.na(max.count) & max.count.type != "by.cell.line") 
        max.count <- min(findMaxCounts(dat$time, dat$cell.count, 
            paste(dat$cell.line, dat$well, sep = "_"), verbose = verbose))
    for (cl in unique(dat$cell.line)) {
        if (is.na(max.count) & max.count.type == "by.cell.line") {
            dfm <- dat[dat$Cell.line == cl & dat$drug1 == "control", 
                ]
            max.count <- min(findMaxCounts(dfm$time, dfm$cell.count, 
                dfm$well, verbose = verbose))
        }
        d <- prepCountData(dat[dat$cell.line == cl, ], max.cell.count = max.count)
        fn <- paste(expt.date, cl, "growth curves.pdf")
        if (!toFile) 
            dev.new(width = 11, height = 8.5)
        if (toFile) 
            pdf(file = fn, width = 11, height = 8.5)
        par(mfrow = c(3, 5), mar = c(2, 3, 1.5, 0.5))
        for (dr in unique(d$drug1)) {
            if (dr == "control") 
                next
            dtp <- d[d$drug1 == dr | d$drug1 == "control", ]
            if (plotFit) {
                try(plotGC_DIPfit2(data.frame(x = dtp$time, y = dtp$cell.count, 
                  uid = dtp$well, rep = dtp$drug1.conc), main = paste(cl, 
                  dr), ylim = c(-2, 6), xlab = NA, ylab = NA))
            }
            else {
                plotGC(x = dtp$time, y = dtp$cell.count, uid = dtp$well, 
                  rep = dtp$drug1.conc, main = paste(cl, dr), 
                  ylim = c(-2, 6), xlab = NA, ylab = NA)
            }
            out <- rbind(out, d)
        }
        if (toFile) 
            dev.off()
    }
    invisible(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
