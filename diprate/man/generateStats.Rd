\name{generateStats}
\alias{generateStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generateStats(dat = DS.corr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
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
function (dat = DS.corr) 
{
    stats <- aggregate(val ~ subl + var, data = dat, mean)
    colnames(stats)[3] <- "mean"
    stats <- cbind(stats, aggregate(val ~ subl + var, data = dat, 
        length)["val"])
    colnames(stats)[4] <- "n"
    temp <- as.data.frame(aggregate(val ~ subl + var, data = dat, 
        FUN = function(x) as.numeric(t.test(x)$conf.int))[, 3])
    dev <- (temp[, 2] - temp[, 1])/2
    stats <- cbind(stats, dev)
    ec50.mean <- sapply(allDRM, FUN = function(x) ED(x, 50, interval = "delta", 
        display = FALSE)[1])
    ec50.dev <- (sapply(allDRM, FUN = function(x) ED(x, 50, interval = "delta", 
        display = FALSE)[3]) - sapply(allDRM, FUN = function(x) ED(x, 
        50, interval = "delta", display = FALSE)[2]))/2
    ec50.n <- sapply(allDRM, FUN = function(x) length(unique(x$origData$exptDate)))
    ec50 <- cbind(mean = ec50.mean, dev = ec50.dev, n = ec50.n)
    ec50[, 1:2] <- apply(ec50[, 1:2], c(1, 2), FUN = function(x) signif(x * 
        1e+09, 3))
    stats <- rbind(stats, data.frame(subl = rownames(ec50), var = "EC50", 
        ec50))
    stats <- stats[order(stats$var, stats$subl), ]
    stats[, 3:4] <- apply(stats[, 3:4], c(1, 2), FUN = function(x) signif(x, 
        digits = 3))
    stats <- stats[!stats$subl \%in\% c("DS8", "DS9", "PC9"), ]
    stats$dev <- sapply(stats$dev, FUN = function(x) signif(x, 
        3))
    rm(list = c("ec50.mean", "ec50.dev", "ec50.n"))
    type <- as.character(unique(stats$var))
    stats
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
