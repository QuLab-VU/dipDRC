\name{filterCtrlData}
\alias{filterCtrlData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filterCtrlData(times, counts, ids, min.ar2 = 0.99, verbose = FALSE)
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
  \item{min.ar2}{
%%     ~~Describe \code{min.ar2} here~~
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
function (times, counts, ids, min.ar2 = 0.99, verbose = FALSE) 
{
    dip <- numeric()
    dtk.times <- integer()
    dtk.counts <- integer()
    dtk.ids <- character()
    for (i in unique(ids)) {
        if (verbose) 
            message(paste("Processing", i))
        ttf <- times[ids == i]
        ctf <- log2(counts[ids == i])
        if (length(ttf) < 5) {
            message(paste("Need at least 5 data points for ", 
                i))
            next
        }
        mintimept <- ifelse(length(ttf)/2 > 5, floor(length(ttf)/2), 
            5)
        m.ad <- lm(ctf ~ ttf)
        ar2 <- summary(m.ad)$adj.r.squared
        for (l in length(ttf):5) {
            x <- ttf[seq(l)]
            y <- ctf[seq(l)]
            if (!is.na(summary(lm(y ~ x))$adj.r.squared)) {
                pass <- summary(lm(y ~ x))$adj.r.squared >= min.ar2
            }
            else {
                pass <- FALSE
            }
            if (pass) {
                x <- head(x, -1)
                y <- head(y, -1)
                break
            }
        }
        if (!pass) {
            if (verbose) 
                message(paste("The first 5 pts of ids =", i, 
                  "is not linear within ar2 =", min.ar2))
            next
        }
        else {
            dtk.times <- append(dtk.times, x)
            dtk.counts <- append(dtk.counts, floor(2^y))
            dtk.ids <- append(dtk.ids, rep(i, length(x)))
        }
    }
    data.frame(times = dtk.times, counts = dtk.counts, ids = dtk.ids)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
