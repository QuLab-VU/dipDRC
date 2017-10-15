\name{plotRaw}
\alias{plotRaw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotRaw(dat = ad, id.name = "cell.line", toFile = FALSE, fn = "RawGC", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{id.name}{
%%     ~~Describe \code{id.name} here~~
}
  \item{toFile}{
%%     ~~Describe \code{toFile} here~~
}
  \item{fn}{
%%     ~~Describe \code{fn} here~~
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
function (dat = ad, id.name = "cell.line", toFile = FALSE, fn = "RawGC", 
    ...) 
{
    for (id in unique(dat[, id.name])) {
        dtp <- dat[dat[, id.name] == id, ]
        dtp <- dtp[order(dtp$drug1.conc, dtp$time), ]
        prepMultiGCplot(toFile = toFile, fn = paste(fn, "_", 
            id, ".pdf", sep = ""))
        plotAllDrugs(dtp, ylim = c(-2, 7), ...)
        if (toFile) 
            dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
