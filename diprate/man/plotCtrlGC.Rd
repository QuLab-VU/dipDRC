\name{plotCtrlGC}
\alias{plotCtrlGC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotCtrlGC(dat = d, type = c("raw", "lin")[1], uniq = "plate.name", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{uniq}{
%%     ~~Describe \code{uniq} here~~
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
function (dat = d, type = c("raw", "lin")[1], uniq = "plate.name", 
    ...) 
{
    prepMultiGCplot()
    for (u in unique(dat[, uniq])) {
        dtp <- d[d[, uniq] == u & d$drug1 == "control", ]
        dtp <- dtp[order(dtp$well, dtp$time), ]
        switch(type, lin = controlQC(dtp$time, dtp$cell.count, 
            dtp$well, c("time", "cell.count", "well"), cell.line.name = pn, 
            ...), raw = plotGC(dtp$time, dtp$cell.count, dtp$uid, 
            dtp$uid, ...))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
