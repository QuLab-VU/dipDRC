\name{plotGC}
\alias{plotGC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotGC(x, y, uid, rep = uid, y.type = "count", color = TRUE, leg = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{uid}{
%%     ~~Describe \code{uid} here~~
}
  \item{rep}{
%%     ~~Describe \code{rep} here~~
}
  \item{y.type}{
%%     ~~Describe \code{y.type} here~~
}
  \item{color}{
%%     ~~Describe \code{color} here~~
}
  \item{leg}{
%%     ~~Describe \code{leg} here~~
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
function (x, y, uid, rep = uid, y.type = "count", color = TRUE, 
    leg = TRUE, ...) 
{
    ifelse(y.type == "count", plot.y <- log2norm(y, ids = uid), 
        plot.y <- y)
    urep <- unique(rep)
    ifelse(color, rep.col <- gplots::colorpanel(n = length(urep), 
        low = "blue", mid = "orange", high = "red"), rep.col <- rep("black", 
        length(urep)))
    plot(x, plot.y, type = "n", ...)
    for (myrep in urep) for (id in unique(uid)) lines(x[uid == 
        id & rep == myrep], plot.y[uid == id & rep == myrep], 
        col = rep.col[match(myrep, urep)])
    if (leg) 
        legend("topleft", legend = urep, col = rep.col, lwd = 1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
