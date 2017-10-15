\name{plot.dipDRC}
\alias{plot.dipDRC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.dipDRC(drmo, plot.type = "confidence", showEC50 = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drmo}{
%%     ~~Describe \code{drmo} here~~
}
  \item{plot.type}{
%%     ~~Describe \code{plot.type} here~~
}
  \item{showEC50}{
%%     ~~Describe \code{showEC50} here~~
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
function (drmo, plot.type = "confidence", showEC50 = TRUE, ...) 
{
    if (is.na(drmo[1])) {
        message("DRM object not available to plot")
        return(NA)
    }
    param <- list(...)
    cell.line <- unique(drmo$origData$cell.line)[1]
    drug <- unique(drmo$origData$drug1)[1]
    if (!("main" \%in\% names(param))) {
        tit <- paste(cell.line, drug, sep = "_")
        param[["main"]] <- tit
    }
    if (!("type" \%in\% names(param))) 
        param[["type"]] <- plot.type
    myargs <- append(list(x = drmo), param)
    do.call(plot, args = myargs)
    if (showEC50) 
        abline(v = ED(drmo, 50, interval = "delta", display = FALSE)[1], 
            col = "red")
    abline(h = 0, col = grey(0.5), lty = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
