\name{getWellRates}
\alias{getWellRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getWellRates(raw, time.range = c(70, 120))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{
%%     ~~Describe \code{raw} here~~
}
  \item{time.range}{
%%     ~~Describe \code{time.range} here~~
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
function (raw, time.range = c(70, 120)) 
{
    timeName <- colnames(raw)[grep("[Tt]ime", colnames(raw))]
    wellName <- colnames(raw)[grep("[Ww]ell", colnames(raw))]
    dateName <- colnames(raw)[grep("[Dd]ate", colnames(raw))]
    if (length(wellName) > 1) 
        wellName <- wellName[nchar(wellName) == 4]
    f <- formula(paste("nl2 ~ ", timeName, " * ", wellName))
    m <- lm(f, data = raw[raw[, timeName] > time.range[1] & raw[, 
        timeName] < time.range[2], ])
    wells <- unique(raw[, wellName])
    rates <- coef(m)[grep(timeName, names(coef(m)))]
    rates <- c(rates[1], rates[-1] + rates[1])
    cl <- unique(raw$cellLine)
    expt <- ifelse(is.null(unique(raw[, dateName])), "unknown date", 
        unique(raw[, dateName]))
    out <- data.frame(Well = wells, DIP = rates, cellLine = cl, 
        Date = expt)
    rownames(out) <- NULL
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
