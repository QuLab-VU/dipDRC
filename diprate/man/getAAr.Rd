\name{getAAr}
\alias{getAAr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getAAr(p, minmax = c(1e-12, 1e-05), RespRatio = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{minmax}{
%%     ~~Describe \code{minmax} here~~
}
  \item{RespRatio}{
%%     ~~Describe \code{RespRatio} here~~
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
function (p, minmax = c(1e-12, 1e-05), RespRatio = TRUE) 
{
    if (class(p) == "drc") 
        p <- coef(p)
    if (!(class(p) == "numeric" & length(p) == 4)) {
        message("argument to getAAr must be drm object or its 4 coefficients (assuming LL.4)")
        return(NA)
    }
    names(p) <- c("h", "Emax", "E0", "EC50")
    out <- 1/p["h"] * log10((p["EC50"]^p["h"] + minmax[2]^p["h"])/(p["EC50"]^p["h"] + 
        minmax[1]^p["h"]))
    if (RespRatio) 
        out <- out * (p["E0"] - p["Emax"])/p["E0"]
    names(out) <- "AAr"
    out[out < 0] <- 0
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
