\name{drcDIP}
\alias{drcDIP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
drcDIP(dtf, cl = cl, drg = drg, PIP.range = PIP.range, norm = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtf}{
%%     ~~Describe \code{dtf} here~~
}
  \item{cl}{
%%     ~~Describe \code{cl} here~~
}
  \item{drg}{
%%     ~~Describe \code{drg} here~~
}
  \item{PIP.range}{
%%     ~~Describe \code{PIP.range} here~~
}
  \item{norm}{
%%     ~~Describe \code{norm} here~~
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
function (dtf, cl = cl, drg = drg, PIP.range = PIP.range, norm = T) 
{
    rates <- data.frame(cellLine = character(), conc = integer(), 
        DIPrate = numeric(), drug = character())
    dtp <- dtf[dtf$cellLine == cl & dtf$drug == drg, ]
    Uconc <- unique(dtp$conc)
    for (co in 2:length(Uconc)) {
        rates <- rbind(rates, data.frame(cellLine = cl, conc = Uconc[co], 
            DIPrate = getDIP(dtp[dtp$conc == Uconc[co] & dtp$time >= 
                PIP.range[1] & dtp$time <= PIP.range[2], ]), 
            drug = drg))
    }
    rates$normDIP <- rates$DIPrate/rates[rates$conc == min(rates$conc), 
        "DIPrate"]
    if (norm) {
        drm(normDIP ~ conc, data = rates, fct = LL.4())
    }
    else {
        drm(DIPrate ~ conc, data = rates, fct = LL.4())
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
