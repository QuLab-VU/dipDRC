\name{effectAtMeanIC50}
\alias{effectAtMeanIC50}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
effectAtMeanIC50(drmlist = drc, drugname = "erlotinib", sig = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drmlist}{
%%     ~~Describe \code{drmlist} here~~
}
  \item{drugname}{
%%     ~~Describe \code{drugname} here~~
}
  \item{sig}{
%%     ~~Describe \code{sig} here~~
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
function (drmlist = drc, drugname = "erlotinib", sig = 3) 
{
    drms <- drmByDrug(drmlist, drug_name = drugname)
    ic50s <- sapply(drms, function(x) ED(x, coef(x)["d:(Intercept)"]/2, 
        type = "absolute", interval = "delta", display = FALSE)[1])
    conc <- signif(10^mean(log10(ic50s), na.rm = TRUE), 3)
    effect <- sapply(drms, function(x) signif(PR(x, conc)/PR(x, 
        0), sig))
    cl <- sapply(names(drms), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
    out <- data.frame(cell.line = cl, drug1 = drugname, mean.IC50 = conc, 
        rEffect.at.meanIC50 = effect)
    rownames(out) <- NULL
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
