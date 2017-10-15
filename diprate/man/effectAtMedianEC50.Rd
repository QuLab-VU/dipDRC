\name{effectAtMedianEC50}
\alias{effectAtMedianEC50}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
effectAtMedianEC50(drmlist = drc, drugname = "erlotinib", sig = 3)
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
    med.ec50 <- signif(median(unlist(sapply(drms, function(x) coef(x)["e:(Intercept)"], 
        simplify = FALSE))), sig)
    effect <- sapply(drms, function(x) signif(PR(x, med.ec50)/PR(x, 
        0), sig))
    cl <- sapply(names(drms), function(x) gsub(paste0(".", drugname), 
        "", x))
    out <- data.frame(cell.line = cl, drug1 = drugname, median.EC50 = med.ec50, 
        rEffect.at.medianEC50 = effect)
    out <- out[order(out$cell.line), ]
    rownames(out) <- NULL
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
