\name{plotAllDrugs}
\alias{plotAllDrugs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotAllDrugs(dat, drug.col = "drug1", uid = "well", output = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{drug.col}{
%%     ~~Describe \code{drug.col} here~~
}
  \item{uid}{
%%     ~~Describe \code{uid} here~~
}
  \item{output}{
%%     ~~Describe \code{output} here~~
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
function (dat, drug.col = "drug1", uid = "well", output = FALSE, 
    ...) 
{
    cl <- unique(dat$cell.line)
    if (length(cl) > 1) {
        message("prepDRCplot needs data.frame containing single cell line")
        return(invisible(NA))
    }
    out <- list()
    drugs <- sort(unique(dat[, drug.col]))
    drugs <- setdiff(drugs, "control")
    for (dr in drugs) {
        df.name <- paste(cl, dr, sep = "_")
        dtp <- dat[dat[, drug.col] == dr | dat[, drug.col] == 
            "control", ]
        dtp <- dtp[order(dtp$drug1.conc, dtp[, uid], dtp$time), 
            ]
        plotGC(dtp$time, dtp$cell.count, dtp[, uid], signif(dtp$drug1.conc, 
            2), main = paste(cl, dr), ...)
        out[[df.name]] <- dtp
    }
    if (output) 
        invisible(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
