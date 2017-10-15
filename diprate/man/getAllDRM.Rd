\name{getAllDRM}
\alias{getAllDRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getAllDRM(dat, plotIt = TRUE, toFile = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{plotIt}{
%%     ~~Describe \code{plotIt} here~~
}
  \item{toFile}{
%%     ~~Describe \code{toFile} here~~
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
function (dat, plotIt = TRUE, toFile = FALSE, ...) 
{
    all.drm <- lapply(names(dat), FUN = function(x) {
        cl <- unlist(lapply(strsplit(x, "-", fixed = TRUE), function(z) z[1]))
        drugs <- unique(dat[[x]]$drug1)[!grepl("control", unique(dat[[x]]$drug1))]
        prepMultiDRCplot(toFile = toFile, fn = paste("multiDRM_", 
            x, ".pdf", sep = ""))
        temp <- sapply(drugs, function(dr) {
            dtf <- ssCLD(dat[[x]], cl, dr)
            n <- unique(paste(cl, dr, sep = "_"))
            as.list(dipDRC(dtf, plotIt = plotIt, main = n, ...))
        }, simplify = FALSE)
        if (toFile) 
            dev.off()
        temp
    })
    names(all.drm) <- names(dat)
    invisible(all.drm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
