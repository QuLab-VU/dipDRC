\name{assembleCountData}
\alias{assembleCountData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
assembleCountData(filepaths)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepaths}{
%%     ~~Describe \code{filepaths} here~~
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
function (filepaths) 
{
    readFiles <- function(fp) {
        out <- read.csv(fp, as.is = TRUE)
        out[, !grepl("X", colnames(out))]
    }
    out <- sapply(filepaths, readFiles, simplify = FALSE)
    if (!all(sapply(out, ncol) == ncol(out[[1]]))) {
        cn <- colnames(out[sapply(out, ncol) == max(sapply(out, 
            ncol))][[1]])
        temp1 <- do.call(rbind, out[sapply(out, ncol) != length(cn)])
        temp2 <- do.call(rbind, out[sapply(out, ncol) == length(cn)])
        newcn <- setdiff(cn, colnames(temp1))
        temp3 <- data.frame(matrix(ncol = length(newcn), nrow = nrow(temp1)))
        colnames(temp3) <- newcn
        temp1 <- cbind(temp1, temp3)
        out <- rbind(temp1, temp2)
        out <- out[, cn]
    }
    else {
        out <- do.call(rbind, out)
    }
    rownames(out) <- NULL
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
