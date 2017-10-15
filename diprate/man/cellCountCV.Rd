\name{cellCountCV}
\alias{cellCountCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cellCountCV(rawCVdata, normPos = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawCVdata}{
%%     ~~Describe \code{rawCVdata} here~~
}
  \item{normPos}{
%%     ~~Describe \code{normPos} here~~
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
function (rawCVdata, normPos = 1) 
{
    log2norm <- function(count, ids, norm_idx = 1, zero = log2(0.999)) {
        l2 <- log2(count)
        l2[is.infinite(l2)] <- zero
        norm <- numeric()
        group <- as.character(unique(ids))
        for (i in group) {
            d <- l2[ids == i]
            norm <- append(norm, d - d[norm_idx])
        }
        norm
    }
    d <- rawCVdata
    d <- d[, c("Column", "Row", "Cell.Nucleus")]
    colnames(d)[3] <- "Cell.count"
    d$Row <- as.character(d$Row)
    d$Column <- as.character(d$Column)
    idx <- grep("Timespan", d[, 1])
    times <- as.numeric(d[, 2][idx])
    time.idx <- seq(length(times))
    numWells <- idx[2] - idx[1] - 2
    for (tp in time.idx) {
        temp <- d[(idx[tp] + 1):(idx[tp] + numWells), ]
        temp$Time <- times[tp]
        ifelse(tp == time.idx[1], a <- temp, a <- rbind(a, temp))
    }
    rownames(a) <- seq(nrow(a))
    a[nchar(a$Column) == 1, "Column"] <- paste0("0", a[nchar(a$Column) == 
        1, "Column"])
    a$Well <- paste0(a$Row, a$Column)
    a$Column <- as.integer(a$Column)
    a <- a[order(a$Well), ]
    a <- a[!is.na(a$Cell.count), ]
    a$nl2 <- log2norm(a$Cell.count, a$Well, norm_idx = normPos)
    rownames(a) <- NULL
    a
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
