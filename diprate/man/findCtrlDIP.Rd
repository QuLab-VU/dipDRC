\name{findCtrlDIP}
\alias{findCtrlDIP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findCtrlDIP(times, counts, ids, col.names = c("time", "cell.counts", "well"), type = "mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{counts}{
%%     ~~Describe \code{counts} here~~
}
  \item{ids}{
%%     ~~Describe \code{ids} here~~
}
  \item{col.names}{
%%     ~~Describe \code{col.names} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (times, counts, ids, col.names = c("time", "cell.counts", 
    "well"), type = "mean") 
{
    dip <- numeric()
    dip.temp <- numeric()
    if (type == "mean") {
        m <- list()
        for (i in unique(ids)) {
            x <- times[ids == i]
            y <- log2(counts[ids == i])
            m[[i]] <- lm(y ~ x)
            dip.temp <- append(dip.temp, coef(m[[i]])[2])
            names(dip.temp)[length(names(dip.temp))] <- i
        }
        dip <- mean(dip.temp)
        dip <- append(dip, sd(dip.temp)/sqrt(length(dip.temp)))
        names(dip) <- c("mean.dip", "std.err")
        count.t0 <- floor(mean(sapply(unique(ids), FUN = function(x) head(counts[ids == 
            x], 1))))
        dts.times <- unique(times)
        dts.counts <- c(count.t0, floor(count.t0 * 2^(dip[1] * 
            dts.times[-1])))
    }
    if (type == "max") {
        n <- nEach(ids)
        n.max <- max(n)
        max.ids <- names(n)[n == n.max]
        dts.times <- times[ids \%in\% max.ids]
        dts.counts <- counts[ids \%in\% max.ids]
        a <- aggregate(dts.counts ~ dts.times, FUN = sum)
        m <- lm(log2(dts.counts) ~ dts.times, data = a)
        dip <- coef(m)[2]
        dip <- append(dip, summary(m)$coefficients[2, 2])
        names(dip) <- c("max.dip", "std.err")
        dts.times <- a$dts.times
        dts.counts <- a$dts.counts
    }
    if (type == "sum") {
        n <- nEach(ids)
        n.min <- min(n)
        dts.times <- as.numeric(sapply(unique(ids), FUN = function(x) head(times[ids == 
            x], n.min)))
        dts.counts <- as.numeric(sapply(unique(ids), FUN = function(x) head(counts[ids == 
            x], n.min)))
        a <- aggregate(dts.counts ~ dts.times, FUN = sum)
        m <- lm(log2(dts.counts) ~ dts.times, data = a)
        dip <- coef(m)[2]
        dip <- append(dip, summary(m)$coefficients[2, 2])
        names(dip) <- c("sum.dip", "std.err")
        dts.times <- a$dts.times
        dts.counts <- a$dts.counts
    }
    out.data <- data.frame(dts.times, dts.counts)
    colnames(out.data) <- col.names[1:2]
    list(data = out.data, dip = dip, model = m)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
