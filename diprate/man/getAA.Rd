\name{getAA}
\alias{getAA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getAA(p, drugconcrange = c(1e-12, 1e-05), minval = -1, norm = TRUE, removeNE = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{drugconcrange}{
%%     ~~Describe \code{drugconcrange} here~~
}
  \item{minval}{
%%     ~~Describe \code{minval} here~~
}
  \item{norm}{
%%     ~~Describe \code{norm} here~~
}
  \item{removeNE}{
%%     ~~Describe \code{removeNE} here~~
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
function (p, drugconcrange = c(1e-12, 1e-05), minval = -1, norm = TRUE, 
    removeNE = FALSE) 
{
    if (class(p) == "drc") 
        p <- coef(p)
    if (!(class(p) == "numeric" & length(p) == 4)) {
        message("argument to getAAr must be drm object or its 4 coefficients (assuming LL.4)")
        return(NA)
    }
    names(p) <- letters[2:5]
    if (norm) {
        p[["c"]] <- p[["c"]]/p[["d"]]
        p[["d"]] <- 1
    }
    xvals <- 10^seq(log10(drugconcrange[1]), log10(drugconcrange[2]), 
        0.1)
    yvals <- do.call(myll4, args = append(list(x = xvals), as.list(p))) - 
        1
    yvals <- sapply(yvals, function(x) ifelse(x < minval, minval, 
        x))
    if (removeNE) {
        sapply(seq(length(p)), function(x) ifelse(p[["b"]][x] < 
            0, 0, (-sum(yvals)/length(yvals))[x]))
    }
    else {
        -sum(yvals)/length(yvals)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
