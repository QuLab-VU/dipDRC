\name{addLL4curve}
\alias{addLL4curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
addLL4curve(drmodel, fromval = 1e-12, toval = 1e-05, norm = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drmodel}{
%%     ~~Describe \code{drmodel} here~~
}
  \item{fromval}{
%%     ~~Describe \code{fromval} here~~
}
  \item{toval}{
%%     ~~Describe \code{toval} here~~
}
  \item{norm}{
%%     ~~Describe \code{norm} here~~
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
function (drmodel, fromval = 1e-12, toval = 1e-05, norm = FALSE, 
    ...) 
{
    param <- coef(drmodel)
    names(param) <- letters[2:5]
    if (norm) {
        param[2] <- param[2]/param[3]
        param[3] <- 1
        if (param[2] > param[3]) {
            temp <- param[2]
            param[2] <- param[3]
            param[3] <- temp
        }
    }
    curve(do.call(myll4, args = append(list(x), as.list(param))), 
        from = fromval, to = toval, add = TRUE, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
