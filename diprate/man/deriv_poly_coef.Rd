\name{deriv_poly_coef}
\alias{deriv_poly_coef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
deriv_poly_coef(co)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{co}{
%%     ~~Describe \code{co} here~~
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
function (co) 
{
    stopifnot(all(sapply(names(co), FUN = function(x) x \%in\% 
        c("int", "b1", "b2", "b3", "b4", "b5"))))
    new.co <- c(co["b1"], 2 * co["b2"], 3 * co["b3"], 4 * co["b4"], 
        5 * co["b5"], 0)
    new.co[is.na(new.co)] <- 0
    names(new.co) <- c("int", "b1", "b2", "b3", "b4", "b5")
    new.co
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
