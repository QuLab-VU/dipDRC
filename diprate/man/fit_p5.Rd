\name{fit_p5}
\alias{fit_p5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fit_p5(dtf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtf}{
%%     ~~Describe \code{dtf} here~~
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
function (dtf) 
{
    x <- colnames(dtf)[1]
    y <- colnames(dtf)[2]
    form <- formula(paste(y, "~ p5(", x, ", int, b1, b2, b3, b4, b5)"))
    m.p5 <- nls(form, data = dtf, start = list(int = 1, b1 = 1, 
        b2 = 1, b3 = 1, b4 = 1, b5 = 1))
    m.p5.coef.1stderiv <- deriv_poly_coef(coef(m.p5))
    m.p5.coef.2ndderiv <- deriv_poly_coef(m.p5.coef.1stderiv)
    list(m = m.p5, coef.1stderiv = m.p5.coef.1stderiv, coef.2ndderiv = m.p5.coef.2ndderiv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
