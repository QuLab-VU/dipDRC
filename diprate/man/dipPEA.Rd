\name{dipPEA}
\alias{dipPEA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dipPEA(k_on, k_off, DIP0, DIPmax, drug)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k_on}{
%%     ~~Describe \code{k_on} here~~
}
  \item{k_off}{
%%     ~~Describe \code{k_off} here~~
}
  \item{DIP0}{
%%     ~~Describe \code{DIP0} here~~
}
  \item{DIPmax}{
%%     ~~Describe \code{DIPmax} here~~
}
  \item{drug}{
%%     ~~Describe \code{drug} here~~
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
function (k_on, k_off, DIP0, DIPmax, drug) 
{
    (1/log(2)) * (k_off * DIP0 + (k_on * drug * DIPmax))/(k_on * 
        drug + k_off)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
