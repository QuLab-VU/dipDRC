\name{makeDil}
\alias{makeDil}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeDil(max.conc = 4000, DF = 4, n = 8, sig = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{max.conc}{
%%     ~~Describe \code{max.conc} here~~
}
  \item{DF}{
%%     ~~Describe \code{DF} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{sig}{
%%     ~~Describe \code{sig} here~~
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
function (max.conc = 4000, DF = 4, n = 8, sig = 2) 
{
    signif(DF^(seq(from = log(max.conc)/log(DF), to = log(max.conc)/log(DF) - 
        (n - 1), by = -1)), sig)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
