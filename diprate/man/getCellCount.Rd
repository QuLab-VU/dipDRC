\name{getCellCount}
\alias{getCellCount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCellCount(parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
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
function (parameters) 
{
    for (dc in drug.conc) {
        p <- append(parameters, dc)
        names(p)[5] <- "Drug"
        temp <- as.data.frame(ode(y = states, times = times, 
            func = twoStateModel, parms = p))
        Cell.tot <- as.numeric(temp$Cell + temp$CellPrime)
        ifelse(dc == drug.conc[1], out <- data.frame(times, Cell.tot), 
            out <- cbind(out, Cell.tot))
    }
    colnames(out) <- c("time", as.character(drug.conc))
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
