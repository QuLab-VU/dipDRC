\name{utMAlin}
\alias{utMAlin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
utMAlin(values, times, tau = 24)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{values}{
%%     ~~Describe \code{values} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
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
function (values, times, tau = 24) 
{
    if (class(times) == "difftime") 
        times <- as.numeric(times)
    out <- rep(0, length(values))
    out[1] <- values[1]
    left <- 1
    roll_area <- left_area <- values[1] * tau
    for (right in 2:length(values)) {
        roll_area <- roll_area + (values[right - 1] + values[right])/2 * 
            (times[right] - times[right - 1])
        roll_area = roll_area - left_area
        t_left_new = times[right] - tau
        while (times[left] <= t_left_new) {
            roll_area = roll_area - (values[left] + values[left + 
                1])/2 * (times[left + 1] - times[left])
            left = left + 1
        }
        left_area <- trapezoid(times[max(1, left - 1)], t_left_new, 
            times[left], values[max(1, left - 1)], values[left])
        roll_area <- roll_area + left_area
        out[right] <- roll_area/tau
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
