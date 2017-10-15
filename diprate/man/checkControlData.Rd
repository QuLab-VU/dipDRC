\name{checkControlData}
\alias{checkControlData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
checkControlData(ctrl = ad[ad$drug1 == "control", ], plotIt = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ctrl}{
%%     ~~Describe \code{ctrl} here~~
}
  \item{plotIt}{
%%     ~~Describe \code{plotIt} here~~
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
function (ctrl = ad[ad$drug1 == "control", ], plotIt = TRUE) 
{
    ctrl.fd <- data.frame()
    for (cl in unique(ctrl$cell.line)) {
        dfc <- ctrl[ctrl$cell.line == cl, ]
        temp <- controlQC(dfc$time, dfc$cell.count, dfc$well, 
            xlim = c(0, 150), ylim = c(0, 7), cell.line.name = cl, 
            plotIt = plotIt)
        temp$cell.line <- cl
        ctrl.fd <- rbind(ctrl.fd, temp)
    }
    ctrl.fd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
