\name{scatterPlotErr}
\alias{scatterPlotErr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scatterPlotErr(dtp, x.var, y.var, add.line = FALSE, xlab.txt = x.var, ylab.txt = y.var, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtp}{
%%     ~~Describe \code{dtp} here~~
}
  \item{x.var}{
%%     ~~Describe \code{x.var} here~~
}
  \item{y.var}{
%%     ~~Describe \code{y.var} here~~
}
  \item{add.line}{
%%     ~~Describe \code{add.line} here~~
}
  \item{xlab.txt}{
%%     ~~Describe \code{xlab.txt} here~~
}
  \item{ylab.txt}{
%%     ~~Describe \code{ylab.txt} here~~
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
function (dtp, x.var, y.var, add.line = FALSE, xlab.txt = x.var, 
    ylab.txt = y.var, ...) 
{
    x <- dtp[dtp$var == x.var, ]
    y <- dtp[dtp$var == y.var, ]
    plotCI(x$mean, y$mean, uiw = y$dev, err = "y", pch = NA, 
        ...)
    plotCI(x$mean, y$mean, uiw = x$dev, err = "x", pch = NA, 
        add = TRUE, ...)
    points(x$mean, y$mean, pch = 19)
    mtext(side = 1, xlab.txt, font = 2, line = 2)
    mtext(side = 2, ylab.txt, font = 2, line = 2)
    my.fit <- lm(y$mean ~ x$mean)
    ar2 <- round(summary(my.fit)$adj.r.squared, 2)
    mtext(as.expression(bquote(Adj ~ R^2 ~ "=" ~ .(ar2))), side = 1, 
        line = -1.25, cex = 0.75, adj = 0.95)
    if (add.line) 
        abline(coef(my.fit), lty = 2, lwd = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
