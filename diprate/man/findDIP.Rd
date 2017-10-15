\name{findDIP}
\alias{findDIP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findDIP(dtf, name = "unknown", all.models = FALSE, metric = c("opt", "ar2", "rmse")[1], o = 0, dat.type = "cell.counts", print.dip = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtf}{
%%     ~~Describe \code{dtf} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{all.models}{
%%     ~~Describe \code{all.models} here~~
}
  \item{metric}{
%%     ~~Describe \code{metric} here~~
}
  \item{o}{
%%     ~~Describe \code{o} here~~
}
  \item{dat.type}{
%%     ~~Describe \code{dat.type} here~~
}
  \item{print.dip}{
%%     ~~Describe \code{print.dip} here~~
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
function (dtf, name = "unknown", all.models = FALSE, metric = c("opt", 
    "ar2", "rmse")[1], o = 0, dat.type = "cell.counts", print.dip = FALSE) 
{
    if (nrow(dtf) < 4 | ncol(dtf) < 2) {
        message("findDIP requires data frame with nrows >= 4 and ncol >= 2")
        return(NA)
    }
    x <- dtf[, 1]
    if (dat.type == "cell.counts") {
        y <- log2(dtf[, 2])
    }
    else {
        y <- dtf[, 2]
    }
    n <- nrow(dtf) - 2
    m <- list()
    rmse <- numeric()
    ar2 <- numeric()
    p <- numeric()
    b <- y
    a <- x
    i <- 1
    for (i in seq(n)) {
        m[[i]] <- lm(b ~ a)
        rmse <- append(rmse, rmsd(m[[i]]$residuals))
        ar2 <- append(ar2, summary(m[[i]])$adj.r.squared)
        p <- append(p, summary(m[[i]])$coefficients[2, 4])
        a <- a[-1]
        b <- b[-1]
    }
    eval.times <- x[seq(n)]
    rmse.p5 <- tryCatch({
        fit_p5(data.frame(x = eval.times, rmse = rmse))
    }, error = function(cond) {
        NA
    })
    rmse.p5.coef <- tryCatch({
        coef(rmse.p5$m)
    }, error = function(cond) {
        NA
    })
    rmse.p5.1std.coef <- tryCatch({
        rmse.p5$coef.1stderiv
    }, error = function(cond) {
        NA
    })
    f <- function(..., offset) p5(...) + offset
    opt <- ar2 * (1 - rmse)^2 * (length(eval.times) - seq(eval.times))^0.25
    idx <- switch(metric, ar2 = ifelse(n <= 5, match(max(ar2), 
        ar2), match(max(ar2[seq(floor(n * 0.9))]), ar2)), rmse = tryCatch({
        floor(uniroot(f, interval = c(0, 72), tol = 1e-06, extendInt = "upX", 
            offset = o, int = rmse.p5.1std.coef[1], b1 = rmse.p5.1std.coef[2], 
            b2 = rmse.p5.1std.coef[3], b3 = rmse.p5.1std.coef[4], 
            b4 = rmse.p5.1std.coef[5], b5 = rmse.p5.1std.coef[6])$root)
    }, error = function(cond) {
        message("Could not find root of 1st deriv of RMSE")
        message(cond)
        1
    }), opt = match(max(opt), opt))
    dip <- coef(m[[idx]])[2]
    ci <- diff(confint(m[[idx]])[2, ])/2
    names(ci) <- "\302\261"
    if (print.dip) 
        print(paste("DIP =", round(dip, 4), "starting at", round(x[idx], 
            2), "with", n + 2, "data points"))
    out = list(data = data.frame(Time_h = x, l2 = y), model = m, 
        metric.used = metric, n = n, idx = idx, best.model = m[[idx]], 
        opt = opt, eval.times = x[seq(n)], rmse = rmse, ar2 = ar2, 
        p = p, start.time = x[idx], dip = dip, dip.95ci = ci, 
        rmse.p5.1std.coef = rmse.p5.1std.coef)
    if (!all.models) {
        out[["model"]] <- NULL
    }
    invisible(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
