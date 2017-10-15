\name{plotGC_DIPfit}
\alias{plotGC_DIPfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotGC_DIPfit(dtp, tit = "unknown", toFile = FALSE, newDev = TRUE, add.line.met = "none", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtp}{
%%     ~~Describe \code{dtp} here~~
}
  \item{tit}{
%%     ~~Describe \code{tit} here~~
}
  \item{toFile}{
%%     ~~Describe \code{toFile} here~~
}
  \item{newDev}{
%%     ~~Describe \code{newDev} here~~
}
  \item{add.line.met}{
%%     ~~Describe \code{add.line.met} here~~
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
function (dtp, tit = "unknown", toFile = FALSE, newDev = TRUE, 
    add.line.met = "none", ...) 
{
    stuff <- list(...)
    if ("o" \%in\% names(stuff) & tit == "rmse") 
        tit <- paste(tit, "with o =", stuff[["o"]])
    dip <- findDIP(dtp, ...)
    add.line <- FALSE
    if (add.line.met != "none") {
        add.line <- TRUE
        dip2 <- findDIP(dtp, met = add.line.met)
    }
    fn <- paste(tit, nrow(dtp), "points.pdf")
    if ("metric" \%in\% names(stuff)) 
        tit <- paste(tit, stuff[["metric"]])
    if (newDev & !toFile) 
        dev.new(width = 7.5, height = 3)
    if (newDev & toFile) 
        pdf(file = fn, width = 7.5, height = 3)
    if (newDev) 
        par(mfrow = c(1, 3), oma = c(0, 0, 1, 0))
    plot(dip$data, main = NA, xlab = NA, ylab = NA)
    mtext(side = 1, "Time (h)", font = 2, line = 2)
    mtext(side = 2, "log2(cell number)", font = 2, line = 2)
    dip.val <- round(dip$dip, 4)
    dip.95conf <- round(abs(dip.val - confint(dip$best.model)[2, 
        1]), 5)
    legend("bottomright", c(paste("DIP =", dip.val), paste0("  \302\261", 
        dip.95conf), paste0("start =", round(dip$start.time, 
        1))), bty = "n", pch = "")
    curve(coef(dip$best.model)[1] + coef(dip$best.model)[2] * 
        x, from = 0, to = 150, add = TRUE, col = "red", lwd = 3)
    try(polygon(x = c(dip$start.time, 150, 150), y = c(coef(dip$best.model)[1] + 
        coef(dip$best.model)[2] * dip$start.time, coef(dip$best.model)[1] + 
        confint(dip$best.model)[2, 1] * 150, coef(dip$best.model)[1] + 
        confint(dip$best.model)[2, 2] * 150), col = adjustcolor("gray", 
        alpha.f = 0.4), density = NA))
    abline(v = dip$start.time, lty = 2)
    if (add.line) 
        abline(v = dip2$start.time, lty = 3, col = "red")
    plot(dip$eval.times, dip$ar2, ylim = c(0, 1), xlab = NA, 
        ylab = NA, main = NA)
    mtext(side = 1, "Time (h)", font = 2, line = 2)
    mtext(side = 2, "adj R2", font = 2, line = 2)
    abline(v = dip$start.time, lty = 2)
    if (add.line) 
        abline(v = dip2$start.time, lty = 3, col = "red")
    plot(dip$eval.times, dip$rmse, ylim = c(0, 0.5), xlab = NA, 
        ylab = NA, main = NA)
    mtext(side = 1, "Time (h)", font = 2, line = 2)
    mtext(side = 2, "RMSE", font = 2, line = 2)
    abline(v = dip$start.time, lty = 2)
    if (add.line) 
        abline(v = dip2$start.time, lty = 3, col = "red")
    r1d <- coef(fit_p5(data.frame(as.numeric(dip$eval.times), 
        dip$rmse))$m)
    curve(p5(x, r1d[1], r1d[2], r1d[3], r1d[4], r1d[5], r1d[6]), 
        from = 0, to = 120, col = "blue", lwd = 3, add = TRUE)
    if (newDev) 
        mtext(tit, outer = TRUE, side = 3, font = 2, line = -1.5)
    if (newDev & toFile) 
        dev.off()
    invisible(dip)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
