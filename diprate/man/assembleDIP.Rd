\name{assembleDIP}
\alias{assembleDIP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
assembleDIP(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
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
function (dat) 
{
    mean.rates <- data.frame(subline = sort(unique(dat$Subline)), 
        DIP = 0, DIP.dev = 0)
    all.rates <- data.frame()
    for (sl in sort(unique(dat$Subline))) {
        dtf <- dat[dat$Subline == sl, ]
        m <- lm(log2(Cell.count) ~ Time_h * factor(UID), data = dtf[dtf$Time_h > 
            72, ])
        dip <- coef(m)[grep("Time_h", names(coef(m)))]
        dip <- c(dip[1], dip[-1] + dip[1])
        names(dip) <- gsub("Time_h:factor\\(UID\\)", "", names(dip))
        names(dip)[1] <- setdiff(unique(dtf$UID), names(dip))
        all.rates <- rbind(all.rates, data.frame(Subline = sl, 
            DIP = dip, UID = names(dip)))
        mean.val <- t.test(dip)[["estimate"]]
        dev.val <- (t.test(dip)["conf.int"][[1]][2] - t.test(dip)["conf.int"][[1]][1])/2
        mean.rates[mean.rates$subline == sl, 2:3] <- c(mean.val, 
            dev.val)
    }
    list(all.rates = all.rates, mean.rates = mean.rates)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
