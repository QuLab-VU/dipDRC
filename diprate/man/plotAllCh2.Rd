\name{plotAllCh2}
\alias{plotAllCh2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotAllCh2(datlist = pd, toFile = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datlist}{
%%     ~~Describe \code{datlist} here~~
}
  \item{toFile}{
%%     ~~Describe \code{toFile} here~~
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
function (datlist = pd, toFile = FALSE) 
{
    lapply(names(datlist), function(cl) {
        prepMultiGCplot(toFile = toFile, fn = paste("Ch2PosPlots_", 
            cl, ".pdf", sep = ""))
        drugs <- setdiff(unique(datlist[[cl]]$drug1), "control")
        lapply(drugs, function(dr) {
            dat <- ssCLD(datlist[[cl]], cl = cl, drug = dr)
            dat <- dat[order(dat$drug1.conc), ]
            if (is.null(dat$ch2posfrac)) 
                dat$ch2posfrac <- dat$ch2.pos/dat$cell.count
            ylim <- c(0, 1)
            main <- paste(cl, substr(dr, 1, 15), sep = "+")
            ylab <- expression(Fraction ~ of ~ Ch2^"+" ~ cells)
            xlab <- "Time (h)"
            myargs <- list(dat = dat, ylim = ylim, main = main, 
                xlab = xlab, ylab = ylab)
            do.call(plotFrac, myargs)
            return()
        })
        if (toFile) 
            dev.off()
        return()
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
