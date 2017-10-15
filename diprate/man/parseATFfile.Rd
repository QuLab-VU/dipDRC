\name{parseATFfile}
\alias{parseATFfile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
parseATFfile(path.to.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path.to.file}{
%%     ~~Describe \code{path.to.file} here~~
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
function (path.to.file) 
{
    fn <- path.to.file
    h <- read.delim(fn, nrows = 5, as.is = TRUE)[4, 1]
    h <- unlist(strsplit(h, "_"))
    if (length(h) == 6) {
        names(h) <- c("expt.name", "state", "date", "time", "instrument", 
            "plate.id")
    }
    else {
        names(h) <- c("expt.name", "state", "date", "time", "plate.name", 
            "instrument", "plate.id")
    }
    d <- read.csv(fn, skip = 6, sep = "\t", header = TRUE)
    d <- d[, 1:7]
    renameCN <- function(coln) {
        new.name <- switch(coln, Plate.ID = "plate.id", Well.Name = "well", 
            Run.Settings.ID = "run.settings.id", MEASUREMENT.SET.ID = "measurement.id", 
            Total.Cells..MultiWaveScoring. = "cell.count", Positive.W2..MultiWaveScoring. = "s_g2.count")
        ifelse(is.null(new.name), coln, new.name)
    }
    colnames(d) <- sapply(colnames(d), renameCN)
    pcn <- c("plate.id", "well", "run.settings.id", "measurement.id", 
        "cell.count", "s_g2.count")
    d <- d[, colnames(d) \%in\% pcn]
    d$plate.name <- h["plate.name"]
    d$image.time <- paste(h["date"], paste(substr(h["time"], 
        1, 2), substr(h["time"], 3, 4), substr(h["time"], 5, 
        6), sep = ":"))
    d$image.time <- strptime(d$image.time, format = "\%Y\%m\%d \%T")
    d
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
