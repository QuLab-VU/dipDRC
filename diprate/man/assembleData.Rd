\name{assembleData}
\alias{assembleData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
assembleData(data.loc = getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.loc}{
%%     ~~Describe \code{data.loc} here~~
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
function (data.loc = getwd()) 
{
    file.names <- dir(data.loc, full.names = TRUE)
    map.names <- file.names[grepl("[Mm]ap", file.names)]
    file.names <- setdiff(file.names[grepl("csv", file.names)], 
        map.names)
    for (fn in file.names) {
        map.pre <- unlist(strsplit(fn, ".csv"))
        if (is.na(map.pre)) {
            message(paste("Could not find a matching plate map for: ", 
                fn))
            return(NA)
        }
        else {
            map.name <- map.names[grepl(map.pre, map.names, fixed = TRUE)]
        }
        temp <- cellCountCV(read.csv(fn, as.is = TRUE))
        temp <- addMapInfo3(temp, map.name)
        ifelse(fn == file.names[1], out <- temp, out <- rbind(out, 
            temp))
    }
    out$uid <- paste(out$cell.line, out$drug1, sep = "_")
    out$nl2 <- NULL
    colnames(out) <- tolower(colnames(out))
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
