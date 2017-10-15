\name{addMapInfo}
\alias{addMapInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
addMapInfo(dfa, path.to.map)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfa}{
%%     ~~Describe \code{dfa} here~~
}
  \item{path.to.map}{
%%     ~~Describe \code{path.to.map} here~~
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
function (dfa, path.to.map) 
{
    map <- getMapInfo(path.to.map)
    cn <- colnames(map)
    if (any(grepl("[Dd]ate", cn))) {
        date <- unique(map[, grep("[Dd]ate", cn)[1]])
    }
    else date <- substr(path.to.map, 1, 10)
    date <- gsub("-", "", date)
    date <- date[!grepl("null", date)]
    wellName <- cn[grep("[Ww]ell", cn)]
    wellNameDFA <- colnames(dfa)[grep("[Ww]ell", colnames(dfa))]
    if (length(wellName) > 1) 
        wellName <- wellName[nchar(wellName) == 4]
    dfa$expt.date <- date
    dfa$cell.line <- map[match(dfa[, wellNameDFA], map[, wellName]), 
        grep("[cC]ell.line", cn)]
    dfa$drug1 <- map[match(dfa[, wellNameDFA], map[, wellName]), 
        grep("[dD]rug1", cn)]
    dfa$drug1.conc <- map[match(dfa[, wellNameDFA], map[, wellName]), 
        grepl("[cC]onc1", cn) & !grepl("units", cn)]
    dfa$drug1.units <- map[match(dfa[, wellNameDFA], map[, wellName]), 
        grepl("[cC]onc1.units", cn)]
    dfa$drug2 <- map[match(dfa[, wellNameDFA], map[, wellName]), 
        grep("[dD]rug2", cn)]
    dfa$drug2.conc <- map[match(dfa[, wellNameDFA], map[, wellName]), 
        grepl("[cC]onc2", cn) & !grepl("units", cn)]
    dfa$drug2.units <- map[match(dfa[, wellNameDFA], map[, wellName]), 
        grepl("[cC]onc2.units", cn)]
    dfa
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
