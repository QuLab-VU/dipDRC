\name{parseMBK}
\alias{parseMBK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
parseMBK(mbk_file_path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mbk_file_path}{
%%     ~~Describe \code{mbk_file_path} here~~
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
function (mbk_file_path) 
{
    fn <- mbk_file_path
    h <- readLines(fn, 4, warn = FALSE)
    m <- read.csv(fn, sep = "\t", skip = 4, header = FALSE, as.is = TRUE)
    cn <- unlist(strsplit(h[3], ","))
    cn[1] <- unlist(strsplit(cn[1], "[Columns]=\"", fixed = TRUE))[2]
    cn[length(cn)] <- unlist(strsplit(tail(cn, 1), "\"", fixed = TRUE))
    colnames(m) <- cn
    m[, sapply(m[1, ], class) == "character"] <- apply(m[, sapply(m[1, 
        ], class) == "character"], 2, function(x) gsub("\\\\", 
        "|", x))
    m
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
