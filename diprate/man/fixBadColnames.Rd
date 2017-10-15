\name{fixBadColnames}
\alias{fixBadColnames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fixBadColnames(dat, cn = c("cell_count", "file_name", "plate_id", "well"), overwrite = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{cn}{
%%     ~~Describe \code{cn} here~~
}
  \item{overwrite}{
%%     ~~Describe \code{overwrite} here~~
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
function (dat, cn = c("cell_count", "file_name", "plate_id", 
    "well"), overwrite = FALSE) 
{
    if (class(dat) == "character") {
        file_path <- dat
        if (file.exists(file_path)) 
            temp <- read.csv(file_path, as.is = TRUE)
    }
    else if (class(dat) == "data.frame") 
        temp <- dat
    if (colnames(temp)[1] != "cell_count") {
        message("Bad colnames in ", tail(unlist(strsplit(file_path, 
            "/")), 1))
        colnames(temp) <- cn
        temp <- temp[apply(temp, 1, function(x) !any(is.na(x))), 
            ]
        temp <- temp[order(temp$well), ]
        rownames(temp) <- NULL
        if (overwrite & class(dat) == "character") 
            write.csv(temp, file = file_path, row.names = FALSE)
    }
    invisible(temp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
