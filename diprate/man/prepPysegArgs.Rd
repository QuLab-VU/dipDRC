\name{prepPysegArgs}
\alias{prepPysegArgs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
prepPysegArgs(finfo, rootdir, expdir = "", ch2_name = "FITC", savefi = FALSE, momlogpath = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{finfo}{
%%     ~~Describe \code{finfo} here~~
}
  \item{rootdir}{
%%     ~~Describe \code{rootdir} here~~
}
  \item{expdir}{
%%     ~~Describe \code{expdir} here~~
}
  \item{ch2_name}{
%%     ~~Describe \code{ch2_name} here~~
}
  \item{savefi}{
%%     ~~Describe \code{savefi} here~~
}
  \item{momlogpath}{
%%     ~~Describe \code{momlogpath} here~~
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
function (finfo, rootdir, expdir = "", ch2_name = "FITC", savefi = FALSE, 
    momlogpath = "", ...) 
{
    if (finfo == "new") {
        source("getMXfileInfo.r")
        message("generating imageFileInfo file")
        fi <- getMXfileInfo(rootdir)
        if (momlogpath != "" & dir.exists(momlogpath) & length(dir(momlogpath)) > 
            0) 
            fi <- tryCatch({
                addPlateNames(fi, momlogpath)
            }, error = function(cond) {
                message("could not add plate_names")
                return(fi)
            })
        if (savefi) {
            bn <- basename(rootdir)
            write.csv(fi, file = paste(rootdir, "/", bn, "imageFileInfo.csv", 
                sep = ""), row.names = FALSE)
        }
    }
    else {
        if (file.exists(normalizePath(finfo))) {
            fi <- read.csv(finfo, as.is = TRUE)
        }
        else {
            message("could not find fileInfo")
            return(NA)
        }
    }
    makeTaskArgs(fi, root_dir = rootdir, defargs = as.list(data.frame(ch2_name = ch2_name, 
        expdir = expdir)), ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
