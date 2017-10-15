\name{parseImageJdata}
\alias{parseImageJdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
parseImageJdata(dataFol, fileName, mapName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataFol}{
%%     ~~Describe \code{dataFol} here~~
}
  \item{fileName}{
%%     ~~Describe \code{fileName} here~~
}
  \item{mapName}{
%%     ~~Describe \code{mapName} here~~
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
function (dataFol, fileName, mapName) 
{
    source("log2norm.R")
    xl <- FALSE
    if (grepl(".xl", fileName)) 
        xl <- TRUE
    ifelse(xl, d <- gdata::read.xls(paste0(dataFol, fileName), 
        header = FALSE), d <- read.csv(paste0(dataFol, fileName), 
        as.is = TRUE, header = FALSE))
    d <- d[, 1:2]
    d[, 1] <- as.character(d[, 1])
    colnames(d) <- c("name", "cell.count")
    d$acq.time <- strptime(substr(d$name, 1, 14), "\%Y\%m\%d\%H\%M\%S")
    well.temp <- unlist(strsplit(d$name, "-"))
    row.temp <- LETTERS[as.numeric(substr(well.temp[grep("R", 
        well.temp)], 2, 3))]
    col.temp <- substr(well.temp[grep("C", well.temp)], 2, 3)
    d$well <- paste0(row.temp, col.temp)
    d$rel.time <- 0
    for (w in unique(d$well)) {
        d[d$well == w, "rel.time"] <- difftime(d[d$well == w, 
            "acq.time"], d[d$well == w, "acq.time"][1], units = "hours")
    }
    a <- aggregate(cell.count ~ well + rel.time, data = d, FUN = sum)
    a <- a[order(a$well), ]
    a$nl2 <- log2norm(a$cell.count, a$well)
    map <- getMapInfo(paste0(dataFol, mapName))
    a$cellLine <- map[match(a$well, map$Well), "Description"]
    a$drug <- "PLX4720"
    a$conc <- map[match(a$well, map$Well), grep("PLX", colnames(map))]
    a
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
