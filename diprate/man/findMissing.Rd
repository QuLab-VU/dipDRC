\name{findMissing}
\alias{findMissing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findMissing(fi, datapath = "", verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fi}{
%%     ~~Describe \code{fi} here~~
}
  \item{datapath}{
%%     ~~Describe \code{datapath} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (fi, datapath = "", verbose = FALSE) 
{
    if (class(fi) != "data.frame") {
        fipath <- ifelse(class(fi) == "character", fi, "")
        if (file.exists(fipath)) {
            fin <- tail(unlist(strsplit(fipath, "/")), 1)
            message(paste("Reading", fin, "as imageFileInfo file"))
            root_dir <- sub(fin, "", fipath)
            fi <- read.csv(fipath, as.is = TRUE)
        }
    }
    if (datapath == "" | !any(grepl("cellcount.csv", dir(datapath)))) {
        message("No data found. Must specify <datapath> containing cellcount.csv files")
        return(NA)
    }
    else if (datapath != "" & dir.exists(datapath) & any(grepl("cellcount.csv", 
        dir(datapath)))) 
        datafns <- dir(datapath)[grepl("cellcount.csv", dir(datapath))]
    if (verbose) 
        message(paste("Found", length(datafns), "data files"))
    plate_ids <- unique(fi$plate_id)
    missing_pids <- setdiff(as.integer(substr(datafns, 6, 10)), 
        plate_ids)
    fixWellName <- function(wellName) sapply(wellName, FUN = function(x) ifelse(nchar(x) == 
        2, paste0(substr(x, 1, 1), 0, substr(x, 2, 2)), x))
    wells <- sapply(2:23, function(col) sapply(LETTERS[2:15], 
        function(row) fixWellName(paste0(row, col))))
    dimnames(wells) <- NULL
    missing_data <- data.frame()
    pb <- txtProgressBar(1, length(datafns), style = 3)
    for (fn in datafns) {
        if (!verbose) 
            setTxtProgressBar(pb, match(fn, datafns))
        pid <- as.integer(substr(fn, 6, 10))
        temp <- read.csv(file.path(datapath, fn), as.is = TRUE)
        if (nrow(temp) < 308) {
            missing_wells <- wells[sapply(wells, function(x) !(x \%in\% 
                temp$well))]
            if (verbose) 
                message(paste(length(missing_wells), "wells without cell count data in ", 
                  tail(unlist(strsplit(fn, "/")), 1)))
            for (w in missing_wells) {
                im_fn <- fi[fi$plate_id == pid & fi$well == w, 
                  "file_name"]
                missing_data <- rbind(missing_data, data.frame(file_name = im_fn, 
                  plate_id = pid, well = w))
            }
        }
    }
    close(pb)
    message(paste("Found", nrow(missing_data), "images without cell counts"))
    invisible(missing_data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
