\name{makeMissingTaskArgs}
\alias{makeMissingTaskArgs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeMissingTaskArgs(misdat, fi, defargs = list(), root_dir = "/mnt/quaranta/VDIPRR/HTS006")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{misdat}{
%%     ~~Describe \code{misdat} here~~
}
  \item{fi}{
%%     ~~Describe \code{fi} here~~
}
  \item{defargs}{
%%     ~~Describe \code{defargs} here~~
}
  \item{root_dir}{
%%     ~~Describe \code{root_dir} here~~
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
function (misdat, fi, defargs = list(), root_dir = "/mnt/quaranta/VDIPRR/HTS006") 
{
    ch2_ip <- "None"
    gp <- "/mnt/quaranta/VDIPRR/Control_images/MXgain_new.tif"
    rp <- "False"
    wats <- "True"
    for (n in names(defargs)[names(defargs) \%in\% c("ch2_ip", 
        "rp", "wats")]) {
        if (n == "ch2_ip") 
            ch2_ip <- defargs[[n]]
        if (n == "gp") 
            gp <- defargs[[n]]
        if (n == "rp") 
            rp <- defargs[[n]]
        if (n == "wats") 
            wats <- defargs[[n]]
    }
    fns <- fi$file_name
    if (!all(misdat$file_name \%in\% fns)) {
        message("Not all missing data files found in fileInfo file. Quitting")
        return(NA)
    }
    nip <- paste(root_dir, "/Plate", as.integer(misdat$plate_id), 
        "/", as.character(misdat$file_name), sep = "")
    pid <- as.integer(misdat$plate_id)
    sp <- paste0(root_dir, "/Segmentation")
    well <- as.character(misdat$well)
    out <- data.frame(ch2_im_path = as.character(rep(ch2_ip, 
        length(nip))), gain_path = as.character(gp), nuc_im_path = as.character(nip), 
        plate_id = as.integer(pid), regprops = as.character(rp), 
        save_path = as.character(sp), well = as.character(well), 
        ws = as.character(wats))
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
