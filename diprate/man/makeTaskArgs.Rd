\name{makeTaskArgs}
\alias{makeTaskArgs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeTaskArgs(fi, root_dir, defargs = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fi}{
%%     ~~Describe \code{fi} here~~
}
  \item{root_dir}{
%%     ~~Describe \code{root_dir} here~~
}
  \item{defargs}{
%%     ~~Describe \code{defargs} here~~
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
function (fi, root_dir, defargs = list(), ...) 
{
    rp <- "FALSE"
    ch1_name <- "Cy3"
    ch2_name <- "None"
    expdir <- ""
    sp <- file.path(root_dir, "Segmentation")
    owrite <- "TRUE"
    for (n in names(defargs)[names(defargs) \%in\% c("rp", "regprops", 
        "ch1_name", "ch2_name", "expdir", "gp", "sp", "overwrite", 
        "owrite")]) {
        if (n \%in\% c("rp", "regprops")) 
            rp <- defargs[[n]]
        if (n == "ch1_name") 
            ch1_name <- defargs[[n]]
        if (n == "ch2_name") 
            ch2_name <- defargs[[n]]
        if (n == "expdir") 
            expdir <- defargs[[n]]
        if (n == "sp") 
            sp <- defargs[[n]]
        if (n \%in\% c("owrite", "overwrite")) 
            owrite <- defargs[[n]]
    }
    args <- as.list(substitute(list(...)))[-1L]
    nip <- paste(root_dir, expdir, "/Plate", as.integer(fi[fi$channel == 
        ch1_name, "plate_id"]), "/", as.character(fi[fi$channel == 
        ch1_name, "file_name"]), sep = "")
    if (ch2_name != "None") 
        ch2_ip <- paste(root_dir, expdir, "/Plate", as.integer(fi[fi$channel == 
            ch2_name, "plate_id"]), "/", as.character(fi[fi$channel == 
            ch2_name, "file_name"]), sep = "")
    else ch2_ip <- rep("None", length(nip))
    pid <- as.integer(fi[fi$channel == ch1_name, "plate_id"])
    well <- as.character(fi[fi$channel == ch1_name, "well"])
    out <- data.frame(ch2_im_path = as.character(ch2_ip), nuc_im_path = as.character(nip), 
        overwrite = as.character(owrite), plate_id = as.integer(pid), 
        regprops = as.character(rp), save_path = as.character(sp), 
        well = as.character(well))
    if (length(args) >= 1) {
        message(paste(length(args), "extra args included in output"))
        for (i in seq(length(args))) {
            out <- cbind(out, args[i])
            colnames(out)[ncol(out)] <- names(args)[i]
        }
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
