\name{plotMultiCurve}
\alias{plotMultiCurve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotMultiCurve(drm_list, norm = FALSE, leg.scale = 0.75, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drm_list}{
%%     ~~Describe \code{drm_list} here~~
}
  \item{norm}{
%%     ~~Describe \code{norm} here~~
}
  \item{leg.scale}{
%%     ~~Describe \code{leg.scale} here~~
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
function (drm_list, norm = FALSE, leg.scale = 0.75, ...) 
{
    if (class(drm_list) != "list" | !all(lapply(drm_list, class) == 
        "drc")) {
        message("plotMultiCurve needs list of drm objects")
        return(invisible(NA))
    }
    pargs <- as.list(substitute(list(...)))[-1L]
    drm_cl <- unlist(lapply(drm_list, function(x) unique(x$origData$cell.line)))
    drm_dr <- unlist(lapply(drm_list, function(x) unique(x$origData$drug1)[unique(x$origData$drug1) != 
        "control"]))
    if (!"main" \%in\% names(pargs)) {
        pargs["main"] <- ifelse(all(drm_dr == drm_dr[1]), drm_dr[1], 
            "")
    }
    drm_names <- paste(drm_cl, drm_dr)
    if (all(drm_dr == drm_dr[1])) 
        drm_names <- drm_cl
    if (!"ylim" \%in\% names(pargs)) {
        if (norm) {
            pargs <- append(pargs, list(ylim = c(-0.2, 1.2)))
        }
        else {
            pargs <- append(pargs, list(ylim = c(-0.02, 0.07)))
        }
    }
    par(mar = c(4, 4, 2, 0.5))
    n <- length(drm_list)
    mycolors <- c("red", "orange", "yellow", "green", "blue", 
        "purple", "brown", "black")
    if (n <= 8) {
        rep.col <- mycolors[seq(n)]
    }
    else {
        rep.col <- mycolors[c(1:8, (seq(n)[8:n]\%\%8) + 1)]
    }
    rep.lty <- rep(seq(ceiling(n/8)), each = 8)[1:n]
    allargs <- append(list(x = drm_list[[1]], xlim = c(0, 1e-05), 
        log = "x", type = "none", xlab = "[drug], M", ylab = "Response ratio", 
        lwd = 0, col = "white"), pargs)
    do.call(plot, allargs)
    lapply(seq(length(drm_list)), function(x) {
        addLL4curve(drm_list[[x]], col = rep.col[x], lty = rep.lty[x], 
            norm = norm, ...)
        invisible(return())
    })
    legend("bottomleft", legend = drm_names, col = rep.col, lty = rep.lty, 
        lwd = leg.scale, cex = leg.scale)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
