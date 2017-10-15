\name{prepCountData}
\alias{prepCountData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
prepCountData(dat, time.count.id.names = c("time", "cell.count", "well"), max.cell.count = 3000, cell.line.name = unique(dat$cell.line), ctrl.type = "mean", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{time.count.id.names}{
%%     ~~Describe \code{time.count.id.names} here~~
}
  \item{max.cell.count}{
%%     ~~Describe \code{max.cell.count} here~~
}
  \item{cell.line.name}{
%%     ~~Describe \code{cell.line.name} here~~
}
  \item{ctrl.type}{
%%     ~~Describe \code{ctrl.type} here~~
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
function (dat, time.count.id.names = c("time", "cell.count", 
    "well"), max.cell.count = 3000, cell.line.name = unique(dat$cell.line), 
    ctrl.type = "mean", ...) 
{
    tryCatch({
        ctrl <- dat[dat$drug1 == "control", ]
        dat <- dat[dat$drug1 != "control", ]
        cn <- colnames(dat)
        if (length(cell.line.name) > 1) {
            message("WARNING: Length of prepCountData cell.line.name > 1, using first value only")
            cell.line.name <- cell.line.name[1]
        }
        dat[dat[time.count.id.names[2]] > max.cell.count, time.count.id.names[2]] <- NA
        dat <- dat[!is.na(dat[, time.count.id.names[2]]), ]
        ctrl <- tryCatch({
            controlQC(ctrl[, time.count.id.names[1]], ctrl[, 
                time.count.id.names[2]], ctrl[, time.count.id.names[3]], 
                col.names = time.count.id.names, plotIt = FALSE, 
                ctrl.type = ctrl.type, ...)
        }, error = function(cond) {
            message(paste("error in controlQC parsing data from:", 
                cell.line.name))
            message(cond)
            return(NA)
        })
        if (length(ctrl) != 1) {
            temp <- tail(dat, nrow(ctrl))
            temp[, time.count.id.names[1:2]] <- ctrl
            temp[, !(colnames(temp) \%in\% c(time.count.id.names, 
                "cell.line", "expt.date"))] <- NA
            temp$uid <- paste0(cell.line.name, "_control")
            if ("ucd" \%in\% cn) 
                temp$ucd <- paste0(cell.line.name, "_control")
            if ("plate.name" \%in\% cn & all(dat$plate.name == 
                dat$plate.name[1])) 
                temp$plate.name <- dat$plate.name[1]
            temp$drug1 <- "control"
            temp$drug1.conc <- 0
            temp[, time.count.id.names[3]] <- "CTRL"
            return(rbind(dat, temp))
        }
    }, error = function(cond) {
        message(paste("error in prepCountData processing:", cell.line.name))
        message(cond)
        cat("\n")
        return(NA)
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
