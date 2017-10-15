\name{getParam}
\alias{getParam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getParam(drmod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drmod}{
%%     ~~Describe \code{drmod} here~~
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
function (drmod) 
{
    if (class(drmod) != "drc") {
        message("getParam() requires drm object")
        return(invisible(NA))
    }
    yname <- drmod$dataList$names$orName
    xname <- drmod$dataList$names$dName
    ci <- as.data.frame(confint(drmod))
    ci <- cbind(est = coef(drmod), ci)
    rownames(ci) <- c("slope", "Emax", "E0", "EC50")
    colnames(ci) <- c("est", "lower", "upper")
    rrEmax <- ci["Emax", ]/ci["E0", "est"]
    names(rrEmax) <- c("est", "lower", "upper")
    e_halfmax <- predict(drmod, data.frame(drug1.conc = ci["EC50", 
        "est"]), interval = "confidence")
    names(e_halfmax) <- c("est", "lower", "upper")
    rr_e_halfmax <- (ci["Emax", ]/ci["E0", "est"]) + ((ci["E0", 
        "est"] - ci["Emax", "est"])/(2 * ci["E0", "est"]))
    names(rr_e_halfmax) <- c("est", "lower", "upper")
    ic10 <- ED(drmod, ci["E0", "est"] * 0.9, type = "absolute", 
        interval = "delta", display = FALSE)
    rownames(ic10) <- "IC10"
    ic50 <- ED(drmod, ci["E0", "est"]/2, type = "absolute", interval = "delta", 
        display = FALSE)
    rownames(ic50) <- "IC50"
    ic100 <- ED(drmod, 0, type = "absolute", interval = "delta", 
        display = FALSE)
    rownames(ic100) <- "IC100"
    aa <- getAAr(drmod)
    dat <- drmod$origData
    emaxobs <- min(dat[dat[, xname] >= max(dat[, xname])/11, 
        yname])
    rremaxobs <- ifelse(yname == "dip", min(dat[dat[, xname] >= 
        max(dat[, xname])/11, "norm.dip"]), NA)
    out <- rbind(ci, rrEmax, e_halfmax, rr_e_halfmax, ic10[c(1, 
        3, 4)], ic50[c(1, 3, 4)], ic100[c(1, 3, 4)], emaxobs, 
        rremaxobs, aa)
    rownames(out) <- c(rownames(out)[1:4], "rrEmax", "E50", "rrE50", 
        "IC10", "IC50", "IC100", "Emaxobs", "rrEmaxobs", "AAr")
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
