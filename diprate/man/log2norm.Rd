\name{log2norm}
\alias{log2norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
log2norm(count, ids, norm_type = c("idx", "ref")[1], norm_id = "0", norm_vals, norm_idx = 1, zero = log2(0.999))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count}{
%%     ~~Describe \code{count} here~~
}
  \item{ids}{
%%     ~~Describe \code{ids} here~~
}
  \item{norm_type}{
%%     ~~Describe \code{norm_type} here~~
}
  \item{norm_id}{
%%     ~~Describe \code{norm_id} here~~
}
  \item{norm_vals}{
%%     ~~Describe \code{norm_vals} here~~
}
  \item{norm_idx}{
%%     ~~Describe \code{norm_idx} here~~
}
  \item{zero}{
%%     ~~Describe \code{zero} here~~
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
function (count, ids, norm_type = c("idx", "ref")[1], norm_id = "0", 
    norm_vals, norm_idx = 1, zero = log2(0.999)) 
{
    l2 <- log2(count)
    l2[is.infinite(l2)] <- zero
    norm <- numeric()
    group <- as.character(unique(ids))
    if (norm_type == "idx") {
        for (i in group) {
            d <- l2[ids == i]
            norm <- append(norm, d - d[norm_idx])
        }
        norm
    }
    else {
        for (i in group) {
            norm_pos <- match(norm_id, as.character(norm_vals)[ids == 
                i])
            d <- l2[ids == i]
            norm <- append(norm, d - d[norm_pos])
        }
        norm
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
