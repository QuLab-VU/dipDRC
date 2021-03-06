% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dipQC.r
\name{getGCargs}
\alias{getGCargs}
\title{Get growth curve arguments}
\usage{
getGCargs(
  dat,
  arg.name = c("time", "cell.count", "ids"),
  dat.col = c("time", "cell.count", "uid")
)
}
\arguments{
\item{dat}{data.frame containing colnames in \emph{dat.col}}

\item{arg.name}{character of names for output \code{list} of arguments; default is
c('time','cell.count','ids')}

\item{dat.col}{character of colnames present in \emph{dat}; default is
c('time','cell.count','uid')}
}
\value{
list of arguments with names from \emph{arg.name}
}
\description{
Converts data from \code{data.frame} into a \code{list} of arguments that can be
 passed to other functions.
}
