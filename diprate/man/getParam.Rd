% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dipDRC.r
\name{getParam}
\alias{getParam}
\title{Extract parameters from dose-response models (drm objects)}
\usage{
getParam(drmod)
}
\arguments{
\item{drmod}{\emph{drm} object}
}
\value{
data.frame of extracted parameters and metrics (with confidence intervals)
}
\description{
Function to extract relevant parameter values and other extracted metrics from
 dose-response model (drm) objects generated by the \code{drc} package used by
 \code{diprate}.
}
