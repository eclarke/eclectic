% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eclectic.R
\name{na_zeros}
\alias{na_zeros}
\title{Set zeros to NA}
\usage{
na_zeros(p)
}
\arguments{
\item{p}{a ggplot object with a `fill' aesthetic}
}
\description{
Set zeros in the data mapped to the fill aesthetic to NA, allowing them to
be given a specific color in the heatmap independent of the fill scale.
}
