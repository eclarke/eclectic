% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eclectic.R
\name{counts_matrix}
\alias{counts_matrix}
\title{Creates a counts matrix from an agglomerated data frame.}
\usage{
counts_matrix(agg, additional.columns)
}
\arguments{
\item{agg}{agglomerated data frame}

\item{additional.columns}{any columns to include as columns in the matrix}
}
\value{
default: a numeric matrix with otu as rows and samples as columns. If
  additional columns specified, a data.frame with otus as columns and samples
  as rows
}
\description{
Optionally, metadata columns can be included (useful for some functions), in
which case the output is a data.frame, not a matrix.
}
