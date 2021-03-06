% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps.R
\name{saturated_rainbow_cts}
\alias{saturated_rainbow_cts}
\title{Saturated rainbow for counts}
\usage{
saturated_rainbow_cts(..., na.value = "white", threshold = 0.4, scale = 1)
}
\arguments{
\item{...}{arguments passed to scale_fill_gradientn}

\item{na.value}{the color to assign to missing taxa}

\item{threshold}{the proportional limit where the high-range color appears}

\item{scale}{the the maximum of the proportion (i.e. is it 0-1 scaled or 0-100?)}
}
\description{
Adapted from qiimer::saturated_rainbow, this provides a ggplot2
  scale_fill_gradientn that adopts the conventions from that scale. Namely:
  The lowest-abundance taxa are in dark blue, empty taxa are white, and
  highly-abundant taxa (> threshold) are in red.
}
