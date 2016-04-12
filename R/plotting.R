#' @title Saturated rainbow for proportions
#' @description
#'   Adapted from qiimer::saturated_rainbow, this provides a ggplot2
#'   ale_fill_gradientn that adopts the conventions from that scale. Namely:
#'   The lowest-abundance taxa are in dark blue, empty taxa are white, and
#'   highly-abundant taxa (> threshold) are in red.
#' @param ... arguments passed to scale_fill_gradientn
#' @param threshold the limit where the high-range color appears
#' @param na.value the color to assign to missing taxa
#' @param scale the the maximum of the proportion (i.e. is it 0-1 scaled or 0-100?)
#' @export
saturated_rainbow_pct <- function(..., na.value="white", threshold=.4, scale=1) {
  rainbow_colors <- rev(rainbow(100*threshold, start=0, end=0.6))
  last_color <- tail(rainbow_colors, n=1)
  colors <- c(rainbow_colors, rep(last_color, 100*(1-threshold)))
  # colors[1] <- "white"
  ggplot2::scale_fill_gradientn(
    na.value=na.value,
    colors = colors,
    limits = c(0,1)*scale,
    breaks = seq(0.1, 0.9, by=0.2)*scale,
    labels = paste0(seq(10, 90, by=20), "%"),
    ...)
}


#' @title Saturated rainbow for counts
#' @description
#'   Adapted from qiimer::saturated_rainbow, this provides a ggplot2
#'   scale_fill_gradientn that adopts the conventions from that scale. Namely:
#'   The lowest-abundance taxa are in dark blue, empty taxa are white, and
#'   highly-abundant taxa (> threshold) are in red.
#' @param ... arguments passed to scale_fill_gradientn
#' @param na.value the color to assign to missing taxa
#' @param threshold the proportional limit where the high-range color appears
#' @param scale the the maximum of the proportion (i.e. is it 0-1 scaled or 0-100?)
#' @export
saturated_rainbow_cts <- function(..., na.value="white", threshold=.4, scale=1) {
  rainbow_colors <- rev(rainbow(100*threshold, start=0, end=0.6))
  last_color <- tail(rainbow_colors, n=1)
  colors <- c(rainbow_colors, rep(last_color, 100*(1-threshold)))
  # colors[1] <- "white"
  ggplot2::scale_fill_gradientn(
    na.value=na.value,
    colors=colors,
    # limits = c(0,1)*scale,
    # breaks = seq(0.1, 0.9, by=0.2)*scale,
    # labels = paste0(seq(10, 90, by=20), "%"),
    ...)
}

#' Find the dimensions of a ggplot2 heatmap.
#'
#' This works by identifying the number of unique points in the variables mapped
#' to the x and y axes.
#'
#' @param p a ggplot2 heatmap (or any ggplot with x and y mappings)
#' @return a list with the number of rows and columns
#' @export
heatmap_dims <- function(p) {
  .x <- as.character(p$mapping$x)
  .y <- as.character(p$mapping$y)
  ncols <- length(unique(p$data[[.x]]))
  nrows <- length(unique(p$data[[.y]]))
  return(list(ncols=ncols, nrows=nrows))
}

#' Uses the number of columns and number of rows in a heatmap to fix the
#' aspect ratio such that the cells are perfectly square.
#'
#' @param p a ggplot object with discrete x and y aesthetics set
#' @param fudge a fudge factor multiplying the aspect ratio (< 1 = wider,
#'   > 1 = taller). This does not appear to work right now.
#' @return the same object, with the aspect ratio fixed to be the (number of
#' rows)/(number of columns)
#' @export
make_square <- function(p, fudge=1) {
  dims <- heatmap_dims(p)
  p + ggplot2::theme(aspect.ratio = (dims$nrows/dims$ncols)*fudge)
}
