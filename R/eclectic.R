#' Default BLAST outfmt 6 column names
#'
#' @export
blast6_colnames <- c(
  "qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend",
  "sstart","send","evalue","bitscore")


#' Reads an output file from BLAST in the outfmt=6 format (tab delimited)
#'
#' @param fp file path to blast output
#' @param col.names the column names of the BLAST output. The "std" columns are
#' defined in `blast6_colnames` (default)
#' @export
read_blast6 <- function(fp, col.names=blast6_colnames) {
  infile <- read.delim(fp, col.names=col.names, stringsAsFactors = FALSE)
  return(infile)
}


#' Create a named vector using the values from one column in the second.
#'
#' The levels of `a` are used to pick the corresponding values in `b` so that
#' the named vector returned can be used in e.g. ggplot2::scale label
#' parameters.
#'
#' @param .df the dataframe to reference the columns from
#' @param a the column whose levels will be used as names
#' @param b the values of the named vector
#' @return a named vector
#' @export
named_vector <- function(.df, a, b) {
  la <- levels(as.factor(.df[[a]]))
  lb <- .df[[b]][match(la, .df[[a]])]
  setNames(as.character(lb), la)
}


#' Subset a (counts) matrix based on the column from a dataframe.
#'
#' This is useful for subsetting a counts matrix based on the sample IDs present
#' in the sample data dataframe. It then prunes the matrix to remove any rows
#' whose sums are 0.
#'
#' @param s a sample data dataframe
#' @param mat a counts matrix
#' @param colname the name of the column to use for subsetting
#' @return the subsetted matrix
#' @export
subset_matrix <- function(s, mat, colname="SampleID") {
  .mat <- mat[, s[[colname]]]
  .mat[rowSums(.mat) > 0, ]
}

#' Saturated rainbow for proportions
#'
#' Adapted from qiimer::saturated_rainbow, this provides a ggplot2
#' scale_fill_gradientn that adopts the conventions from that scale. Namely:
#' The lowest-abundance taxa are in dark blue, empty taxa are white, and
#' highly-abundant taxa (> threshold) are in red.
#'
#' @param ... arguments passed to scale_fill_gradientn
#' @param threshold the limit where the high-range color appears
#' @param na.value the color to assign to missing taxa
#' @param scale the the maximum of the proportion (i.e. is it 0-1 scaled or 0-100?)
#' @import ggplot2
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


#' Saturated rainbow for counts
#'
#' Adapted from qiimer::saturated_rainbow, this provides a ggplot2
#' scale_fill_gradientn that adopts the conventions from that scale. Namely:
#' The lowest-abundance taxa are in dark blue, empty taxa are white, and
#' highly-abundant taxa (> threshold) are in red.
#'
#' @param ... arguments passed to scale_fill_gradientn
#' @param na.value the color to assign to missing taxa
#' @param threshold the proportional limit where the high-range color appears
#' @param scale the the maximum of the proportion (i.e. is it 0-1 scaled or 0-100?)
#' @import ggplot2
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

#' Uses the number of columns and number of rows in a heatmap to fix the
#' aspect ratio such that the cells are perfectly square.
#'
#' @param p a ggplot object, ideally with discrete x and y aesthetics set
#' (i.e. a heatmap)
#' @param fudge a fudge factor multiplying the aspect ratio (< 1 = wider, > 1 = taller)
#'
#' @return the same object, with the aspect ratio fixed to be the (number of
#' rows)/(number of columns)
#' @import ggplot2
#' @export
make_square <- function(p, fudge=1) {
  .x <- as.character(p$mapping$x)
  .y <- as.character(p$mapping$y)
  ncols <- length(unique(p$data[[.x]]))
  nrows <- length(unique(p$data[[.y]]))
  print(nrows)
  p + ggplot2::theme(aspect.ratio = (nrows/ncols)*fudge)
}

#' Set zeros to NA
#'
#' Set zeros in the data mapped to the fill aesthetic to NA, allowing them to
#' be given a specific color in the heatmap independent of the fill scale.
#'
#' @param p a ggplot object with a `fill' aesthetic
#' @export
na_zeros <- function(p) {
  .fill <- as.character(p$mapping$fill)
  fill <- p$data[[.fill]]
  fill[fill==0] <- NA
  p$data[[.fill]] <- fill
  p
}
