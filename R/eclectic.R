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


#' Creates a chunk for a plot with a given figure height and width.
#' @param p the plot to insert into the chunk
#' @param width the width of the device
#' @param height the height of the device
#' @param name the name of the chunk
#' @param extras a string to add to the end of the chunk definition, for kicks
#' @export
dynamic_chunk <- function(p, width, height, name=NULL, extras=NULL) {
  chunk_template = "```{r {{name}}, fig.height={{height}}, fig.width={{width}}, echo=FALSE, {{extras}}}
  p
  ```"
  chunk_text = knitr::knit(
    text=knitr::knit_expand(text=chunk_template, width=width,
                     height=height, name=name, extras=extras),
    quiet=TRUE)
  chunk_text
}


#' Agglomerates the sample metadata, counts, and taxa info into one melted df.
#'
#' @param s appropriately-subsetted sample metadata df
#' @param cts counts matrix (colnames should be present in column 'SampleID' in 's')
#' @param taxa taxonomic data (must have column named 'otu')
#' @return the agglomerated df, invisibly
#' @import dplyr
#' @importFrom tidyr gather
#' @export
agglomerate <- function(s, cts, taxa) {
  .c <- cts[, colnames(cts) %in% s$SampleID]
  .c <- as.data.frame(.c[rowSums(.c) > 0, ])
  .c$otu <- rownames(.c)

  agg <- gather(.c, SampleID, count, -otu) %>%
    mutate_each("as.factor", SampleID, otu) %>%
    filter(count > 0) %>%
    merge(s, all.x=TRUE) %>% merge(taxa, all.x=TRUE) %>%
    mutate_each("as.factor", Kingdom:Species)

  invisible(agg)
}

#' Find the mode of a logical or numeric vector
#' @param x a logical or numeric vector
#' @return the mode
#' @export
Mode <- function(x) {
  if (all(is.na(x))) return(NA)
  ux <- unique(na.omit(x))
  ux[which.max(tabulate(match(x, ux)))]
}



