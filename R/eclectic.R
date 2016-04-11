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


#' Uses a Dirichlet multinomial distribution to test for OTUs in the first
#' vector that are enriched.
#'
#' @param counts1 the first paired sample (and the target of the one-sided test)
#' @param counts2 the second paired sample (i.e. prewash)
#' @param p.adjust.method the method to adjust p.values for multiple testing
#' @return a vector of adjusted p values
#' @export
polyatest <- function(counts1, counts2, p.adjust.method="fdr") {
  if (!requireNamespace("polyafit")) {
    stop("Please install the polyafit package from https://github.com/kylebittinger/polyafit")
  }
  if (length(counts1) == 1 | length(counts2) == 1) {
    warning("Only one observation; returning NA for both conditions")
    return(NaN)
  }
  mat <- matrix(c(as.integer(counts1), as.integer(counts2)), nrow=2, byrow=TRUE)
  # mat <- as.integer(mat)
  fit <- polyafit::optim_polya(mat)
  p.values <- 1 - polyafit::ppolya_marginal(mat[1,], fit$par, log.p=FALSE)
  p.adjust(p.values, method=p.adjust.method)
}


#' Test for an OTU's (inverse) correlation between its abundance and amplicon concentration.
#'
#' @param .df the data frame
#' @param ... the grouping variables (unquoted), i.e. ExtractionType and SampleType
#' @param otu the column with unique otu names
#' @param freq the column denoting OTU proportional abundance (as string)
#' @param conc the column containing amplicon concentration data
#' @return the data frame with correlation results added in for each OTU
#' @import dplyr
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @export
contam_test <- function(.df, ..., otu="otu", freq="freq", conc="PostPCRConc") {
  .df %>% group_by_(.dots=lazyeval::lazy_dots(...)) %>%
    group_by_(.dots = otu, add=TRUE) %>%
    mutate(otuCount = n()) %>% ungroup() %>% filter(otuCount > 1) %>%
    group_by_(.dots=lazyeval::lazy_dots(...)) %>% group_by_(.dots=otu, add=TRUE) %>%
    nest() %>%
    mutate(model = purrr::map(data, ~cor.test(x=.[[freq]], y=.[[conc]], alt="less", method="spearman"))) %>%
    unnest(model %>% purrr::map(broom::glance)) %>%
    mutate(adj.p.value = p.adjust(p.value, method="fdr")) %>%
    unnest(data)
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

#' Finds the lowest taxonomic rank that isn't NA, stopping at `end`.
#' @param otus the otus to look for
#' @param taxonomy the taxa annotations, with otus as rownames (important!)
#' @param end the lowest desired taxonomic rank
#' @param label label the lowest rank with [kpcofgs]?
#' @param sep separator to use between rank label and rank value
#' @return the lowest non-NA taxonomic rank assignment
#' @export
tax_climber <- function(otus, taxonomy, end="Genus", label=FALSE, sep=":") {
  end_idx <- which(qiimer::taxonomic_ranks==end)
  otus <- as.character(otus)
  if (all(is.na(otus)))
    return(rep(NA, length(otus)))
  taxa <- taxonomy[otus, 1:end_idx]
  min.ranks <- colnames(taxa)[apply(taxa, 1, function(x) max(which(!is.na(x))))]
  lowest <- taxa[cbind(otus, min.ranks)]
  if (label)
    paste(tolower(substr(min.ranks, 1, 1)), lowest, sep=sep)
  else
    lowest
}
