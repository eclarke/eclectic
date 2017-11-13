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

#' Prepare agglomerated dataframe for heatmap display.
#'
#' Aggregates counts by minimum rank, completing missing cases, and clustering
#' rows. Expects input \code{agg} to have columns called \code{SampleID} and
#' \code{count} at minimum.
#' @param agg an agglomerated df from \link{agglomerate}.
#' @param s sample metadata dataframe
#' @param min.rank.1 the most specific rank to aggregate by (determines y-axis)
#'   Uses \link{min_rank} to select next most specific rank if not available.
#' @param min.rank.2 (optional) higher rank to prefix min.rank.1 (e.g. Phylum)
#'   Uses \link{min_rank} to select next most specific rank if not available.
#' @param min.samples a row must appear more in than this number of samples to
#'   be shown
#' @param hclust_method a character vector specifying the method to use for
#'   clustering y-axis (from \link{hclust})
#' @param ... additional samples passed to the min_rank calls
#' @export
MakeHeatmapData <- function(
  agg, s, min.rank.1 = "Genus", min.rank.2 = NA, min.samples=1,
  hclust_method="average",  ...) {
  # Create the y-axes for the data via the desired taxonomic rank level(s)
  .agg <- agg %>% mutate(
    MinRank1 = eclectic::min_rank(agg, end=min.rank.1, ...))
  if (!is.na(min.rank.2)) {
    .agg <- .agg %>% mutate(
      MinRank2 = eclectic::min_rank(agg, end=min.rank.2, ...),
      MinRank1 = paste(MinRank2, MinRank1)
    )
  }

  .agg <- .agg %>%
    # Aggregate counts by minimum available rank
    group_by(SampleID, MinRank1) %>%
    summarize(count = sum(count)) %>%
    # Update proportions
    group_by(SampleID) %>%
    mutate(proportion=count/sum(count))

  # Convert to matrix form
  .mat <- reshape2::dcast(
    .agg, MinRank1 ~ SampleID, value.var="proportion", fill = 0) %>%
    filter(!is.na(MinRank1))
  rownames(.mat) <- .mat$MinRank1
  .mat$MinRank1 <- NULL

  # Cluster and pull out row order
  minrank_order <- (dist(.mat, method="euclidean") %>%
                      hclust(method=hclust_method))$order

  .agg %>%
    # Reorder by clustering order
    mutate(MinRank1 = factor(MinRank1, levels=rownames(.mat)[minrank_order])) %>%
    # Complete missing cases: fills missing combinations in with NA
    # instead of just omitting the row entirely (necessary to show blank cells)
    ungroup() %>%
    tidyr::complete(SampleID, nesting(MinRank1)) %>%
    select(SampleID, MinRank1, count, proportion) %>%
    distinct() %>%
    # Filter taxa that appear in fewer than required number of samples
    group_by(MinRank1) %>%
    filter(sum(proportion > 0, na.rm=TRUE) > min.samples) %>%
    # Re-add metadata
    left_join(s, by="SampleID", all.y=TRUE) %>%
    ungroup()
}

#' Plot a standard rainbow heatmap with faceting.
#'
#' Uses the output of \link{MakeHeatmapData} to display a heatmap with (by default)
#' SampleIDs on the x-axis and a minimum taxonomic rank of Genus on the y-axis.
#' @param heatmap.data the dataframe from \link{MakeHeatmapData}
#' @param use.reads if TRUE, plot the raw readcounts; if FALSE, use proportions.
#'   Raw readcounts are subject to biases from uneven library size but can be
#'   preferable for low-biomass or extremely sparse samples, as proportions can
#'   be distorting.
#' @param threshold at what relative proportion the color scale goes to red.
#'   This should be adjusted such that there is only one red cell per row.
#' @param x.text display labels on the x-axis?
#' @param x.axis name of factor to be used to define the x-axis (must exist in
#'   \code{heatmap.data})
#' @export
PlotHeatmap <- function(
  heatmap.data, use.reads=TRUE, threshold=0.6, x.text=TRUE, x.axis='SampleID') {

  if (use.reads) {
    fill.var <- "count"
    scale.name <- "Reads"
    fill.scale <- saturated_rainbow_cts(threshold = threshold, name=scale.name)
  } else {
    fill.var <- "proportion"
    scale.name <- "Abundance"
    fill.scale <- saturated_rainbow_pct(threshold = threshold, name=scale.name)
  }

  p <- ggplot(heatmap.data, aes_string(x.axis, "MinRank1", fill=fill.var)) +
    geom_tile(color="grey50", size=0.4) +
    facet_grid(. ~ StudyGroup + SampleType, space="free", scales="free") +
    fill.scale +
    theme_grey() +
    theme(
      strip.text.y = element_text(angle=0, vjust=0),
      strip.text.x = element_text(angle=90, vjust=0),
      panel.border = element_blank())
  if (x.text) {
    p <- p + theme(
      axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  } else{
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  }
  p
}
