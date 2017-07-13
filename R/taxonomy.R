#' @title Reorder levels of subtaxa by parent taxa
#' @description
#'    Reorders the levels of each taxonomic rank to be arranged by the rank
#'    above it, so that e.g. each class level is ordered by the phylum it
#'    belongs to, and so on.
#' @param
#'    taxonomy the taxonomy data frame, with columns for each [kpcofgs] rank
#' @return
#'    the taxonomy with each column converted to a factor, and the levels
#'    reordered
#' @export
reorder_taxa <- function(taxonomy, ranks=qiimer::taxonomic_ranks) {
  stopifnot(length(ranks) > 1)
  ftax <- dplyr::mutate_each(taxonomy, 'as.factor')
  for (rank in 2:length(ranks)) {
    prev.rank <- ranks[[rank-1]]
    this.rank <- ranks[[rank]]

    ftax[[this.rank]] <- reorder(ftax[[this.rank]], as.numeric(ftax[[prev.rank]]))
  }
  ftax
}

#' Finds the lowest taxonomic rank that isn't NA, stopping at `end`.
#' @param otus the otus to look for
#' @param taxonomy the taxa annotations, with otus as rownames (important!)
#' @param end the lowest desired taxonomic rank
#' @param label label the lowest rank with [kpcofgs]?
#' @param sep separator to use between rank label and rank value
#' @param ranks a character vector of taxonomic ranks, in order
#' @return the lowest non-NA taxonomic rank assignment
#' @export
tax_climber <- function(otus, taxonomy, end="Genus", label=FALSE, sep=":", ranks=qiimer::taxonomic_ranks) {
  end_idx <- which(ranks==end)
  if (is.na(end_idx))
    stop("End rank not found in ranks")
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

#' Creates a MinRank vector from an agglomerated data frame.
#'
#' Version of \link{tax_climber} for use inside \link{dplyr::mutate} or directly
#' assign to col. MinRank columns are the most specific taxonomic assignments
#' (up to a threshold)
#' @param agg \link{agglomerate}d data frame (with taxonomic ranks)
#' @param end.rank the most specific taxonomic rank if available
#' @param rank a character vector of taxonomic ranks present in agg
#' @param ... additional arguments passed to \link{tax_climber}
#' @export
min_rank <- function(agg, end.rank, ranks=qiimer::taxonomic_ranks, ...) {
  .agg <- filter(agg, !is.na(otu))
  md <- .agg[,c(ranks, "otu")]
  md <- data.frame(distinct(md))
  rownames(md) <- md$otu
  minrank <- tax_climber(md$otu, md, end=end.rank, ranks=ranks, ...)
  left_join(select(agg, otu), data.frame(otu=md$otu, MinRank=minrank))$MinRank
}
