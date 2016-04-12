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
reorder_taxa <- function(taxonomy) {
  ftax <- dplyr::mutate_each(taxonomy, 'as.factor')
  for (rank in 2:7) {
    prev.rank <- qiimer::taxonomic_ranks[[rank-1]]
    this.rank <- qiimer::taxonomic_ranks[[rank]]

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
