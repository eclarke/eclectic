#' Prunes a tree to only include specified OTUs.
#' @param tree an object of class "\code{phylo}"
#' @param otus a list of OTUs (or tip names) to keep
#' @return pruned version of \code{tree}
#' @export
prune_tree <- function(tree, otus) {
  to.drop <- tree$tip.label[!(tree$tip.label %in% otus)]
  ape::drop.tip(tree, to.drop)
}

#' Root a tree using a random node as the root.
#' @param tree a phylogenetic tree (class "phylo")
#' @return rooted tree
#' @export
root_tree <- function(tree, otus) {
  # Root the tree using a random node as the root
  .rtree <- ape::root(tree, sample(tree$tip.label, 1), resolve.root = TRUE)
  stopifnot(ape::is.rooted(.rtree))
  return(.rtree)
}

#' Calculates the generalized UniFrac distances for an agglomerated dataframe.
#'
#' This involves a random rooting of the tree and rarefaction, which are
#' non-deterministic. For reproducibility, recommend calling set.seed before this
#' @param agg agglomerated dataframe
#' @param s sample metadata dataframe
#' @param tree a tree (will root using random node)
#' @param alpha the alpha parameter for GUniFrac
#' @return the generalized UniFrac distance matrix for the given \code{alpha}
#' @export
MakeUnifracData <- function(agg, s, tree, rarefy.depth, alpha=0.5) {

  # Root tree
  .rtree <- root_tree(tree, agg$otu)
  # Create sample matrix
  .mat <- counts_matrix(agg)
  # Keep only otus that appear in the tree
  .mat <- .mat[, colnames(.mat) %in% .rtree$tip.label]

  # Rarefy to specified depth
  if (!missing(rarefy.depth)) {
    .mat.rf <- vegan::rrarefy(.mat, sample=rarefy.depth)
  } else {
    .mat.rf <- .mat
  }

  # Calculate generalized UniFrac
  .unifrac <- GUniFrac::GUniFrac(.mat.rf, .rtree, alpha=alpha)$unifracs
  # Return the generalized UniFrac matrix at specified alpha level
  .unifrac[,,paste0("d_",alpha)]
}
