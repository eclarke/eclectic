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
