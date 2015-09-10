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
