#' @title Grouped color palette for factors
#' @description
#'    Picks a starting hue for each level of the factor, then returns colors
#'    with that hue, but varying luminosity for each item within that factor.
#' @param x the vector or factor grouped by f (i.e. Genus)
#' @param f the vector or factor that groups the values of x (i.e. Order)
#' @param h the range of hues to use
#' @param c chroma value to use
#' @param l the base luminance [0,100]
#' @param l.range the luminance above and below the base to use for subcolors
#' @param h.start the starting hue
#' @param direction 1=clockwise, 0=counterclockwise around color wheel
#' @return
#'    a list containing `major`: the base color for each level in f; and
#'    `minor`: the color for each item in x
#' @export
subcolor_pal <- function (x, f, h = c(0, 360) + 15, c = 100, l=65, l.range=20, h.start = 0, direction = 1, na.value="grey")
{
  stopifnot(length(x) == length(f))
  stopifnot(l-l.range >= 0)
  stopifnot(l+l.range <= 100)

  .f <- droplevels(as.factor(f))
  n = length(levels(.f))
  if ((diff(h)%%360) < 1) {
    h[2] <- h[2] - 360/n
  }
  rotate <- function(x) (x + h.start)%%360 * direction
  hues <- rotate(seq(h[1], h[2], length.out = n))

  levels(.f) <- hues
  major.colors <- hcl(hues, c, l)
  subcolors <- lapply(split(.f, .f), function(subgroup) {
    hue = as.integer(as.character(unique(subgroup)))
    stopifnot(length(hue) == 1)

    sub_n <- length(subgroup)

    if (sub_n > 1)
      l <- seq(l-l.range, l+l.range, length.out = sub_n)

    grDevices::hcl(hue, c, l)
  })

  minor.colors <- unsplit(subcolors, .f)
  minor.colors[is.na(minor.colors)] <- na.value
  names(minor.colors) <- as.character(x)

  levels(.f) <- major.colors
  major.colors <- as.character(.f)
  major.colors[is.na(major.colors)] <- na.value
  names(major.colors) <- as.character(f)

  list(major=major.colors, minor=minor.colors)
}
