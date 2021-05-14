#' Fast entropy estimation of multi-dimensional data.
#'
#' Non-parametric estimator for the differential entropy of a multidimensional distribution, given a limited set of data points, by a recursive rectilinear partitioning.
#'
#' @param X data
#' @param z Z-score threshold
#'
#' @return Differential entropy estimate.
#'
#' @references
#' D. Stowell and M. D. Plumbley
#' Fast multidimensional entropy estimation by k-d partitioning.
#' IEEE Signal Processing Letters 16 (6), 537--540, June 2009.
#' http://dx.doi.org/10.1109/LSP.2009.2017346
#'
#' @examples
#' Xu <- matrix(runif(1000 * 100), ncol=100)
#' kdpee(Xu)
#'
#' Xn <- matrix(rnorm(1000 * 100), ncol=100)
#' kdpee(Xn)
#'
#' @export
kdpee <- function(X, z=1.96) {
  .Call(do_kdpee, X, z)
}
