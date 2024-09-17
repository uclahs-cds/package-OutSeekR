#' k-means fraction
#'
#' Given a vector of cluster assigments from `quantify.outliers()` run with `method = 'kmeans'`, compute the fraction of observations belonging to the smaller of the two clusters.
#'
#' @details This function only considers clusters 1 and 2 even if `quantify.outliers()` was run with `exclude.zero = TRUE`.  In that case, zeros are effectively excluded from the counts used to define the k-means fraction.  See examples.
#'
#' @param x A numeric vector.
#' @return A number.
#' @export
#' @examples
#' x <- c(1, 1, 2, 2, 2, 2, 2, 2, 2, 2);
#' names(x) <- letters[1:length(x)];
#' kmeans.fraction(x);
kmeans.fraction <- function(x) {
    if (1 == length(unique(x))) {
        fraction <- NA;
        }
    else {
        n.cluster.1 <- length(x[1 == x]);
        n.cluster.2 <- length(x[2 == x]);
        total.n.cluster <- n.cluster.1 + n.cluster.2;
        n.smaller.cluster <- min(n.cluster.1, n.cluster.2);
        fraction <- round(
            x = n.smaller.cluster / total.n.cluster,
            digits = 4
            );
        }
    fraction;
    }
