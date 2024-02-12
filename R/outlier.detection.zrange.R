#' Add range of z-scores to data
#'
#' Compute the range of a vector of z-scores, and add it to the end of the vector.
#'
#' @param x A numeric vector
#' @return `x`, with the range of `x` appended to the end.  The result retains `names(x)`, and the added element is named `'zrange'`.
#'
#' @noRd
outlier.detection.zrange <- function(x) {
    zrange <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE);
    x.with.zrange <- c(x, zrange);
    names(x.with.zrange) <- c(names(x), 'zrange');
    x.with.zrange;
    }
