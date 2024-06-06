#' Range of z-scores
#'
#' Compute the range of a vector of z-scores.
#'
#' @param x A numeric vector
#' @return A number.
#'
#' @examples
#' set.seed(1234);
#' x <- rnorm(
#'     n = 10
#'     );
#' zrange(
#'     x = x
#'     );
#'
#' @noRd
zrange <- function(x) {
    if (all(is.na(x))) {
        return(NA);
        }
    zrange <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE);
    zrange;
    }
