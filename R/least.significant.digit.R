#' Least significant digit
#'
#' Return the least significant digit from a vector of numbers.
#'
#' @param x A numeric vector.
#' @return The least significant digit among the numbers in `x`, with the digit itself replaced by 1.  See Examples.
#'
#' @examples
#' least.significant.digit(
#'     x = c(0.2, 0.03, 0.004, 0.05)
#'     );
#'
#' @noRd
least.significant.digit <- function(x) {
    decimal.number.max <- sapply(
        X = stats::na.omit(x),
        FUN = function(y) {
            decimal.numbers <- sapply(
                X = y,
                FUN = function(z) {
                    nchar(as.character(z)) - nchar(as.integer(z)) - 1
                    }
                );
            decimal.numbers;
            }
        );
    1 / 10^as.numeric(max(decimal.number.max));
    }
