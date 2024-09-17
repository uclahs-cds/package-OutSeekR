#' Trim a vector of numbers
#'
#' Symmetrically trim a vector of numbers after sorting it.
#'
#' @details If `length(x) <= 10`, the function returns `x[2:(length(x) - 1)]`.
#'
#' @param x A numeric vector.
#' @param trim A number, the fraction of observations to be trimmed from each end of `x`.
#' @return A sorted, trimmed copy of `x`.
#' @export trim.sample
#' @examples
#' trim.sample(
#'     x = 1:20,
#'     trim = 0.05
#'     );
trim.sample <- function(x, trim = 0.05) {
    x <- sort(x);
    if (length(x) <= 10) {
        patient.trim.value <- 2:(length(x) - 1);
        } else {
        trim.sample.number <- length(x) * trim;
        trim.sample.number.integer <- round(trim.sample.number);
        patient.trim.value <- (trim.sample.number.integer + 1):(length(x) - trim.sample.number.integer);
        }
    x[patient.trim.value];
    }
