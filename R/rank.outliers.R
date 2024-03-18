#' Rank outlier statistics
#'
#' Given a matrix of outlier statistics, organized with transcripts on rows and outlier statistics on columns, return the matrix formed by replacing the values of the outlier statistics with their ranks within each column.
#'
#' @details Ranks are computed such that smaller ranks correspond to more extreme values of the outlier statistics.  For z-ranges based on the mean and standard deviation, trimmed mean and trimmed standard deviation, and median and median absolute deviation, this requires first multiplying the values of the statistics by -1.
#'
#' @param outlier.statistics.matrix A matrix of outlier statistics such as those returned by `quantify.outliers()` and `outlier.detection.cosine()`.
#' @return A numeric matrix with the same dimensions as the input.  Row and column names are preserved.
#'
#' @examples
#' # Generate example outlier statistics.
#' x <- matrix(
#'     data = c(
#'          4.9254,  4.4737,  6.0652,  5.9465,  6.0137,  5.3294,
#'          6.6455,  4.9609, 11.6272, 18.6988,  4.7265,  5.7300,
#'          6.6802,  5.6204, 10.3305, 10.8605,  7.1425,  6.9261,
#'          0.1203,  0.4165,  0.1854,  0.0731,  0.4135,  0.2434,
#'          0.9587,  0.9997,  0.9988,  0.9934,  0.9986,  0.9863
#'         ),
#'     nrow = 6,
#'     ncol = 5,
#'     dimnames = list(
#'         c('A', 'B', 'C', 'D', 'E', 'F'),
#'         c(
#'             'zrange.mean', 'zrange.median', 'zrange.trimmean',
#'             'fraction.kmeans', 'cosine.similarity'
#'             )
#'         )
#'     );
#' outlier.rank(x);
#'
#' @noRd
outlier.rank <- function(outlier.statistics.matrix) {
    methods <- c(
        'zrange.mean', 'zrange.median', 'zrange.trimmean',
        'fraction.kmeans', 'cosine.similarity'
        );
    ranks.matrix <- matrix(
        data = NA,
        nrow = nrow(outlier.statistics.matrix),
        ncol = ncol(outlier.statistics.matrix),
        dimnames = list(
            rownames(outlier.statistics.matrix),
            colnames(outlier.statistics.matrix)
            )
        );
    for (method in methods) {
        if (method %in% c('zrange.mean', 'zrange.median', 'zrange.trimmean')) {
            ranks.matrix[, method] <- rank(
                x = -outlier.statistics.matrix[, method],
                na.last = 'keep',
                ties.method = 'max'
                );
            }
        else if (method %in% c('fraction.kmeans', 'cosine.similarity')) {
            ranks.matrix[, method] <- rank(
                x = outlier.statistics.matrix[, method],
                na.last = 'keep',
                ties.method = 'max'
                );
            }
        }
    ranks.matrix;
    }

outlier.rank2 <- function(outlier.statistic, method) {
    if (method %in% c('zrange.mean', 'zrange.median', 'zrange.trimmean')) {
        ranks <- rank(
            x = -outlier.statistic,
            na.last = 'keep',
            ties.method = 'max'
            );
        }
    else if (method %in% c('fraction.kmeans', 'cosine.similarity')) {
        ranks <- rank(
            x = outlier.statistic,
            na.last = 'keep',
            ties.method = 'max'
            );
        }
    ranks;
    }

#' Rank product of outlier statistics
#'
#' Given a matrix of ranks of outlier statistics, organized with transcripts on rows and outlier statistics on columns, return the products of the ranks across each transcript.
#'
#' @details For \eqn{k} nonmissing ranks, the rank product is defined as
#' \deqn{\sqrt[k]{\prod_{i=1}^k Rank_i}}
#'
#' @param ranks.matrix A matrix of ranks of outlier statistics such as that returned by `outlier.rank()`.
#' @param num.allowed.NA The number of allowable missing values in each row of ranks.  If the number of missing values exceeds `num.allowed.NA`, the corresponding rank product will be `NA`.  Default is zero.
#' @return A numeric vector of length equal to `nrow(ranks.matrix)` with names taken from `rownames(ranks.matrix)`.
#'
#' @examples
#' # Generate example outlier statistics.
#' x <- matrix(
#'     data = c(
#'          4.9254,  4.4737,  6.0652,  5.9465,  6.0137,  5.3294,
#'          6.6455,  4.9609, 11.6272, 18.6988,  4.7265,  5.7300,
#'          6.6802,  5.6204, 10.3305, 10.8605,  7.1425,  6.9261,
#'          0.1203,  0.4165,  0.1854,  0.0731,  0.4135,  0.2434,
#'          0.9587,  0.9997,  0.9988,  0.9934,  0.9986,  0.9863
#'         ),
#'     nrow = 6,
#'     ncol = 5,
#'     dimnames = list(
#'         c('A', 'B', 'C', 'D', 'E', 'F'),
#'         c(
#'             'zrange.mean', 'zrange.median', 'zrange.trimmean',
#'             'fraction.kmeans', 'cosine.similarity'
#'             )
#'         )
#'     );
#' ranks <- outlier.rank(x);
#' outlier.rank.product(
#'     ranks.matrix = ranks
#'     );
#'
#' # Allow at most one missing rank.
#' x <- matrix(
#'     data = c(
#'          4.9254,  4.4737,  6.0652,  5.9465,  6.0137,  5.3294,
#'          6.6455,  4.9609, 11.6272, 18.6988,  4.7265,  5.7300,
#'          6.6802,  5.6204, 10.3305, 10.8605,  7.1425,  6.9261,
#'          0.1203,  0.4165,      NA,  0.0731,      NA,      NA,
#'          0.9587,  0.9997,  0.9988,  0.9934,      NA,  0.9863
#'         ),
#'     nrow = 6,
#'     ncol = 5,
#'     dimnames = list(
#'         c('A', 'B', 'C', 'D', 'E', 'F'),
#'         c(
#'             'zrange.mean', 'zrange.median', 'zrange.trimmean',
#'             'fraction.kmeans', 'cosine.similarity'
#'             )
#'         )
#'     );
#' ranks <- outlier.rank(x);
#' outlier.rank.product(
#'     ranks.matrix = ranks,
#'     num.allowed.NA = 1
#'     );
#'
#' @noRd
outlier.rank.product <- function(ranks.matrix, num.allowed.NA = 0) {
    num.NA <- future.apply::future_apply(
        X = ranks.matrix,
        MARGIN = 1,
        FUN = function(x) sum(is.na(x))
        );
    rank.product <- future.apply::future_apply(
        X = ranks.matrix,
        MARGIN = 1,
        FUN = prod,
        na.rm = TRUE
        );
    rank.product[num.NA > num.allowed.NA] <- NA;
    rank.product^(1 / (ncol(ranks.matrix) - num.NA));
    }
