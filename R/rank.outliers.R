#' Rank outlier statistics
#'
#' Given a vector of outlier statistics and the name of the method used to compute them, return the vector of ranks of the statistics.
#'
#' @details Ranks are computed such that smaller ranks correspond to more extreme values of the outlier statistics.  For z-ranges based on the mean and standard deviation, trimmed mean and trimmed standard deviation, and median and median absolute deviation, this requires first multiplying the values of the statistics by -1.
#'
#' @param outlier.statistic A numeric vector of outlier statistics such as that returned by `quantify.outliers()` and `outlier.detection.cosine()`.
#' @param method A string describing the method by which the values in `outlier.statistic` were computed.  One of `'zrange.mean'`, `'zrange.median'`, `'zrange.trimmean'`, `'fraction.kmeans'`, or `'cosine.similarity'`.
#'
#' @return A numeric vector of ranks.  Names are preserved.
#'
#' @examples
#' # Generate example outlier statistics.
#' outlier.stats.matrix <- matrix(
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
#' # Assign ranks within each method.
#' ranks.zrange.mean <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'zrange.mean'],
#'     method = 'zrange.mean'
#'     );
#' ranks.zrange.median <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'zrange.median'],
#'     method = 'zrange.median'
#'     );
#' ranks.zrange.trimmean <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'zrange.trimmean'],
#'     method = 'zrange.trimmean'
#'     );
#' ranks.fraction.kmeans <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'fraction.kmeans'],
#'     method = 'fraction.kmeans'
#'     );
#' ranks.cosine.similarity <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'cosine.similarity'],
#'     method = 'cosine.similarity'
#'     );
#'
#' @noRd
outlier.rank <- function(outlier.statistic, method) {
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
#' Given ranks for each of the five outlier statistics, return the rank product.
#'
#' @details For \eqn{k} nonmissing ranks, the rank product is defined as
#' \deqn{\sqrt[k]{\prod_{i=1}^k Rank_i}}
#'
#' @param zrange.mean A numeric vector of ranks from the range of z-scores based on the mean and standard deviation.
#' @param zrange.median A numeric vector of ranks from the range of z-scores based on the median and median absolute deviation.
#' @param zrange.trimmean A numeric vector of ranks from the range of z-scores based on the trimmed mean and trimmed standard deviation.
#' @param fraction.kmeans A numeric vector of ranks from the cluster assignment from k-means with two clusters.
#' @param cosine.similarity A numeric vector of ranks from the cosine similarity.
#' @param num.allowed.NA The number of allowable missing values in each row of ranks.  If the number of missing values exceeds `num.allowed.NA`, the corresponding rank product will be `NA`.  Default is zero.
#'
#' @return A numeric vector of length equal to `nrow(ranks.matrix)` with names taken from `rownames(ranks.matrix)`.
#'
#' @examples
#'
#' # Generate example outlier statistics.
#' outlier.stats.matrix <- matrix(
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
#' # Assign ranks within each method.
#' ranks.zrange.mean <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'zrange.mean'],
#'     method = 'zrange.mean'
#'     );
#' ranks.zrange.median <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'zrange.median'],
#'     method = 'zrange.median'
#'     );
#' ranks.zrange.trimmean <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'zrange.trimmean'],
#'     method = 'zrange.trimmean'
#'     );
#' ranks.fraction.kmeans <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'fraction.kmeans'],
#'     method = 'fraction.kmeans'
#'     );
#' ranks.cosine.similarity <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'cosine.similarity'],
#'     method = 'cosine.similarity'
#'     );
#'
#' outlier.rank.product(
#'     zrange.mean = ranks.zrange.mean,
#'     zrange.median = ranks.zrange.median,
#'     zrange.trimmean = ranks.zrange.trimmean,
#'     fraction.kmeans = ranks.fraction.kmeans,
#'     cosine.similarity = ranks.cosine.similarity,
#'     num.allowed.NA = 0
#'     );
#'
#' # Allow at most one missing rank.
#' outlier.stats.matrix <- matrix(
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
#'
#' ranks.zrange.mean <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'zrange.mean'],
#'     method = 'zrange.mean'
#'     );
#' ranks.zrange.median <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'zrange.median'],
#'     method = 'zrange.median'
#'     );
#' ranks.zrange.trimmean <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'zrange.trimmean'],
#'     method = 'zrange.trimmean'
#'     );
#' ranks.fraction.kmeans <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'fraction.kmeans'],
#'     method = 'fraction.kmeans'
#'     );
#' ranks.cosine.similarity <- outlier.rank(
#'     outlier.statistic = outlier.stats.matrix[, 'cosine.similarity'],
#'     method = 'cosine.similarity'
#'     );
#'
#' outlier.rank.product(
#'     zrange.mean = ranks.zrange.mean,
#'     zrange.median = ranks.zrange.median,
#'     zrange.trimmean = ranks.zrange.trimmean,
#'     fraction.kmeans = ranks.fraction.kmeans,
#'     cosine.similarity = ranks.cosine.similarity,
#'     num.allowed.NA = 1
#'     );
#'
#' @noRd
outlier.rank.product <- function(
    zrange.mean,
    zrange.median,
    zrange.trimmean,
    fraction.kmeans,
    cosine.similarity,
    num.allowed.NA = 0
    ) {
    num.NA <- future.apply::future_mapply(
        FUN = function(x1, x2, x3, x4, x5) {
            sum(c(is.na(x1), is.na(x2), is.na(x3), is.na(x4), is.na(x5)))
            },
        zrange.mean,
        zrange.median,
        zrange.trimmean,
        fraction.kmeans,
        cosine.similarity
        );
    rank.product <- future.apply::future_mapply(
        FUN = prod,
        zrange.mean,
        zrange.median,
        zrange.trimmean,
        fraction.kmeans,
        cosine.similarity,
        MoreArgs = list(na.rm = TRUE)
        );
    rank.product[num.NA > num.allowed.NA] <- NA;
    rank.product^(1 / (5 - num.NA));
    }
