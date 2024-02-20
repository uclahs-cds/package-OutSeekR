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
