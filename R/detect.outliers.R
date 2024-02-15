#' Detect outliers
#'
#' @param data A matrix or data frame of FPKM values, organized with transcripts on rows and samples on columns.  Transcript identifiers should be stored as `rownames(data)`.
#'
#' @export
detect.outliers <- function(data) {
    # Determine which of the normal, log-normal, exponential, or gamma
    # distributions provides the best fit to each row of values in
    # `data`.
    optimal.distribution.data <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = identify.bic.optimal.data.distribution,
        future.seed = TRUE
        );
    # Compute quantities for outlier detection: (1) z-scores based on
    # the mean / standard deviation, (2) z-scores based on the trimmed
    # mean / trimmed standard deviation, (3) z-scores based on the
    # median / median absolute deviation, and (4) the cluster
    # assignment from k-means with two clusters.
    data.mean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean'
        );
    data.trimmean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean',
        trim = 0.05
        );
    data.median <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'median'
        );
    data.kmeans <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'kmeans',
        nstart = 1000,
        future.seed = TRUE
        );

    list(
        optimal.distribution.data = optimal.distribution.data,
        data.mean = data.mean,
        data.trimmean = data.trimmean,
        data.median = data.median,
        data.kmeans = data.kmeans
        );
    }
