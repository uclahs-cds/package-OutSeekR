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
    data.median <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'median'
        );
    data.trimmean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean',
        trim = 0.05
        );
    data.kmeans <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'kmeans',
        nstart = 1000,
        future.seed = TRUE
        );
    # Transpose so the rows and columns in the computed quantities
    # match those in the original data.
    data.mean <- t(data.mean);
    data.median <- t(data.median);
    data.trimmean <- t(data.trimmean);
    data.kmeans <- t(data.kmeans);
    # Compute the ranges of the z-score statistics.
    data.zrange.mean <- future.apply::future_apply(
        X = data.mean,
        MARGIN = 1,
        FUN = outlier.detection.zrange
        );
    data.zrange.median <- future.apply::future_apply(
        X = data.median,
        MARGIN = 1,
        FUN = outlier.detection.zrange
        );
    data.zrange.trimmean <- future.apply::future_apply(
        X = data.trimmean,
        MARGIN = 1,
        FUN = outlier.detection.zrange
        );
    # Compute the k-means fraction.
    data.fraction.kmeans <- future.apply::future_apply(
        X = data.kmeans,
        MARGIN = 1,
        FUN = outlier.detection.kmeans
        );
    # Transpose so the rows and columns in the computed quantities
    # match those in the original data.
    data.zrange.mean <- t(data.zrange.mean);
    data.zrange.median <- t(data.zrange.median);
    data.zrange.trimmean <- t(data.zrange.trimmean);
    data.fraction.kmeans <- t(data.fraction.kmeans);
    # Compute the cosine similarity.
    data.cosine <- future.apply::future_sapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            outlier.detection.cosine(
                x = as.numeric(data[i, ]),
                distribution = optimal.distribution.data[i]
                );
            }
        );
    names(data.cosine) <- rownames(data);
    # Assemble the statistics from the five methods into a single
    # matrix.
    observed.5method <- cbind(
        zrange.mean = data.zrange.mean[, 'zrange'],
        zrange.median = data.zrange.median[, 'zrange'],
        zrange.trimmean = data.zrange.trimmean[, 'zrange'],
        fraction.kmeans = data.fraction.kmeans[, 'fraction'],
        cosine = data.cosine
        );

    list(
        optimal.distribution.data = optimal.distribution.data,
        data.mean = data.mean,
        data.median = data.median,
        data.trimmean = data.trimmean,
        data.kmeans = data.kmeans,
        data.zrange.mean = data.zrange.mean,
        data.zrange.median = data.zrange.median,
        data.zrange.trimmean = data.zrange.trimmean,
        data.fraction.kmeans = data.fraction.kmeans,
        data.cosine = data.cosine,
        observed.5method = observed.5method
        );
    }
