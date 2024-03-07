#' Detect outliers
#'
#' @param data A matrix or data frame of FPKM values, organized with transcripts on rows and samples on columns.  Transcript identifiers should be stored as `rownames(data)`.
#' @param num.null The number of transcripts to generate when simulating from null distributions.
#'
#' @export
detect.outliers <- function(data, num.null) {
    # Determine which of the normal, log-normal, exponential, or gamma
    # distributions provides the best fit to each row of values in
    # `data`.
    optimal.distribution.data <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = identify.bic.optimal.data.distribution,
        future.seed = TRUE
        );
    # Calculate residuals of the observed data with respect to the
    # optimal distribution.  (We use `as.numeric()` on the input to
    # `calculate.residuals()` in order to handle data frame input.
    observed.residuals <- future.apply::future_lapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            calculate.residuals(
                x = as.numeric(data[i, ]),
                distribution = optimal.distribution.data[i]
                );
            }
        );
    observed.residuals <- do.call(
        what = rbind,
        args = observed.residuals
        );
    rownames(observed.residuals) <- rownames(data);
    # Apply 5% trimming to each row of residuals.
    observed.residuals.trimmed <- future.apply::future_apply(
        X = observed.residuals,
        MARGIN = 1,
        FUN = trim.sample
        );
    observed.residuals.trimmed <- t(observed.residuals.trimmed);
    # Determine which of the normal, log-normal, exponential, or gamma
    # distributions provides the best fit to each row of residuals.
    optimal.distribution.residuals <- future.apply::future_apply(
        X = observed.residuals.trimmed,
        MARGIN = 1,
        FUN = identify.bic.optimal.residuals.distribution
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
    # Compute the ranges of the z-score statistics.
    zrange.mean <- future.apply::future_apply(
        X = data.mean,
        MARGIN = 2,
        FUN = zrange
        );
    zrange.median <- future.apply::future_apply(
        X = data.median,
        MARGIN = 2,
        FUN = zrange
        );
    zrange.trimmean <- future.apply::future_apply(
        X = data.trimmean,
        MARGIN = 2,
        FUN = zrange
        );
    # Compute the k-means fraction.
    fraction.kmeans <- future.apply::future_apply(
        X = data.kmeans,
        MARGIN = 2,
        FUN = kmeans.fraction
        );
    # Compute the cosine similarity.
    cosine.similarity <- future.apply::future_sapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            outlier.detection.cosine(
                x = as.numeric(data[i, ]),
                distribution = optimal.distribution.data[i]
                );
            }
        );
    names(cosine.similarity) <- rownames(data);
    # Assemble the statistics from the five methods into a single
    # matrix.
    observed.5method <- cbind(
        zrange.mean = zrange.mean,
        zrange.median = zrange.median,
        zrange.trimmean = zrange.trimmean,
        fraction.kmeans = fraction.kmeans,
        cosine.similarity = cosine.similarity
        );
    # Assign ranks within each method.
    observed.5method.ranks <- outlier.rank(
        outlier.statistics.matrix = observed.5method
        );
    # Compute the rank product for each transcript.
    observed.5method.rank.product <- outlier.rank.product(
        ranks.matrix = observed.5method.ranks,
        num.allowed.NA = 0
        );

    # Generate a matrix of null transcripts by simulating from their
    # respective optimal distributions.
    sampled.indices <- sample(
        x = nrow(data),
        size = num.null,
        replace = TRUE
        );
    null.data <- future.apply::future_lapply(
        X = sampled.indices,
        FUN = function(i) {
            simulate.null(
                x = as.numeric(data[i, ]),
                x.distribution = optimal.distribution.data[i],
                r = as.numeric(observed.residuals.trimmed[i, ]),
                r.distribution = optimal.distribution.residuals[i]
                );
            },
        future.seed = TRUE
        );
    null.data <- do.call(
        what = rbind,
        args = null.data
        );
    rownames(null.data) <- rownames(data)[sampled.indices];

    list(
        optimal.distribution.data = optimal.distribution.data,
        optimal.distribution.residuals = optimal.distribution.residuals,
        data.mean = data.mean,
        data.median = data.median,
        data.trimmean = data.trimmean,
        data.kmeans = data.kmeans,
        zrange.mean = zrange.mean,
        zrange.median = zrange.median,
        zrange.trimmean = zrange.trimmean,
        fraction.kmeans = fraction.kmeans,
        cosine.similarity = cosine.similarity,
        observed.5method = observed.5method,
        observed.5method.ranks = observed.5method.ranks,
        observed.5method.rank.product = observed.5method.rank.product,
        null.data = null.data
        );
    }
