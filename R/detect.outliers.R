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
    # Determine which of the normal, log-normal, exponential, or gamma
    # distributions provides the best fit to each row of values in
    # `null.data`.
    optimal.distribution.null.data <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = identify.bic.optimal.data.distribution,
        future.seed = TRUE
        );

    # Compute quantities for outlier detection on the null data: (1)
    # z-scores based on the mean / standard deviation, (2) z-scores
    # based on the trimmed mean / trimmed standard deviation, (3)
    # z-scores based on the median / median absolute deviation, and
    # (4) the cluster assignment from k-means with two clusters.
    data.mean <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean'
        );
    data.median <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'median'
        );
    data.trimmean <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean',
        trim = 0.05
        );
    data.kmeans <- future.apply::future_apply(
        X = null.data,
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
        X = seq_len(nrow(null.data)),
        FUN = function(i) {
            outlier.detection.cosine(
                x = as.numeric(null.data[i, ]),
                distribution = optimal.distribution.null.data[i]
                );
            }
        );
    names(cosine.similarity) <- rownames(null.data);
    # Assemble the statistics from the five methods into a single
    # matrix.
    null.5method <- cbind(
        zrange.mean = zrange.mean,
        zrange.median = zrange.median,
        zrange.trimmean = zrange.trimmean,
        fraction.kmeans = fraction.kmeans,
        cosine.similarity = cosine.similarity
        );

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
        null.data = null.data,
        optimal.distribution.null.data = optimal.distribution.null.data,
        null.5method = null.5method
        );
    }

detect.outliers2 <- function(data, num.null) {
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
    observed.residuals <- future.apply::future_apply(
        X = observed.residuals,
        MARGIN = 1,
        FUN = trim.sample
        );
    observed.residuals <- t(observed.residuals);
    # Determine which of the normal, log-normal, exponential, or gamma
    # distributions provides the best fit to each row of residuals.
    optimal.distribution.residuals <- future.apply::future_apply(
        X = observed.residuals,
        MARGIN = 1,
        FUN = identify.bic.optimal.residuals.distribution
        );

    # Compute quantities for outlier detection: (1) z-scores based on
    # the mean / standard deviation, (2) z-scores based on the trimmed
    # mean / trimmed standard deviation, (3) z-scores based on the
    # median / median absolute deviation, and (4) the cluster
    # assignment from k-means with two clusters.
    observed.zrange.mean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean'
        );
    observed.zrange.median <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'median'
        );
    observed.zrange.trimmean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean',
        trim = 0.05
        );
    observed.fraction.kmeans <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'kmeans',
        nstart = 1000,
        future.seed = TRUE
        );
    # Compute the ranges of the z-score statistics.
    observed.zrange.mean <- future.apply::future_apply(
        X = observed.zrange.mean,
        MARGIN = 2,
        FUN = zrange
        );
    observed.zrange.median <- future.apply::future_apply(
        X = observed.zrange.median,
        MARGIN = 2,
        FUN = zrange
        );
    observed.zrange.trimmean <- future.apply::future_apply(
        X = observed.zrange.trimmean,
        MARGIN = 2,
        FUN = zrange
        );
    # Compute the k-means fraction.
    observed.fraction.kmeans <- future.apply::future_apply(
        X = observed.fraction.kmeans,
        MARGIN = 2,
        FUN = kmeans.fraction
        );
    # Compute the cosine similarity.
    observed.cosine.similarity <- future.apply::future_sapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            outlier.detection.cosine(
                x = as.numeric(data[i, ]),
                distribution = optimal.distribution.data[i]
                );
            }
        );
    names(observed.cosine.similarity) <- rownames(data);
    # Assign ranks within each method.
    rank.observed.zrange.mean <- outlier.rank2(
        outlier.statistic = observed.zrange.mean,
        method = 'zrange.mean'
        );
    rank.observed.zrange.median <- outlier.rank2(
        outlier.statistic = observed.zrange.median,
        method = 'zrange.median'
        );
    rank.observed.zrange.trimmean <- outlier.rank2(
        outlier.statistic = observed.zrange.trimmean,
        method = 'zrange.trimmean'
        );
    rank.observed.fraction.kmeans <- outlier.rank2(
        outlier.statistic = observed.fraction.kmeans,
        method = 'fraction.kmeans'
        );
    rank.observed.cosine.similarity <- outlier.rank2(
        outlier.statistic = observed.cosine.similarity,
        method = 'cosine.similarity'
        );
    # Compute the rank product for each transcript.
    rank.product.observed <- outlier.rank.product2(
        zrange.mean = rank.observed.zrange.mean,
        zrange.median = rank.observed.zrange.median,
        zrange.trimmean = rank.observed.zrange.trimmean,
        fraction.kmeans = rank.observed.fraction.kmeans,
        cosine.similarity = rank.observed.cosine.similarity
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
                r = as.numeric(observed.residuals[i, ]),
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
    # Determine which of the normal, log-normal, exponential, or gamma
    # distributions provides the best fit to each row of values in
    # `null.data`.
    optimal.distribution.null.data <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = identify.bic.optimal.data.distribution,
        future.seed = TRUE
        );

    # Compute quantities for outlier detection on the null data: (1)
    # z-scores based on the mean / standard deviation, (2) z-scores
    # based on the trimmed mean / trimmed standard deviation, (3)
    # z-scores based on the median / median absolute deviation, and
    # (4) the cluster assignment from k-means with two clusters.
    null.zrange.mean <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean'
        );
    null.zrange.median <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'median'
        );
    null.zrange.trimmean <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean',
        trim = 0.05
        );
    null.fraction.kmeans <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'kmeans',
        nstart = 1000,
        future.seed = TRUE
        );
    # Compute the ranges of the z-score statistics.
    null.zrange.mean <- future.apply::future_apply(
        X = null.zrange.mean,
        MARGIN = 2,
        FUN = zrange
        );
    null.zrange.median <- future.apply::future_apply(
        X = null.zrange.median,
        MARGIN = 2,
        FUN = zrange
        );
    null.zrange.trimmean <- future.apply::future_apply(
        X = null.zrange.trimmean,
        MARGIN = 2,
        FUN = zrange
        );
    # Compute the k-means fraction.
    null.fraction.kmeans <- future.apply::future_apply(
        X = null.fraction.kmeans,
        MARGIN = 2,
        FUN = kmeans.fraction
        );
    # Compute the cosine similarity.
    null.cosine.similarity <- future.apply::future_sapply(
        X = seq_len(nrow(null.data)),
        FUN = function(i) {
            outlier.detection.cosine(
                x = as.numeric(null.data[i, ]),
                distribution = optimal.distribution.null.data[i]
                );
            }
        );
    names(null.cosine.similarity) <- rownames(null.data);

    list(
        optimal.distribution.data = optimal.distribution.data,
        optimal.distribution.residuals = optimal.distribution.residuals,
        observed.zrange.mean = observed.zrange.mean,
        observed.zrange.median = observed.zrange.median,
        observed.zrange.trimmean = observed.zrange.trimmean,
        observed.fraction.kmeans = observed.fraction.kmeans,
        observed.cosine.similarity = observed.cosine.similarity,
        rank.product.observed = rank.product.observed,
        null.data = null.data,
        optimal.distribution.null.data = optimal.distribution.null.data,
        null.zrange.mean = null.zrange.mean,
        null.zrange.median = null.zrange.median,
        null.zrange.trimmean = null.zrange.trimmean,
        null.fraction.kmeans = null.fraction.kmeans,
        null.cosine.similarity = null.cosine.similarity
        );
    }

detect.outliers3 <- function(
    data,
    num.null,
    p.value.threshold = 0.05,
    kmeans.nstart = 1
    ) {
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
    observed.residuals <- future.apply::future_apply(
        X = observed.residuals,
        MARGIN = 1,
        FUN = trim.sample
        );
    observed.residuals <- t(observed.residuals);
    # Determine which of the normal, log-normal, exponential, or gamma
    # distributions provides the best fit to each row of residuals.
    optimal.distribution.residuals <- future.apply::future_apply(
        X = observed.residuals,
        MARGIN = 1,
        FUN = identify.bic.optimal.residuals.distribution
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
                r = as.numeric(observed.residuals[i, ]),
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
    # Determine which of the normal, log-normal, exponential, or gamma
    # distributions provides the best fit to each row of values in
    # `null.data`.
    optimal.distribution.null.data <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = identify.bic.optimal.data.distribution,
        future.seed = TRUE
        );

    # Calculate p-values.
    p.values <- lapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            calculate.p.values(
                x = as.numeric(data[i, ]),
                x.distribution = optimal.distribution.data[i],
                null.data = null.data,
                null.distributions = optimal.distribution.null.data,
                p.value.threshold = p.value.threshold,
                kmeans.nstart = kmeans.nstart
                );
            }
        );
    p.values <- do.call(
        what = rbind,
        args = p.values
        );
    rownames(p.values) <- rownames(data);

    list(
        optimal.distribution.data = optimal.distribution.data,
        optimal.distribution.residuals = optimal.distribution.residuals,
        null.data = null.data,
        optimal.distribution.null.data = optimal.distribution.null.data,
        p.values
        );
    }
