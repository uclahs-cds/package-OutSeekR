#' Detect outliers
#'
#' @param data A matrix or data frame of FPKM values, organized with transcripts on rows and samples on columns.  Transcript identifiers should be stored as `rownames(data)`.
#' @param num.null The number of transcripts to generate when simulating from null distributions.
#' @param p.value.threshold The p-value threshold for the outlier test; default is 0.05.  Once the p-value for a sample exceeds `p.value.threshold`, testing for that transcript ceases, and all samples (including the one that initially exceeded `p.value.threshold`) are assigned p-values equal to `p.value.threshold`.
#' @param kmeans.nstart The number of random starts when computing k-means fraction; default is 1.  See `?stats::kmeans` for further details.
#'
#' @return A list consisting of the following entries:
#' * `p.values`: a matrix of p-values for the outlier test run on each sample for each transcript in `data`.
#' * `num.outliers`: a vector giving the number of outliers (specifically, the number of samples for which the outlier test yielded a p-value less than `p.value.threshold`)  for each transcript.
#' * `outlier.statistics.matrix.list`: a list of length `max(num.outliers) + 1` containing entries `outlier.statistics.matrix.N`, where `N` is between zero and `max(num.outliers)`.  `outlier.statistics.matrix.N` is the matrix of outlier statistics after excluding the `N`th outlier sample, with `outlier.statistics.matrix.0` being for the full data set.  A transcript will only appear in the matrices up to the total number of outliers in that transcript; all transcripts appear in `outlier.statistics.matrix.0`.
#' * `distributions`: a numeric vector indicating the optimal distribution for each transcript.  Possible values are 1 (normal), 2 (log-normal), 3 (exponential), and 4 (gamma).
#'
#' @export
detect.outliers <- function(
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

    # Calculate p-values.  The result is a list of length equal to
    # `nrow(data)`, with each sublist containing the results of
    # `calculate.p.values()` for a transcript in the observed data.
    # See the documentation for `calculate.p.values()` for a
    # description of its return value.
    outlier.test.results <- lapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            x <- as.numeric(data[i, ]);
            names(x) <- colnames(data);
            calculate.p.values(
                x = x,
                x.distribution = optimal.distribution.data[i],
                null.data = null.data,
                null.distributions = optimal.distribution.null.data,
                p.value.threshold = p.value.threshold,
                kmeans.nstart = kmeans.nstart
                );
            }
        );
    #
    # Prepare output
    #
    # Assemble a matrix of p-values for each transcript and sample.
    # Currently, the p-values are stored in the 'p.values' entry of
    # each sublist of `outlier.test.results`, with each sublist
    # corresponding to a different transcript.  We extract these
    # entries and then `rbind()` them to form a matrix.
    p.values <- do.call(
        what = rbind,
        args = lapply(
            X = outlier.test.results,
            FUN = function(x) {
                getElement(
                    object = x,
                    name = 'p.values'
                    );
                }
            )
        );
    rownames(p.values) <- rownames(data);
    # Get counts of the number of outliers per transcript.
    num.outliers <- apply(
        X = p.values,
        MARGIN = 1,
        FUN = function(x) sum(x < p.value.threshold)
        );
    # Create a list of outlier statistics matrices.  The list will be
    # of length `max(num.outliers) + 1` and will contain entries
    # `outlier.statistics.matrix.N`, where `N` is between zero and
    # `max(num.outliers)`.  `outlier.statistics.matrix.N` is the
    # matrix of outlier statistics after excluding the `N`th outlier
    # sample, with `outlier.statistics.matrix.0` being for the full
    # data set.  A transcript will only appear in the matrices up to
    # the total number of outliers in that transcript; all transcripts
    # appear in `outlier.statistics.matrix.0`.
    outlier.statistics.matrix.list <- list();
    for (i in 0:max(num.outliers)) {
        next.name <- paste0('outlier.statistics.matrix.', i);
        outlier.statistics.matrix.list[[next.name]] <- list();
        for (j in seq_along(outlier.test.results)) {
            outlier.statistics.matrix.list[[next.name]][[j]] <- getElement(
                object = outlier.test.results[[j]],
                name = c(
                    'outlier.statistics.list',
                    paste0('outlier.statistics.', i)
                    )
                );
            }
        outlier.statistics.matrix.list[[next.name]] <- do.call(
            what = rbind,
            args = outlier.statistics.matrix.list[[next.name]]
            );
        rownames(outlier.statistics.matrix.list[[next.name]]) <- names(num.outliers[num.outliers >= i]);
        }

    list(
        p.values = p.values,
        num.outliers = num.outliers,
        outlier.statistics.matrix.list = outlier.statistics.matrix.list,
        distributions = optimal.distribution.data
        );
    }

#' Detect outliers, version 2
#'
#' @param data A matrix or data frame of FPKM values, organized with transcripts on rows and samples on columns.  Transcript identifiers should be stored as `rownames(data)`.
#' @param num.null The number of transcripts to generate when simulating from null distributions.
#' @param p.value.threshold The p-value threshold for the outlier test; default is 0.05.  Once the p-value for a sample exceeds `p.value.threshold`, testing for that transcript ceases, and all samples (including the one that initially exceeded `p.value.threshold`) are assigned p-values equal to `p.value.threshold`.
#' @param kmeans.nstart The number of random starts when computing k-means fraction; default is 1.  See `?stats::kmeans` for further details.
#'
#' @return A list consisting of the following entries:
#' * `p.values`: a matrix of p-values for the outlier test run on each sample for each transcript in `data`.
#' * `num.outliers`: a vector giving the number of outliers (specifically, the number of samples for which the outlier test yielded a p-value less than `p.value.threshold`)  for each transcript.
#' * `outlier.statistics.matrix.list`: a list of length `max(num.outliers) + 1` containing entries `outlier.statistics.matrix.N`, where `N` is between zero and `max(num.outliers)`.  `outlier.statistics.matrix.N` is the matrix of outlier statistics after excluding the `N`th outlier sample, with `outlier.statistics.matrix.0` being for the full data set.  A transcript will only appear in the matrices up to the total number of outliers in that transcript; all transcripts appear in `outlier.statistics.matrix.0`.
#' * `distributions`: a numeric vector indicating the optimal distribution for each transcript.  Possible values are 1 (normal), 2 (log-normal), 3 (exponential), and 4 (gamma).
#'
#' @export
detect.outliers2 <- function(
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

    # Calculate p-values.  The result is a list of length equal to
    # `nrow(data)`, with each sublist containing the results of
    # `calculate.p.values()` for a transcript in the observed data.
    # See the documentation for `calculate.p.values()` for a
    # description of its return value.
    outlier.test.results <- future.apply::future_lapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            x <- as.numeric(data[i, ]);
            names(x) <- colnames(data);
            calculate.p.values(
                x = x,
                x.distribution = optimal.distribution.data[i],
                null.data = null.data,
                null.distributions = optimal.distribution.null.data,
                p.value.threshold = p.value.threshold,
                kmeans.nstart = kmeans.nstart
                );
            },
        future.seed = TRUE
        );
    #
    # Prepare output
    #
    # Assemble a matrix of p-values for each transcript and sample.
    # Currently, the p-values are stored in the 'p.values' entry of
    # each sublist of `outlier.test.results`, with each sublist
    # corresponding to a different transcript.  We extract these
    # entries and then `rbind()` them to form a matrix.
    p.values <- do.call(
        what = rbind,
        args = lapply(
            X = outlier.test.results,
            FUN = function(x) {
                getElement(
                    object = x,
                    name = 'p.values'
                    );
                }
            )
        );
    rownames(p.values) <- rownames(data);
    # Get counts of the number of outliers per transcript.
    num.outliers <- apply(
        X = p.values,
        MARGIN = 1,
        FUN = function(x) sum(x < p.value.threshold)
        );
    # Create a list of outlier statistics matrices.  The list will be
    # of length `max(num.outliers) + 1` and will contain entries
    # `outlier.statistics.matrix.N`, where `N` is between zero and
    # `max(num.outliers)`.  `outlier.statistics.matrix.N` is the
    # matrix of outlier statistics after excluding the `N`th outlier
    # sample, with `outlier.statistics.matrix.0` being for the full
    # data set.  A transcript will only appear in the matrices up to
    # the total number of outliers in that transcript; all transcripts
    # appear in `outlier.statistics.matrix.0`.
    outlier.statistics.matrix.list <- list();
    for (i in 0:max(num.outliers)) {
        next.name <- paste0('outlier.statistics.matrix.', i);
        outlier.statistics.matrix.list[[next.name]] <- list();
        for (j in seq_along(outlier.test.results)) {
            outlier.statistics.matrix.list[[next.name]][[j]] <- getElement(
                object = outlier.test.results[[j]],
                name = c(
                    'outlier.statistics.list',
                    paste0('outlier.statistics.', i)
                    )
                );
            }
        outlier.statistics.matrix.list[[next.name]] <- do.call(
            what = rbind,
            args = outlier.statistics.matrix.list[[next.name]]
            );
        rownames(outlier.statistics.matrix.list[[next.name]]) <- names(num.outliers[num.outliers >= i]);
        }

    list(
        p.values = p.values,
        num.outliers = num.outliers,
        outlier.statistics.matrix.list = outlier.statistics.matrix.list,
        distributions = optimal.distribution.data
        );
    }

#' Detect outliers, version 3
#'
#' @param data A matrix or data frame of FPKM values, organized with transcripts on rows and samples on columns.  Transcript identifiers should be stored as `rownames(data)`.
#' @param num.null The number of transcripts to generate when simulating from null distributions.
#' @param p.value.threshold The p-value threshold for the outlier test; default is 0.05.  Once the p-value for a sample exceeds `p.value.threshold`, testing for that transcript ceases, and all samples (including the one that initially exceeded `p.value.threshold`) are assigned p-values equal to `p.value.threshold`.
#' @param kmeans.nstart The number of random starts when computing k-means fraction; default is 1.  See `?stats::kmeans` for further details.
#'
#' @return A list consisting of the following entries:
#' * `p.values`: a matrix of p-values for the outlier test run on each sample for each transcript in `data`.
#' * `num.outliers`: a vector giving the number of outliers (specifically, the number of samples for which the outlier test yielded a p-value less than `p.value.threshold`)  for each transcript.
#' * `outlier.statistics.matrix.list`: a list of length `max(num.outliers) + 1` containing entries `outlier.statistics.matrix.N`, where `N` is between zero and `max(num.outliers)`.  `outlier.statistics.matrix.N` is the matrix of outlier statistics after excluding the `N`th outlier sample, with `outlier.statistics.matrix.0` being for the full data set.  A transcript will only appear in the matrices up to the total number of outliers in that transcript; all transcripts appear in `outlier.statistics.matrix.0`.
#' * `distributions`: a numeric vector indicating the optimal distribution for each transcript.  Possible values are 1 (normal), 2 (log-normal), 3 (exponential), and 4 (gamma).
#'
#' @export
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

    # Compute quantities for outlier detection: (1) z-scores based on
    # the mean / standard deviation, (2) z-scores based on the trimmed
    # mean / trimmed standard deviation, (3) z-scores based on the
    # median / median absolute deviation, and (4) the cluster
    # assignment from k-means with two clusters.
    data.zrange.mean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean'
        );
    data.zrange.median <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'median'
        );
    data.zrange.trimmean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean',
        trim = 0.05
        );
    data.fraction.kmeans <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'kmeans',
        nstart = kmeans.nstart,
        future.seed = TRUE
        );
    # Compute the ranges of the z-score statistics.
    data.zrange.mean <- future.apply::future_apply(
        X = data.zrange.mean,
        MARGIN = 2,
        FUN = zrange
        );
    data.zrange.median <- future.apply::future_apply(
        X = data.zrange.median,
        MARGIN = 2,
        FUN = zrange
        );
    data.zrange.trimmean <- future.apply::future_apply(
        X = data.zrange.trimmean,
        MARGIN = 2,
        FUN = zrange
        );
    # Compute the k-means fraction.
    data.fraction.kmeans <- future.apply::future_apply(
        X = data.fraction.kmeans,
        MARGIN = 2,
        FUN = kmeans.fraction
        );
    # Compute the cosine similarity.
    data.cosine.similarity <- future.apply::future_sapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            outlier.detection.cosine(
                x = as.numeric(data[i, ]),
                distribution = optimal.distribution.data[i]
                );
            }
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
        nstart = kmeans.nstart,
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

    # Calculate p-values.  The result is a list of length equal to
    # `nrow(data)`, with each sublist containing the results of
    # `calculate.p.values()` for a transcript in the observed data.
    # See the documentation for `calculate.p.values()` for a
    # description of its return value.
    outlier.test.results <- future.apply::future_lapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            x <- as.numeric(data[i, ]);
            names(x) <- colnames(data);
            calculate.p.values2(
                x = x,
                x.distribution = optimal.distribution.data[i],
                x.zrange.mean = data.zrange.mean[i],
                x.zrange.median = data.zrange.median[i],
                x.zrange.trimmean = data.zrange.trimmean[i],
                x.fraction.kmeans = data.fraction.kmeans[i],
                x.cosine.similarity = data.cosine.similarity[i],
                null.data = null.data,
                null.distributions = optimal.distribution.null.data,
                null.zrange.mean = null.zrange.mean,
                null.zrange.median = null.zrange.median,
                null.zrange.trimmean = null.zrange.trimmean,
                null.fraction.kmeans = null.fraction.kmeans,
                null.cosine.similarity = null.cosine.similarity,
                p.value.threshold = p.value.threshold,
                kmeans.nstart = kmeans.nstart
                );
            },
        future.seed = TRUE
        );
    #
    # Prepare output
    #
    # Assemble a matrix of p-values for each transcript and sample.
    # Currently, the p-values are stored in the 'p.values' entry of
    # each sublist of `outlier.test.results`, with each sublist
    # corresponding to a different transcript.  We extract these
    # entries and then `rbind()` them to form a matrix.
    p.values <- do.call(
        what = rbind,
        args = lapply(
            X = outlier.test.results,
            FUN = function(x) {
                getElement(
                    object = x,
                    name = 'p.values'
                    );
                }
            )
        );
    rownames(p.values) <- rownames(data);
    # Get counts of the number of outliers per transcript.
    num.outliers <- apply(
        X = p.values,
        MARGIN = 1,
        FUN = function(x) sum(x < p.value.threshold)
        );
    # Create a list of outlier statistics matrices.  The list will be
    # of length `max(num.outliers) + 1` and will contain entries
    # `outlier.statistics.matrix.N`, where `N` is between zero and
    # `max(num.outliers)`.  `outlier.statistics.matrix.N` is the
    # matrix of outlier statistics after excluding the `N`th outlier
    # sample, with `outlier.statistics.matrix.0` being for the full
    # data set.  A transcript will only appear in the matrices up to
    # the total number of outliers in that transcript; all transcripts
    # appear in `outlier.statistics.matrix.0`.
    outlier.statistics.matrix.list <- list();
    for (i in 0:max(num.outliers)) {
        next.name <- paste0('outlier.statistics.matrix.', i);
        outlier.statistics.matrix.list[[next.name]] <- list();
        for (j in seq_along(outlier.test.results)) {
            outlier.statistics.matrix.list[[next.name]][[j]] <- getElement(
                object = outlier.test.results[[j]],
                name = c(
                    'outlier.statistics.list',
                    paste0('outlier.statistics.', i)
                    )
                );
            }
        outlier.statistics.matrix.list[[next.name]] <- do.call(
            what = rbind,
            args = outlier.statistics.matrix.list[[next.name]]
            );
        rownames(outlier.statistics.matrix.list[[next.name]]) <- names(num.outliers[num.outliers >= i]);
        }

    list(
        p.values = p.values,
        num.outliers = num.outliers,
        outlier.statistics.matrix.list = outlier.statistics.matrix.list,
        distributions = optimal.distribution.data
        );
    }

#' Detect outliers, precompute version
#'
#' @param data A matrix or data frame of FPKM values, organized with transcripts on rows and samples on columns.  Transcript identifiers should be stored as `rownames(data)`.
#' @param num.null The number of transcripts to generate when simulating from null distributions.
#' @param num.null.precompute The number of outlier statistics to precompute.  TODO elaborate
#' @param p.value.threshold The p-value threshold for the outlier test; default is 0.05.  Once the p-value for a sample exceeds `p.value.threshold`, testing for that transcript ceases, and all samples (including the one that initially exceeded `p.value.threshold`) are assigned p-values equal to `p.value.threshold`.
#' @param kmeans.nstart The number of random starts when computing k-means fraction; default is 1.  See `?stats::kmeans` for further details.
#'
#' @return A list consisting of the following entries:
#' * `p.values`: a matrix of p-values for the outlier test run on each sample for each transcript in `data`.
#' * `num.outliers`: a vector giving the number of outliers (specifically, the number of samples for which the outlier test yielded a p-value less than `p.value.threshold`)  for each transcript.
#' * `outlier.statistics.matrix.list`: a list of length `max(num.outliers) + 1` containing entries `outlier.statistics.matrix.N`, where `N` is between zero and `max(num.outliers)`.  `outlier.statistics.matrix.N` is the matrix of outlier statistics after excluding the `N`th outlier sample, with `outlier.statistics.matrix.0` being for the full data set.  A transcript will only appear in the matrices up to the total number of outliers in that transcript; all transcripts appear in `outlier.statistics.matrix.0`.
#' * `distributions`: a numeric vector indicating the optimal distribution for each transcript.  Possible values are 1 (normal), 2 (log-normal), 3 (exponential), and 4 (gamma).
#'
#' @export
detect.outliers.precompute <- function(
    data,
    num.null,
    num.null.precompute = 10,
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

    # Compute quantities for outlier detection: (1) z-scores based on
    # the mean / standard deviation, (2) z-scores based on the trimmed
    # mean / trimmed standard deviation, (3) z-scores based on the
    # median / median absolute deviation, and (4) the cluster
    # assignment from k-means with two clusters.
    data.zrange.mean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean'
        );
    data.zrange.median <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'median'
        );
    data.zrange.trimmean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean',
        trim = 0.05
        );
    data.fraction.kmeans <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'kmeans',
        nstart = kmeans.nstart,
        future.seed = TRUE
        );
    # Compute the ranges of the z-score statistics.
    data.zrange.mean <- future.apply::future_apply(
        X = data.zrange.mean,
        MARGIN = 2,
        FUN = zrange
        );
    data.zrange.median <- future.apply::future_apply(
        X = data.zrange.median,
        MARGIN = 2,
        FUN = zrange
        );
    data.zrange.trimmean <- future.apply::future_apply(
        X = data.zrange.trimmean,
        MARGIN = 2,
        FUN = zrange
        );
    # Compute the k-means fraction.
    data.fraction.kmeans <- future.apply::future_apply(
        X = data.fraction.kmeans,
        MARGIN = 2,
        FUN = kmeans.fraction
        );
    # Compute the cosine similarity.
    data.cosine.similarity <- future.apply::future_sapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            outlier.detection.cosine(
                x = as.numeric(data[i, ]),
                distribution = optimal.distribution.data[i]
                );
            }
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

    # Precompute outlier statistics for the complete null data set and
    # after excluding up to `num.null.precompute` samples (i.e.,
    # columns).
    #
    # Initialize a list to store lists of vectors of outlier
    # statistics.  `null.statistics.list[['0']]` will store the
    # outlier statistics for the complete null data in a named list
    # with entries `null.zrange.mean`, `null.zrange.median`,
    # `null.zrange.trimmean, `null.fraction.kmeans`, and
    # `null.cosine.similarity`.  `null.statistics.list[['1']] will
    # store the outlier statistics for the null data after excluding
    # the first column in a similar named list, and so on.
    null.statistics.list <- list();
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
        nstart = kmeans.nstart,
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
    null.statistics.list[['0']] <- list(
        null.zrange.mean = null.zrange.mean,
        null.zrange.median = null.zrange.median,
        null.zrange.trimmean = null.zrange.trimmean,
        null.fraction.kmeans = null.fraction.kmeans,
        null.cosine.similarity = null.cosine.similarity
        );
    for (i in seq_len(num.null.precompute)) {
        null.zrange.mean <- future.apply::future_apply(
            X = null.data[, -1 * seq_len(i), drop = FALSE],
            MARGIN = 1,
            FUN = quantify.outliers,
            method = 'mean'
            );
        null.zrange.median <- future.apply::future_apply(
            X = null.data[, -1 * seq_len(i), drop = FALSE],
            MARGIN = 1,
            FUN = quantify.outliers,
            method = 'median'
            );
        null.zrange.trimmean <- future.apply::future_apply(
            X = null.data[, -1 * seq_len(i), drop = FALSE],
            MARGIN = 1,
            FUN = quantify.outliers,
            method = 'mean',
            trim = 0.05
            );
        null.fraction.kmeans <- future.apply::future_apply(
            X = null.data[, -1 * seq_len(i), drop = FALSE],
            MARGIN = 1,
            FUN = quantify.outliers,
            method = 'kmeans',
            nstart = kmeans.nstart,
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
                    x = as.numeric(null.data[i, -1 * seq_len(i), drop = FALSE]),
                    distribution = optimal.distribution.null.data[i]
                    );
                }
            );
        null.statistics.list[[as.character(i)]] <- list(
            null.zrange.mean = null.zrange.mean,
            null.zrange.median = null.zrange.median,
            null.zrange.trimmean = null.zrange.trimmean,
            null.fraction.kmeans = null.fraction.kmeans,
            null.cosine.similarity = null.cosine.similarity
            );
        }

    # Calculate p-values.  The result is a list of length equal to
    # `nrow(data)`, with each sublist containing the results of
    # `calculate.p.values()` for a transcript in the observed data.
    # See the documentation for `calculate.p.values()` for a
    # description of its return value.
    outlier.test.results <- future.apply::future_lapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            x <- as.numeric(data[i, ]);
            names(x) <- colnames(data);
            calculate.p.values2(
                x = x,
                x.distribution = optimal.distribution.data[i],
                x.zrange.mean = data.zrange.mean[i],
                x.zrange.median = data.zrange.median[i],
                x.zrange.trimmean = data.zrange.trimmean[i],
                x.fraction.kmeans = data.fraction.kmeans[i],
                x.cosine.similarity = data.cosine.similarity[i],
                null.data = null.data,
                null.distributions = optimal.distribution.null.data,
                null.zrange.mean = null.zrange.mean,
                null.zrange.median = null.zrange.median,
                null.zrange.trimmean = null.zrange.trimmean,
                null.fraction.kmeans = null.fraction.kmeans,
                null.cosine.similarity = null.cosine.similarity,
                p.value.threshold = p.value.threshold,
                kmeans.nstart = kmeans.nstart
                );
            },
        future.seed = TRUE
        );
    #
    # Prepare output
    #
    # Assemble a matrix of p-values for each transcript and sample.
    # Currently, the p-values are stored in the 'p.values' entry of
    # each sublist of `outlier.test.results`, with each sublist
    # corresponding to a different transcript.  We extract these
    # entries and then `rbind()` them to form a matrix.
    p.values <- do.call(
        what = rbind,
        args = lapply(
            X = outlier.test.results,
            FUN = function(x) {
                getElement(
                    object = x,
                    name = 'p.values'
                    );
                }
            )
        );
    rownames(p.values) <- rownames(data);
    # Get counts of the number of outliers per transcript.
    num.outliers <- apply(
        X = p.values,
        MARGIN = 1,
        FUN = function(x) sum(x < p.value.threshold)
        );
    # Create a list of outlier statistics matrices.  The list will be
    # of length `max(num.outliers) + 1` and will contain entries
    # `outlier.statistics.matrix.N`, where `N` is between zero and
    # `max(num.outliers)`.  `outlier.statistics.matrix.N` is the
    # matrix of outlier statistics after excluding the `N`th outlier
    # sample, with `outlier.statistics.matrix.0` being for the full
    # data set.  A transcript will only appear in the matrices up to
    # the total number of outliers in that transcript; all transcripts
    # appear in `outlier.statistics.matrix.0`.
    outlier.statistics.matrix.list <- list();
    for (i in 0:max(num.outliers)) {
        next.name <- paste0('outlier.statistics.matrix.', i);
        outlier.statistics.matrix.list[[next.name]] <- list();
        for (j in seq_along(outlier.test.results)) {
            outlier.statistics.matrix.list[[next.name]][[j]] <- getElement(
                object = outlier.test.results[[j]],
                name = c(
                    'outlier.statistics.list',
                    paste0('outlier.statistics.', i)
                    )
                );
            }
        outlier.statistics.matrix.list[[next.name]] <- do.call(
            what = rbind,
            args = outlier.statistics.matrix.list[[next.name]]
            );
        rownames(outlier.statistics.matrix.list[[next.name]]) <- names(num.outliers[num.outliers >= i]);
        }

    list(
        p.values = p.values,
        num.outliers = num.outliers,
        outlier.statistics.matrix.list = outlier.statistics.matrix.list,
        distributions = optimal.distribution.data
        );
    }
